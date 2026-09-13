"""
Exon-level table — one row per (UniProt isoform, coding-transcript exon).

For every isoform in the MAIN resolved set (mapping/resolve_isoform_transcripts.py:
resolved_unique + resolved_collapsed_equivalent) we take its CHOSEN/RESOLVED ENST
and emit one row per exon of that transcript, annotated with genomic + CDS +
amino-acid coordinates, IDR content, gene-wide constitutive/inclusion statistics,
and principal-transcript status (MANE / Ensembl-canonical / APPRIS).

Alternative-exon definition used by the paper: `exon_usage`, copied from s9
(analysis/exon_idr_annotation.csv, keyed by isoform + exon_id). Among a gene's
mapped transcripts (genes with >=2), an exon present in all of them is
"constitutive", in some but not all "alternative"; genes with one mapped transcript
are "single_tx_gene"; isoforms outside the s9 table are null.

`is_constitutive` and `transcript_inclusion_fraction` are kept for reference only:
they are computed over **all GENCODE protein-coding transcripts of the gene**, under
which ~92% of coding exons count as non-constitutive. Requires a fresh GTF pass
over the relevant genes.

Inputs (all already produced by the pipeline):
    isoform_mapping/isoform_transcript_mapping_resolved.parquet   (02)  chosen ENST + status
    interim/ensembl_transcript_exon_coordinates.parquet          (s2)  per-exon genomic/CDS/phase
    interim/transcript_cds_to_protein_coordinate_map.parquet     (s2)  aa coords per coding exon
    isoform_mapping/exon_region_class.parquet                    (02)  idr_frac + region_class
    analysis/exon_idr_annotation.csv                             (s9)  exon_usage (paper definition)
    processed/tf_idr_isoform_master.parquet                      (s5)  scaffold flag, canonical
    cache/gencode.v{REL}.annotation.gtf.gz                       (s1)  gene-wide tx set + tags

Outputs:
    processed/exon_level_table.parquet
    processed/exon_level_table.csv
    processed/exon_level_table_DATA_DICTIONARY.md
"""
import sys, gzip, re
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C

OUT_PARQUET = C.PROCESSED / "exon_level_table.parquet"
OUT_CSV     = C.PROCESSED / "exon_level_table.csv"
OUT_DICT    = C.PROCESSED / "exon_level_table_DATA_DICTIONARY.md"

MAIN = {"resolved_unique", "resolved_collapsed_equivalent"}

_GID   = re.compile(r'gene_id "([^"]+)"')
_TID   = re.compile(r'transcript_id "([^"]+)"')
_TTYPE = re.compile(r'transcript_type "([^"]+)"')
_EID   = re.compile(r'exon_id "([^"]+)"')
_TAG   = re.compile(r'tag "([^"]+)"')


def parse_gtf(target_ensg):
    """Gene-wide transcript catalog for the target genes, from the GENCODE GTF.

    Returns:
      tx_meta        {ENST: {gene, ttype, is_mane, is_canonical, appris_level}}
      gene_coding_tx {ENSG: set(protein_coding ENST)}
      exon_tx        {(ENSG, ENSE): set(ENST containing that exon)}
    """
    tx_meta, gene_coding_tx, exon_tx = {}, {}, {}
    with gzip.open(C.GENCODE_GTF_GZ, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            feat = f[2]
            if feat not in ("transcript", "exon"):
                continue
            attrs = f[8]
            g = _GID.search(attrs)
            if not g:
                continue
            gene = g.group(1).split(".")[0]
            if gene not in target_ensg:
                continue
            t = _TID.search(attrs)
            enst = t.group(1).split(".")[0] if t else None
            if enst is None:
                continue
            if feat == "transcript":
                tt = _TTYPE.search(attrs)
                ttype = tt.group(1) if tt else ""
                tags = set(_TAG.findall(attrs))
                # APPRIS principal level = trailing digits of an appris_principal_* tag
                # (robust to non-numeric suffixes; None if not principal)
                appris_lvls = [int(m.group(1)) for x in tags
                               if x.startswith("appris_principal")
                               for m in [re.search(r"(\d+)$", x)] if m]
                lvl = min(appris_lvls, default=None)
                tx_meta[enst] = {
                    "gene": gene, "ttype": ttype,
                    "is_mane": "MANE_Select" in tags,
                    "is_canonical": "Ensembl_canonical" in tags,
                    "appris_level": lvl,
                }
                if ttype == "protein_coding":
                    gene_coding_tx.setdefault(gene, set()).add(enst)
            else:  # exon
                e = _EID.search(attrs)
                if not e:
                    continue
                exon_id = e.group(1).split(".")[0]
                exon_tx.setdefault((gene, exon_id), set()).add(enst)
    return tx_meta, gene_coding_tx, exon_tx


def main():
    print("[load] resolved mapping (MAIN set) ...")
    res = pd.read_parquet(C.ISOFORM_MAPPING / "isoform_transcript_mapping_resolved.parquet")
    res = res[res["resolved_mapping_status"].isin(MAIN) & res["ENST_resolved"].notna()].copy()
    res["ENST"] = res["ENST_resolved"].astype(str).str.split(".").str[0]
    # unversion ENSG/ENSP so they key against the unversioned GTF ids below
    res["ENSG"] = res["ENSG"].astype(str).str.split(".").str[0]
    res["ENSP_resolved"] = res["ENSP_resolved"].astype(str).str.split(".").str[0]
    print(f"       {len(res):,} MAIN isoforms with a chosen ENST")

    target_ensg = set(res["ENSG"].dropna())

    # per-exon genomic/CDS/phase for the chosen transcripts (from s2)
    ex = pd.read_parquet(C.INTERIM / "ensembl_transcript_exon_coordinates.parquet")
    chosen = set(res["ENST"])
    ex = ex[ex["ensembl_transcript_id"].isin(chosen)].copy()

    # aa coordinates per coding exon (from s2 CDS->protein map)
    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")[
        ["ensembl_transcript_id", "exon_rank", "aa_start", "aa_end", "aa_length",
         "cds_start_bp_in_cds", "frame"]]
    ex = ex.merge(cds, on=["ensembl_transcript_id", "exon_rank"], how="left")

    # attach the isoform(s) that resolved to each ENST (one row per isoform x exon)
    keep_res = ["isoform_accession", "base_accession", "gene_name", "is_tf",
                "tf_family", "ENSG", "ENST", "ENSP_resolved", "n_candidate_ENST",
                "resolved_mapping_status"]
    df = res[keep_res].merge(ex, left_on="ENST", right_on="ensembl_transcript_id", how="inner")
    print(f"       {len(df):,} (isoform x exon) rows across {df['ENST'].nunique():,} transcripts")

    # IDR content per exon (from 02 exon_region_class), keyed by (isoform, exon_id)
    erc = pd.read_parquet(C.ISOFORM_MAPPING / "exon_region_class.parquet")[
        ["iso", "exon_id", "idr_frac", "region_class", "is_idr_exon"]].drop_duplicates(["iso", "exon_id"])
    df = df.merge(erc, left_on=["isoform_accession", "exon_id"],
                  right_on=["iso", "exon_id"], how="left").drop(columns=["iso"])

    # scaffold validation + canonical flag (from master)
    m = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")[
        ["isoform_accession", "is_canonical", "genomic_scaffold_validated"]].drop_duplicates("isoform_accession")
    df = df.merge(m, on="isoform_accession", how="left")

    # alternative / constitutive exon usage among mapped transcripts (from s9; paper definition)
    usage = pd.read_csv(C.ANALYSIS / "exon_idr_annotation.csv", usecols=["iso", "exon_id", "usage"])
    df = df.merge(usage.rename(columns={"iso": "isoform_accession", "usage": "exon_usage"}),
                  on=["isoform_accession", "exon_id"], how="left")

    # ── gene-wide constitutive / inclusion + principal tags (fresh GTF pass) ──
    print(f"[gtf] scanning {len(target_ensg):,} genes for coding-transcript sets + tags ...")
    tx_meta, gene_coding_tx, exon_tx = parse_gtf(target_ensg)
    n_coding_gene = {g: len(s) for g, s in gene_coding_tx.items()}

    def inclusion(row):
        gene, exon_id = row["ENSG"], row["exon_id"]
        coding = gene_coding_tx.get(gene, set())
        n_gene = len(coding)
        if n_gene == 0:
            return pd.Series([np.nan, np.nan, np.nan, False])
        with_exon = exon_tx.get((gene, exon_id), set()) & coding
        n_with = len(with_exon)
        frac = n_with / n_gene
        return pd.Series([n_gene, n_with, frac, n_with == n_gene])

    df[["n_coding_transcripts_gene", "n_coding_transcripts_with_exon",
        "transcript_inclusion_fraction", "is_constitutive"]] = df.apply(inclusion, axis=1)

    # principal-transcript status of the CHOSEN ENST
    df["is_mane_select"]      = df["ENST"].map(lambda e: tx_meta.get(e, {}).get("is_mane", False))
    df["is_ensembl_canonical"] = df["ENST"].map(lambda e: tx_meta.get(e, {}).get("is_canonical", False))
    df["appris_principal_level"] = df["ENST"].map(lambda e: tx_meta.get(e, {}).get("appris_level"))
    df["is_appris_principal"]  = df["appris_principal_level"].notna()
    df["is_principal"] = df["is_mane_select"] | df["is_appris_principal"]

    # convenience flags
    df["is_first_exon"] = df["exon_rank"] == 1
    maxrank = df.groupby("ENST")["exon_rank"].transform("max")
    df["is_last_exon"] = df["exon_rank"] == maxrank

    # ── tidy column order (requested schema first, extras after) ──────────────
    df = df.rename(columns={
        "ENSG": "gene_id", "ENST": "transcript_id", "ENSP_resolved": "protein_id",
        "chromosome": "chr", "exon_genomic_start": "start", "exon_genomic_end": "end",
        "exon_rank": "exon_number", "cds_genomic_start": "cds_start",
        "cds_genomic_end": "cds_end",
    })
    df["chosen_enst"] = df["transcript_id"]      # explicit alias
    df["mapping_status"] = df["resolved_mapping_status"]

    cols = [
        # identity
        "gene_id", "gene_name", "transcript_id", "protein_id", "exon_id",
        "isoform_accession", "base_accession", "is_tf", "tf_family",
        "mapping_status", "chosen_enst", "n_candidate_ENST",
        # coordinates
        "chr", "start", "end", "strand", "exon_number",
        "cds_start", "cds_end", "cds_phase", "aa_start", "aa_end",
        "is_coding_exon", "exon_length_bp", "cds_length_bp", "aa_length",
        # exon usage: mapped transcripts (paper definition), then GENCODE-wide (reference)
        "exon_usage",
        "is_constitutive", "transcript_inclusion_fraction",
        "n_coding_transcripts_with_exon", "n_coding_transcripts_gene",
        # principal status
        "is_principal", "is_mane_select", "is_ensembl_canonical",
        "is_appris_principal", "appris_principal_level",
        # IDR content
        "idr_frac", "region_class", "is_idr_exon",
        # isoform-level context
        "is_canonical", "genomic_scaffold_validated",
        "is_first_exon", "is_last_exon",
    ]
    cols = [c for c in cols if c in df.columns]
    df = df[cols].sort_values(["gene_name", "isoform_accession", "exon_number"]).reset_index(drop=True)
    df = C.stamp(df)

    df.to_parquet(OUT_PARQUET, index=False)
    df.to_csv(OUT_CSV, index=False)
    write_dictionary(df)

    # ── report ────────────────────────────────────────────────────────────────
    print(f"\n[out] exon_level_table  {len(df):,} rows x {df.shape[1]} cols")
    print(f"      isoforms {df['isoform_accession'].nunique():,} | transcripts {df['transcript_id'].nunique():,} "
          f"| genes {df['gene_id'].nunique():,}")
    cod = df[df["is_coding_exon"] == True]
    print(f"      coding exons: {len(cod):,} | exon_usage: "
          + ", ".join(f"{k} {v:,}" for k, v in cod['exon_usage'].value_counts(dropna=False).items()))
    print(f"      GENCODE-wide constitutive (reference): {int(cod['is_constitutive'].sum()):,} "
          f"({cod['is_constitutive'].mean()*100:.1f}%)")
    print(f"      chosen ENST is principal (MANE|APPRIS): {df.drop_duplicates('transcript_id')['is_principal'].mean()*100:.1f}% of transcripts")
    print(f"      IDR exons (idr_frac>=0.5): {int((cod['idr_frac'] >= 0.5).sum()):,}")
    print(f"      -> {OUT_PARQUET.name} + .csv + {OUT_DICT.name}")


def write_dictionary(df):
    rows = [
        ("gene_id", "Ensembl gene ID (ENSG, unversioned)"),
        ("gene_name", "HGNC gene symbol"),
        ("transcript_id", "Chosen/resolved Ensembl transcript (ENST) for this isoform"),
        ("protein_id", "Ensembl protein ID (ENSP) of the chosen transcript"),
        ("exon_id", "Ensembl exon ID (ENSE, unversioned)"),
        ("isoform_accession", "UniProt isoform accession (e.g. P04637 or P04637-2)"),
        ("base_accession", "UniProt gene-level accession (dash stripped)"),
        ("is_tf", "True if the gene is a curated transcription factor"),
        ("tf_family", "TF DNA-binding-domain family"),
        ("mapping_status", "Resolved mapping status (MAIN set: resolved_unique / resolved_collapsed_equivalent)"),
        ("chosen_enst", "The resolved transcript (== transcript_id); explicit alias"),
        ("n_candidate_ENST", "Number of candidate ENSTs UniProt linked to this isoform"),
        ("chr / start / end / strand", "Exon genomic coordinates (GRCh38/hg38, 1-based inclusive)"),
        ("exon_number", "Exon rank within the transcript (1 = first, 5'->3')"),
        ("cds_start / cds_end", "Coding-portion genomic coordinates (null for non-coding/UTR exons)"),
        ("cds_phase", "CDS phase (0/1/2) at the exon's coding start"),
        ("aa_start / aa_end", "Amino-acid range this exon encodes (1-based; null if non-coding)"),
        ("is_coding_exon", "True if the exon overlaps the CDS"),
        ("exon_length_bp / cds_length_bp / aa_length", "Exon / coding / amino-acid lengths"),
        ("exon_usage", "**Paper definition.** alternative = in some but not all of the gene's mapped transcripts; "
                       "constitutive = in all; single_tx_gene = gene has one mapped transcript; "
                       "null = isoform not in the s9 table (from analysis/exon_idr_annotation.csv)"),
        ("is_constitutive", "Reference only: present in ALL protein-coding transcripts of the gene (GENCODE-wide)"),
        ("transcript_inclusion_fraction", "n_coding_transcripts_with_exon / n_coding_transcripts_gene"),
        ("n_coding_transcripts_with_exon", "# gene's coding transcripts containing this exon (by ENSE)"),
        ("n_coding_transcripts_gene", "# protein-coding transcripts of the gene in GENCODE v" + C.GENCODE_RELEASE),
        ("is_principal", "Chosen ENST is MANE_Select OR APPRIS principal"),
        ("is_mane_select", "Chosen ENST carries the MANE_Select tag"),
        ("is_ensembl_canonical", "Chosen ENST carries the Ensembl_canonical tag"),
        ("is_appris_principal", "Chosen ENST carries an APPRIS principal tag"),
        ("appris_principal_level", "APPRIS principal level (1=best; null if not principal)"),
        ("idr_frac", "Fraction of this exon's residues that fall in a MetaPredict IDR"),
        ("region_class", "entirely_IDR / entirely_nonIDR / crosses_IDR_boundary"),
        ("is_idr_exon", "idr_frac >= 0.5"),
        ("is_canonical", "The isoform is the UniProt canonical sequence"),
        ("genomic_scaffold_validated", "Isoform's CDS round-tripped to protein in s4b"),
        ("is_first_exon / is_last_exon", "Exon is first / last in the transcript"),
        ("data_release_* / pipeline_run_date", "Provenance stamp"),
    ]
    lines = [
        "# `exon_level_table` — data dictionary",
        "",
        f"One row per **(UniProt isoform, exon of its chosen transcript)** for the MAIN "
        f"resolved set. GENCODE v{C.GENCODE_RELEASE} / {C.GENOME_BUILD}. "
        f"{len(df):,} rows, {df['isoform_accession'].nunique():,} isoforms, "
        f"{df['gene_id'].nunique():,} genes.",
        "",
        "| Column | Meaning |",
        "|---|---|",
    ]
    lines += [f"| `{c}` | {d} |" for c, d in rows]
    lines += [
        "",
        "**Notes.** Analyses use `exon_usage` for alternative vs constitutive exons. "
        "`is_constitutive` / `transcript_inclusion_fraction` use **all "
        "GENCODE protein-coding transcripts of the gene** as the denominator (not just "
        "our mapped isoforms) and are kept for reference. Non-coding (UTR) exons have null `cds_*` / `aa_*` / "
        "`idr_frac`. Coordinates are 1-based inclusive (GTF convention); the BED exports "
        "in `processed/*.bed` are 0-based half-open.",
    ]
    OUT_DICT.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
