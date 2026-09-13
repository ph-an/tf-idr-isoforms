"""
Stage D (part 2) — exon/CDS coordinate tables from the GENCODE GTF.

Two outputs, for every ENST referenced in the isoform->transcript mapping:

  interim/ensembl_transcript_exon_coordinates.parquet
      one row per exon: genomic exon interval, CDS interval, coding flag, rank

  interim/transcript_cds_to_protein_coordinate_map.parquet
      one row per CODING exon chunk, carrying the cumulative CDS offset and the
      amino-acid range that chunk encodes (strand-aware). This is what lets us
      project an IDR's aa interval onto specific coding exons in Stage F.

No genome FASTA required — only coordinates. CDS phase is used to keep the
aa<->bp bookkeeping honest across split codons.
"""
import sys, gzip, re
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C

EXON_OUT = C.INTERIM / "ensembl_transcript_exon_coordinates.parquet"
MAP_OUT  = C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet"

def parse_attrs(s):
    """GENCODE GTF attribute column -> dict. Handles quoted string values
    (gene_id "ENSG..") and unquoted numeric values (exon_number 1)."""
    d = {}
    for field in s.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        key, _, val = field.partition(" ")
        d[key] = val.strip().strip('"')
    return d


def wanted_enst():
    cand = pd.read_csv(C.ISOFORM_CANDIDATES)
    m    = pd.read_csv(C.ISOFORM_MAP)
    ids  = pd.concat([cand["ENST"], m["ENST"]]).dropna()
    return set(ids.str.split(".").str[0])


def parse_gtf(keep):
    rows = []
    with gzip.open(C.GENCODE_GTF_GZ, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "CDS"):
                continue
            attrs = parse_attrs(f[8])
            tid = attrs.get("transcript_id", "")
            enst = tid.split(".")[0]
            if enst not in keep:
                continue
            rows.append({
                "ensembl_transcript_id": enst,
                "ensembl_gene_id":  attrs.get("gene_id", "").split(".")[0],
                "ensembl_protein_id": attrs.get("protein_id", "").split(".")[0] or None,
                "feature":   f[2],
                "chromosome": f[0],
                "start":     int(f[3]),
                "end":       int(f[4]),
                "strand":    f[6],
                "phase":     None if f[7] == "." else int(f[7]),
                "exon_rank": int(attrs["exon_number"]) if "exon_number" in attrs else None,
                "exon_id":   attrs.get("exon_id", "").split(".")[0] or None,
                "transcript_biotype": attrs.get("transcript_type"),
            })
    return pd.DataFrame(rows)


def build_exon_table(gtf):
    """Wide per-exon table: exon interval + CDS interval joined on (ENST, rank)."""
    exons = gtf[gtf["feature"] == "exon"].copy()
    cds   = gtf[gtf["feature"] == "CDS"].copy()

    ex = exons.rename(columns={"start": "exon_genomic_start",
                               "end":   "exon_genomic_end"})
    ex["exon_length_bp"] = ex["exon_genomic_end"] - ex["exon_genomic_start"] + 1

    cds_small = cds.rename(columns={"start": "cds_genomic_start",
                                    "end":   "cds_genomic_end",
                                    "phase": "cds_phase"})[
        ["ensembl_transcript_id", "exon_rank",
         "cds_genomic_start", "cds_genomic_end", "cds_phase"]]

    out = ex.merge(cds_small, on=["ensembl_transcript_id", "exon_rank"], how="left")
    out["cds_length_bp"] = (out["cds_genomic_end"] - out["cds_genomic_start"] + 1)
    out["is_coding_exon"] = out["cds_genomic_start"].notna()
    out["genome_build"]   = C.GENOME_BUILD
    out["ensembl_release"] = f"GENCODE_v{C.GENCODE_RELEASE}"
    keep_cols = ["ensembl_gene_id", "ensembl_transcript_id", "ensembl_protein_id",
                 "chromosome", "strand", "exon_id", "exon_rank",
                 "exon_genomic_start", "exon_genomic_end", "exon_length_bp",
                 "cds_genomic_start", "cds_genomic_end", "cds_length_bp",
                 "cds_phase", "is_coding_exon", "transcript_biotype",
                 "genome_build", "ensembl_release"]
    return out[keep_cols].sort_values(["ensembl_transcript_id", "exon_rank"])


def build_cds_protein_map(exon_tbl):
    """Cumulative CDS offset + aa range per coding exon chunk, strand-aware.

    Exons are already ranked 5'->3' by GENCODE (rank 1 = first coding exon in
    transcription direction on either strand), so walking by ascending rank and
    accumulating CDS length gives the coding-sequence offset directly.

    CDS start-phase: GENCODE's GTF phase on the first coding exon is the number of
    bases before the first complete codon. It is 0 for a 5'-complete CDS (nearly
    all sequence-validated isoforms); for a 5'-incomplete CDS (cds_start_NF) it is
    1 or 2. We seed the cursor at -start_phase so aa 1 lands on the first complete
    codon in either case (a no-op when start_phase == 0). We also carry a per-
    transcript cds_bp_total so s4b/QC can assert length coherence against the
    protein (CDS bp == protein_len * 3), Susie's `n_missing_nt == 0` check.
    """
    coding = exon_tbl[exon_tbl["is_coding_exon"]].copy()
    coding = coding.sort_values(["ensembl_transcript_id", "exon_rank"])

    recs = []
    for enst, g in coding.groupby("ensembl_transcript_id", sort=False):
        g_rows = list(g.itertuples())
        p0 = g_rows[0].cds_phase
        start_phase = 0 if (p0 is None or pd.isna(p0)) else int(p0)
        cds_bp_total = int(sum(int(r.cds_length_bp) for r in g_rows))
        cds_bp_cursor = -start_phase                     # 0 for a 5'-complete CDS
        for r in g_rows:
            start_bp = cds_bp_cursor                        # 0-based into CDS
            end_bp   = cds_bp_cursor + int(r.cds_length_bp)  # exclusive
            # aa range this chunk touches (1-based, inclusive); codons may split
            aa_start = start_bp // 3 + 1
            aa_end   = (end_bp - 1) // 3 + 1
            recs.append({
                "ensembl_transcript_id": enst,
                "ensembl_protein_id":    r.ensembl_protein_id,
                "exon_id":   r.exon_id,
                "exon_rank": r.exon_rank,
                "chromosome": r.chromosome,
                "strand":    r.strand,
                "cds_genomic_start": r.cds_genomic_start,
                "cds_genomic_end":   r.cds_genomic_end,
                "cds_start_bp_in_cds": start_bp,
                "cds_end_bp_in_cds":   end_bp,
                "aa_start": aa_start,
                "aa_end":   aa_end,
                "frame":    start_bp % 3,
                "cds_start_phase": start_phase,
                "cds_bp_total":    cds_bp_total,
            })
            cds_bp_cursor = end_bp
    m = pd.DataFrame(recs)
    m["aa_length"] = m["aa_end"] - m["aa_start"] + 1
    return m


if __name__ == "__main__":
    keep = wanted_enst()
    print(f"[gtf] extracting exon/CDS for {len(keep):,} mapped ENSTs")
    gtf = parse_gtf(keep)
    print(f"[gtf] {len(gtf):,} exon+CDS records for {gtf['ensembl_transcript_id'].nunique():,} transcripts")

    exon_tbl = build_exon_table(gtf)
    exon_tbl.to_parquet(EXON_OUT, index=False)
    print(f"[out] exon table            {len(exon_tbl):,} rows -> {EXON_OUT.name}")

    cds_map = build_cds_protein_map(exon_tbl)
    cds_map.to_parquet(MAP_OUT, index=False)
    print(f"[out] cds->protein coord map {len(cds_map):,} rows -> {MAP_OUT.name}")

    # frame-coherence: CDS total should be a whole number of codons
    per_tx = cds_map.drop_duplicates("ensembl_transcript_id")
    n_incoherent = int((per_tx["cds_bp_total"] % 3 != 0).sum())
    n_incomplete = int((per_tx["cds_start_phase"] != 0).sum())
    print(f"[qc]  transcripts with CDS bp not a multiple of 3: {n_incoherent:,}")
    print(f"[qc]  transcripts with 5'-incomplete CDS (start_phase!=0): {n_incomplete:,}")

    print("\nsample exon rows:")
    print(exon_tbl.head(4).to_string(index=False))
    print("\nsample cds->protein rows:")
    print(cds_map.head(4).to_string(index=False))
