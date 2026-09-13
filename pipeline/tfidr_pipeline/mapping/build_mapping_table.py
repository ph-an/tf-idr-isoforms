"""
STEP 2 of the workflow — build the isoform -> transcript mapping table.

Joins our 21,801 UniProt isoforms (IDRisoforms_df_geneage.csv) to the
isoform-resolved Ensembl cross-references parsed in step 1, and assigns a
`mapping_status` to every isoform:

    exact_one_to_one            isoform -> exactly one ENST
    multiple_candidate_transcripts  isoform -> >1 ENST (real: e.g. PAR genes,
                                    readthrough, or several Ensembl models of
                                    the same protein)
    gene_only_mapping           no isoform-specific ENST, but the gene (ENSG)
                                is known -> genomic locus available, transcript
                                identity unresolved
    no_mapping                  base accession absent from UniProt Ensembl xrefs
    sequence_mismatch           reserved; set in step 3 when an ENST translation
                                fails to match the UniProt isoform sequence

Canonical reconciliation:
    Our canonical isoforms are stored bare ("P04637", is_canonical=True) whereas
    UniProt tags them "P04637-1". For canonical rows we therefore accept either
    the bare accession OR "<base>-1" (and the untagged single-isoform form).

Outputs (written next to this script):
    isoform_transcript_map.csv             one row per isoform (primary pick + status)
    isoform_transcript_candidates_long.csv one row per candidate ENST
"""
import sys, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from tfidr_pipeline import config as C

ISO_CSV   = C.ISOFORM_ANNOTATED
PAIRS_CSV = C.CACHE / "uniprot_ensembl_isoform_pairs.csv"

# ── Load inputs ───────────────────────────────────────────────────────────────
iso = pd.read_csv(ISO_CSV)[
    ["isoform_accession", "base_accession", "is_canonical", "gene_name",
     "tf_group", "Length"]
].copy()
iso = iso.rename(columns={"isoform_accession": "uniprot_isoform"})

pairs = pd.read_csv(PAIRS_CSV)

# Fast lookups
#   by isoform tag  -> list of candidate ENST/ENSP/ENSG rows
pairs_by_tag = {tag: g for tag, g in pairs.groupby("isoform_tag")}
#   by base        -> any ENSG (gene-level fallback)
ensg_by_base = (pairs.groupby("base_accession")["ENSG"]
                .agg(lambda s: sorted(s.unique())).to_dict())


def candidate_tags(row):
    """Which UniProt isoform tag(s) could correspond to this dataset isoform."""
    acc  = row["uniprot_isoform"]
    base = row["base_accession"]
    if row["is_canonical"]:
        # canonical stored bare; UniProt tags it as -1 (or untagged single iso)
        return [acc, base, f"{base}-1"]
    return [acc]


# ── Classify every isoform ────────────────────────────────────────────────────
map_rows   = []
long_rows  = []

for _, row in iso.iterrows():
    acc, base = row["uniprot_isoform"], row["base_accession"]

    hit = None
    used_tag = None
    for tag in candidate_tags(row):
        if tag in pairs_by_tag:
            hit = pairs_by_tag[tag]
            used_tag = tag
            break

    if hit is not None:
        ensts = hit.drop_duplicates("ENST")
        n = len(ensts)
        status = "exact_one_to_one" if n == 1 else "multiple_candidate_transcripts"
        primary = ensts.iloc[0]
        map_rows.append({
            "uniprot_isoform": acc, "base_accession": base,
            "is_canonical": row["is_canonical"], "gene_name": row["gene_name"],
            "tf_group": row["tf_group"], "uniprot_length": row["Length"],
            "matched_tag": used_tag,
            "ENSG": primary["ENSG"], "ENST": primary["ENST"], "ENSP": primary["ENSP"],
            "n_candidate_ENST": n,
            "all_candidate_ENST": ";".join(ensts["ENST"].tolist()),
            "mapping_status": status,
        })
        for _, c in ensts.iterrows():
            long_rows.append({
                "uniprot_isoform": acc, "base_accession": base,
                "matched_tag": used_tag,
                "ENSG": c["ENSG"], "ENST": c["ENST"], "ENSP": c["ENSP"],
                "mapping_status": status,
            })
    else:
        # no isoform-specific transcript — is the gene at least known?
        if base in ensg_by_base:
            ensgs = ensg_by_base[base]
            map_rows.append({
                "uniprot_isoform": acc, "base_accession": base,
                "is_canonical": row["is_canonical"], "gene_name": row["gene_name"],
                "tf_group": row["tf_group"], "uniprot_length": row["Length"],
                "matched_tag": None,
                "ENSG": ";".join(ensgs), "ENST": None, "ENSP": None,
                "n_candidate_ENST": 0, "all_candidate_ENST": None,
                "mapping_status": "gene_only_mapping",
            })
        else:
            map_rows.append({
                "uniprot_isoform": acc, "base_accession": base,
                "is_canonical": row["is_canonical"], "gene_name": row["gene_name"],
                "tf_group": row["tf_group"], "uniprot_length": row["Length"],
                "matched_tag": None,
                "ENSG": None, "ENST": None, "ENSP": None,
                "n_candidate_ENST": 0, "all_candidate_ENST": None,
                "mapping_status": "no_mapping",
            })

map_df  = pd.DataFrame(map_rows)
long_df = pd.DataFrame(long_rows)

map_df.to_csv(C.ISOFORM_MAP, index=False)
long_df.to_csv(C.ISOFORM_CANDIDATES, index=False)

# ── Report ────────────────────────────────────────────────────────────────────
print(f"Isoforms processed: {len(map_df):,}\n")
print("mapping_status breakdown:")
vc = map_df["mapping_status"].value_counts()
for k, v in vc.items():
    print(f"  {k:<32} {v:>7,}  ({v/len(map_df)*100:5.1f}%)")

print("\nBy canonical vs alternative:")
print(pd.crosstab(map_df["is_canonical"], map_df["mapping_status"]).to_string())

print("\nBy TF vs non-TF:")
print(pd.crosstab(map_df["tf_group"], map_df["mapping_status"]).to_string())

resolved = map_df["mapping_status"].isin(["exact_one_to_one", "multiple_candidate_transcripts"]).sum()
print(f"\nIsoforms with a transcript-level (ENST) mapping: {resolved:,} "
      f"({resolved/len(map_df)*100:.1f}%)")
print(f"Saved -> isoform_transcript_map.csv  ({len(map_df):,} rows)")
print(f"Saved -> isoform_transcript_candidates_long.csv  ({len(long_df):,} rows)")
