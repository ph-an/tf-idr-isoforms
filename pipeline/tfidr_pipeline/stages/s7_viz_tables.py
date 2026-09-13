"""
Visualization-ready tables (spec section 4) — thin, pre-filtered slices of the
master / gene tables, each shaped for one figure family. Saved as CSV for easy
plotting.

  A viz_tf_vs_nontf_delta_idr.csv          gene-level dIDR, TF vs Non-TF
  B viz_canonical_vs_noncanonical.csv      isoform pct_idr, canonical split
  C viz_isoform_count_controlled.csv       gene dIDR by isoform-count bin
  D viz_gene_age_tf_vs_nontf.csv           gene dIDR by age bin
  E viz_tf_family_summary.csv              per-family dIDR summary
  F viz_idr_exon_boundary.csv              per-segment IDR/non-IDR boundary density
  G viz_length_matched_segments.csv        per-segment, with aa length bin
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C

VIZ = C.ANALYSIS


def main():
    master = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    gene   = pd.read_parquet(C.PROCESSED / "tf_idr_gene_summary.parquet")

    # A — TF vs Non-TF gene-level dIDR
    a = gene[["base_accession", "gene_name", "tf_group", "tf_family",
              "n_isoforms_gene", "pct_change_idr_gene", "mean_pct_idr_gene",
              "gene_age", "age_bin"]].copy()
    a.to_csv(VIZ / "viz_tf_vs_nontf_delta_idr.csv", index=False)

    # B — canonical vs noncanonical isoform pct_idr
    b = master[["isoform_accession", "base_accession", "gene_name", "tf_group",
                "is_canonical", "pct_idr", "protein_length_aa",
                "n_idr_segments"]].copy()
    b.to_csv(VIZ / "viz_canonical_vs_noncanonical.csv", index=False)

    # C — isoform-count-controlled gene dIDR
    c = gene[["base_accession", "gene_name", "tf_group", "n_isoforms_gene",
              "isoform_count_bin", "pct_change_idr_gene", "mean_pct_idr_gene"]].copy()
    c.to_csv(VIZ / "viz_isoform_count_controlled.csv", index=False)

    # D — gene-age stratified dIDR
    d = gene[["base_accession", "gene_name", "tf_group", "gene_age", "age_bin",
              "pct_change_idr_gene", "mean_pct_idr_gene", "n_isoforms_gene"]].copy()
    d.to_csv(VIZ / "viz_gene_age_tf_vs_nontf.csv", index=False)

    # E — TF-family summary
    tf_genes = gene[gene["is_tf"]]
    e = (tf_genes.groupby("tf_family")
         .agg(n_genes=("base_accession", "nunique"),
              mean_pct_change_idr=("pct_change_idr_gene", "mean"),
              median_pct_change_idr=("pct_change_idr_gene", "median"),
              mean_pct_idr=("mean_pct_idr_gene", "mean"),
              median_pct_idr=("mean_pct_idr_gene", "median"),
              mean_n_isoforms=("n_isoforms_gene", "mean"))
         .reset_index().sort_values("median_pct_change_idr", ascending=False))
    e.to_csv(VIZ / "viz_tf_family_summary.csv", index=False)

    # F & G — per-segment boundary density (IDR + non-IDR), length-binned
    idr = pd.read_parquet(C.INTERIM / "idr_exon_overlap_table.parquet")
    non = pd.read_parquet(C.INTERIM / "nonidr_exon_overlap_table.parquet")

    tf_map = master.set_index("isoform_accession")[["tf_group", "tf_family", "gene_name"]]

    idr_seg = idr.rename(columns={
        "idr_length_aa": "segment_length_aa",
        "idr_boundary_density_per_kb": "boundary_density_per_kb",
        "n_overlapping_coding_exons": "n_overlapping_exons"})
    idr_seg["segment_type"] = "IDR"
    idr_seg["segment_id"] = idr_seg["idr_segment_id"]

    non_seg = non.rename(columns={
        "nonidr_length_aa": "segment_length_aa",
        "nonidr_boundary_density_per_kb": "boundary_density_per_kb",
        "n_overlapping_coding_exons": "n_overlapping_exons"})
    non_seg["segment_type"] = "non-IDR"
    non_seg["segment_id"] = non_seg["nonidr_segment_id"]

    cols = ["segment_id", "isoform_accession", "base_accession",
            "segment_type", "segment_length_aa", "boundary_density_per_kb",
            "n_overlapping_exons", "contains_exon_boundary"]
    seg = pd.concat([idr_seg[cols], non_seg[cols]], ignore_index=True)
    seg = seg.join(tf_map, on="isoform_accession")
    seg["segment_length_bp"] = seg["segment_length_aa"] * 3
    seg["segment_length_bin"] = pd.cut(seg["segment_length_aa"],
                                       bins=C.SEGMENT_LEN_AA_BINS,
                                       labels=C.SEGMENT_LEN_AA_LABELS, right=False)

    # F — exon-boundary table
    seg[["isoform_accession", "base_accession", "gene_name", "tf_group",
         "tf_family", "segment_type", "segment_length_aa", "segment_length_bp",
         "n_overlapping_exons", "contains_exon_boundary",
         "boundary_density_per_kb"]].to_csv(VIZ / "viz_idr_exon_boundary.csv", index=False)

    # G — length-matched segment table
    seg[["segment_id", "isoform_accession", "base_accession", "gene_name",
         "tf_group", "segment_type", "segment_length_aa", "segment_length_bin",
         "boundary_density_per_kb", "n_overlapping_exons"]].to_csv(
        VIZ / "viz_length_matched_segments.csv", index=False)

    print("Wrote viz tables:")
    for f in ["viz_tf_vs_nontf_delta_idr", "viz_canonical_vs_noncanonical",
              "viz_isoform_count_controlled", "viz_gene_age_tf_vs_nontf",
              "viz_tf_family_summary", "viz_idr_exon_boundary",
              "viz_length_matched_segments"]:
        n = len(pd.read_csv(VIZ / f"{f}.csv"))
        print(f"  {f}.csv  ({n:,} rows)")

    # quick sanity: length-matched IDR vs non-IDR median density per bin
    print("\nLength-matched median boundary density (per kb), IDR vs non-IDR:")
    piv = seg.pivot_table(index="segment_length_bin", columns="segment_type",
                          values="boundary_density_per_kb", aggfunc="median",
                          observed=True)
    print(piv.to_string())


if __name__ == "__main__":
    main()
