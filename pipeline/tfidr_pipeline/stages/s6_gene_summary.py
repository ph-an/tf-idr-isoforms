"""
Gene-level summary table — one row per base accession (gene/protein).

Aggregates the isoform master into the gene-level dIDR variables that drive the
TF-vs-Non-TF comparisons. The headline variable is:

    pct_change_idr_gene = max(pct_idr across isoforms) - min(pct_idr across isoforms)

Outputs:
    processed/tf_idr_gene_summary.parquet
    processed/tf_idr_gene_summary.csv
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C

OUT_PARQUET = C.PROCESSED / "tf_idr_gene_summary.parquet"
OUT_CSV     = C.PROCESSED / "tf_idr_gene_summary.csv"


def main():
    m = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")

    def agg_gene(g):
        canon = g[g["is_canonical"]]
        noncanon = g[~g["is_canonical"]]
        pct = g["pct_idr"].dropna()
        return pd.Series({
            "gene_name": g["gene_name"].iloc[0],
            "ensembl_gene_id": g["ensembl_gene_id"].dropna().iloc[0] if g["ensembl_gene_id"].notna().any() else None,
            "is_tf": bool(g["is_tf"].iloc[0]),
            "tf_group": g["tf_group"].iloc[0],
            "tf_family": g["tf_family"].iloc[0] if "tf_family" in g else None,
            "gene_age": g["gene_age"].iloc[0],
            "age_bin": g["age_bin"].iloc[0],

            "n_isoforms_gene": g["isoform_accession"].nunique(),
            "n_mapped_isoforms": g["ensembl_transcript_id"].notna().sum(),
            "n_uniquely_mapped_isoforms": g["has_unique_enst"].sum(),
            "n_ambiguous_mapped_isoforms": g["has_multiple_ensts"].sum(),
            "n_gene_only_mapped_isoforms": g["has_gene_only_mapping"].sum(),
            "n_validated_isoforms": g["is_validated_mapping"].sum(),

            "mean_pct_idr_gene": pct.mean(),
            "median_pct_idr_gene": pct.median(),
            "min_pct_idr_gene": pct.min(),
            "max_pct_idr_gene": pct.max(),
            "pct_change_idr_gene": (pct.max() - pct.min()) if len(pct) else np.nan,

            "canonical_pct_idr": canon["pct_idr"].mean() if len(canon) else np.nan,
            "mean_noncanonical_pct_idr": noncanon["pct_idr"].mean() if len(noncanon) else np.nan,

            "mean_protein_length_gene": g["protein_length_aa"].mean(),
            "canonical_length": canon["protein_length_aa"].mean() if len(canon) else np.nan,
            "length_range_gene": (g["protein_length_aa"].max() - g["protein_length_aa"].min()),

            "mean_idr_boundary_density_per_kb": g["idr_boundary_density_per_kb"].mean(),
            "mean_nonidr_boundary_density_per_kb": g["nonidr_boundary_density_per_kb"].mean(),

            "has_any_idr_isoform": bool((g["n_idr_segments"].fillna(0) > 0).any()),
        })

    gene = m.groupby("base_accession", sort=False).apply(agg_gene, include_groups=False).reset_index()

    gene["delta_canonical_vs_mean_noncanonical_pct_idr"] = (
        gene["canonical_pct_idr"] - gene["mean_noncanonical_pct_idr"])
    gene["idr_vs_nonidr_boundary_density_delta"] = (
        gene["mean_idr_boundary_density_per_kb"] - gene["mean_nonidr_boundary_density_per_kb"])
    gene["has_idr_gain_loss"] = gene["pct_change_idr_gene"] > 0
    gene["has_high_delta_idr"] = gene["pct_change_idr_gene"] >= 10

    # isoform-count bin (3,4,...,10+)
    def iso_bin(n):
        if n < C.MIN_ISOFORMS_FOR_GENE_DELTA:
            return f"<{C.MIN_ISOFORMS_FOR_GENE_DELTA}"
        return "10+" if n >= 10 else str(int(n))
    gene["isoform_count_bin"] = gene["n_isoforms_gene"].apply(iso_bin)

    gene = C.stamp(gene)
    gene.to_parquet(OUT_PARQUET, index=False)
    gene.to_csv(OUT_CSV, index=False)

    print(f"[out] gene summary {len(gene):,} genes x {gene.shape[1]} cols -> {OUT_PARQUET.name} + .csv\n")
    print("TF vs Non-TF gene counts:")
    print(gene["tf_group"].value_counts().to_string())
    print("\nmedian pct_change_idr_gene (genes with >=3 isoforms):")
    sub = gene[gene["n_isoforms_gene"] >= C.MIN_ISOFORMS_FOR_GENE_DELTA]
    print(sub.groupby("tf_group")["pct_change_idr_gene"].median().to_string())


if __name__ == "__main__":
    main()
