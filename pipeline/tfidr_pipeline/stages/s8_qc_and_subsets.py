"""
QC report + analysis subsets (spec sections 6-8).

Writes:
  qc/qc_report.md                      human-readable summary
  qc/mapping_status_counts.csv
  qc/missingness_by_column.csv
  qc/ambiguous_mapping_examples.csv
  analysis/analysis_isoforms_all.parquet
  analysis/analysis_isoforms_unique_enst_only.parquet
  analysis/analysis_exon_level_validated_only.parquet
  analysis/analysis_tf_only.parquet
  analysis/analysis_nontf_only.parquet
  analysis/analysis_genes_all.parquet
  analysis/analysis_genes_min3_isoforms.parquet
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from tfidr_pipeline import config as C


def main():
    m    = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    gene = pd.read_parquet(C.PROCESSED / "tf_idr_gene_summary.parquet")
    n = len(m)

    # ── QC CSVs ───────────────────────────────────────────────────────────────
    (m["mapping_level"].value_counts(dropna=False)
     .rename_axis("mapping_level").reset_index(name="n_isoforms")
     .to_csv(C.QC / "mapping_status_counts.csv", index=False))

    miss = (m.isna().mean().mul(100).round(2)
            .rename_axis("column").reset_index(name="pct_missing")
            .sort_values("pct_missing", ascending=False))
    miss.to_csv(C.QC / "missingness_by_column.csv", index=False)

    amb = (m[m["has_multiple_ensts"]]
           .nlargest(50, "n_candidate_ENST")
           [["isoform_accession", "base_accession", "gene_name", "tf_group",
             "n_candidate_ENST", "sequence_match_type", "is_validated_mapping"]])
    amb.to_csv(C.QC / "ambiguous_mapping_examples.csv", index=False)

    # ── analysis subsets ──────────────────────────────────────────────────────
    m.to_parquet(C.ANALYSIS / "analysis_isoforms_all.parquet", index=False)
    m[m["has_unique_enst"] & m["is_validated_mapping"]].to_parquet(
        C.ANALYSIS / "analysis_isoforms_unique_enst_only.parquet", index=False)
    m[m["analysis_include_exon_level"]].to_parquet(
        C.ANALYSIS / "analysis_exon_level_validated_only.parquet", index=False)
    m[m["is_tf"]].to_parquet(C.ANALYSIS / "analysis_tf_only.parquet", index=False)
    m[~m["is_tf"]].to_parquet(C.ANALYSIS / "analysis_nontf_only.parquet", index=False)
    gene.to_parquet(C.ANALYSIS / "analysis_genes_all.parquet", index=False)
    gene[gene["n_isoforms_gene"] >= C.MIN_ISOFORMS_FOR_GENE_DELTA].to_parquet(
        C.ANALYSIS / "analysis_genes_min3_isoforms.parquet", index=False)

    # ── QC metrics ────────────────────────────────────────────────────────────
    def pct(x):
        return f"{x:,} ({x/n*100:.1f}%)"

    has_enst = m["ensembl_transcript_id"].notna()
    lines = [
        "# QC Report — TF-IDR isoform genomic-coordinate pipeline",
        f"_Generated {C.RELEASE_META['pipeline_run_date']} · "
        f"Ensembl {C.RELEASE_META['data_release_ensembl']}_",
        "",
        "## Isoform backbone",
        f"- Total isoforms: **{n:,}**",
        f"- Base accessions (genes): **{m['base_accession'].nunique():,}**",
        f"- Canonical: {pct(m['is_canonical'].sum())}",
        f"- Non-canonical: {pct((~m['is_canonical']).sum())}",
        f"- Fraction with sequence: {m['sequence'].notna().mean()*100:.1f}%",
        "",
        "## Gene classification",
        f"- TF genes: {gene['is_tf'].sum():,}  |  Non-TF genes: {(~gene['is_tf']).sum():,}",
        f"- Isoforms in TF genes: {pct(m['is_tf'].sum())}",
        "",
        "## UniProt → Ensembl mapping",
        f"- Has Ensembl gene ID (ENSG): {pct(m['ensembl_gene_id'].notna().sum())}",
        f"- Unique ENST (exact_one_to_one): {pct(m['has_unique_enst'].sum())}",
        f"- Multiple candidate ENSTs: {pct(m['has_multiple_ensts'].sum())}",
        f"- Gene-only mapping: {pct(m['has_gene_only_mapping'].sum())}",
        f"- Unmapped: {pct((m['mapping_level'] == 'no_mapping').sum())}",
        "",
        "## Sequence validation (GENCODE translation vs UniProt)",
        f"- Validated (exact only) among isoforms with an ENST: "
        f"{m.loc[has_enst,'is_validated_mapping'].sum():,}/{has_enst.sum():,} "
        f"({m.loc[has_enst,'is_validated_mapping'].mean()*100:.1f}%)",
        f"- Exon-level analysis-eligible: {pct(m['analysis_include_exon_level'].sum())}",
        "",
        "## Downstream annotation coverage",
        f"- Gene age present: {m['gene_age'].notna().mean()*100:.1f}%",
        f"- APPRIS annotation present: {m['appris_annotation'].notna().mean()*100:.1f}%",
        f"- Exon-boundary density computed (IDR): {m['idr_boundary_density_per_kb'].notna().mean()*100:.1f}%",
        "",
        "## Caveats",
        "- **APPRIS is gene-level (ENSG), not transcript-level** in the source "
        "annotation, so `is_appris_principal` marks *genes that have a principal "
        "isoform*, applied to every isoform of that gene — it is not a per-isoform "
        "principal call. `canonical_appris_agree` should be read with that in mind.",
        "- Exon-level tables are restricted to sequence-validated isoforms so that "
        "UniProt aa coordinates align with the Ensembl translation the exon map is "
        "built on.",
        "- `multiple_candidate_transcripts` rows use the first UniProt-listed ENST as "
        "primary; validation shows this is the correct sequence "
        f"{m.loc[m['has_multiple_ensts'],'is_validated_mapping'].mean()*100:.1f}% of the time. "
        "Use `analysis_isoforms_unique_enst_only` for the strictest coordinate claims.",
        "",
        "## Files",
        "- `processed/tf_idr_isoform_master.{parquet,csv}` — isoform master",
        "- `processed/tf_idr_gene_summary.{parquet,csv}` — gene summary",
        "- `analysis/viz_*.csv` — figure-ready tables",
        "- `analysis/analysis_*.parquet` — pre-filtered safe subsets",
        "- `interim/*.parquet` — exon coords, CDS↔protein map, overlap tables, validation",
    ]
    (C.QC / "qc_report.md").write_text("\n".join(lines), encoding="utf-8")

    print("Wrote QC report + subsets.\n")
    print("\n".join(lines[:34]))


if __name__ == "__main__":
    main()
