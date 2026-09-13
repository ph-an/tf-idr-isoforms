"""
Stage K — assemble the isoform-level master table.

Reuses the validated annotation already in the project and layers on the new
genomic-coordinate results:

  IDRisoforms_df_geneage.csv        backbone + TF/family + gene age + IDR features
  IDRisoforms_df_appris.csv         APPRIS principal/alternative + ENSG
  isoform_transcript_map.csv        UniProt->Ensembl mapping + candidate counts
  mapping_validation_scaled.parquet at-scale sequence validation
  idr/nonidr_exon_overlap tables    per-isoform exon-boundary aggregates

Adds mapping-confidence flags, APPRIS-vs-canonical agreement, release
provenance, and per-analysis inclusion flags so imperfect rows can be kept but
excluded from analyses that would be invalid for them.

Outputs:
    processed/tf_idr_isoform_master.parquet
    processed/tf_idr_isoform_master.csv
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C

OUT_PARQUET = C.PROCESSED / "tf_idr_isoform_master.parquet"
OUT_CSV     = C.PROCESSED / "tf_idr_isoform_master.csv"


def main():
    # ── backbone + features (already validated annotation) ────────────────────
    base = pd.read_csv(C.ISOFORM_ANNOTATED)
    keep = [c for c in [
        "isoform_accession", "base_accession", "gene_name", "is_canonical",
        "Sequence", "Length", "tf_group", "tf_family",
        "pct_idr", "idr_aaLen", "n_idr_segments", "max_idr_len", "mean_idr_len",
        "pct_change_idr", "gene_age", "age_bin", "n_isoforms",
        "mean_pct_idr_all_isoforms",
    ] if c in base.columns]
    df = base[keep].rename(columns={"Sequence": "sequence",
                                    "Length": "protein_length_aa"})
    df["is_tf"] = df["tf_group"].eq("TF")

    # ── APPRIS (gene/transcript-level principal annotation) ───────────────────
    appr = pd.read_csv(C.ISOFORM_APPRIS,
                       usecols=lambda c: c in {"isoform_accession", "ensg_id",
                                               "appris_annotation", "appris_tier",
                                               "protein_length_appris"})
    df = df.merge(appr, on="isoform_accession", how="left")
    df["is_appris_principal"]  = df["appris_annotation"].astype(str).str.upper().str.contains("PRINCIPAL")
    df["is_appris_alternative"] = df["appris_annotation"].astype(str).str.upper().str.contains("ALTERNATIVE")

    # ── UniProt->Ensembl mapping + candidate counts ───────────────────────────
    m = pd.read_csv(C.ISOFORM_MAP)[
        ["uniprot_isoform", "ENSG", "ENST", "ENSP",
         "n_candidate_ENST", "mapping_status"]
    ].rename(columns={"uniprot_isoform": "isoform_accession",
                      "ENSG": "ensembl_gene_id", "ENST": "ensembl_transcript_id",
                      "ENSP": "ensembl_protein_id", "mapping_status": "mapping_level"})
    df = df.merge(m, on="isoform_accession", how="left")
    df["has_unique_enst"]      = df["mapping_level"].eq("exact_one_to_one")
    df["has_multiple_ensts"]   = df["mapping_level"].eq("multiple_candidate_transcripts")
    df["has_gene_only_mapping"] = df["mapping_level"].eq("gene_only_mapping")

    # ── sequence validation (at scale) + resolved transcript (s3) ──────────────
    val_cols = ["uniprot_isoform", "sequence_match_type", "sequence_identity",
                "is_validated_mapping", "ensembl_translation_length_aa",
                "length_difference_aa", "resolved_transcript",
                "resolved_is_validated", "resolved_by_validation"]
    val = pd.read_parquet(C.INTERIM / "mapping_validation_scaled.parquet")
    val = val[[c for c in val_cols if c in val.columns]].rename(
        columns={"uniprot_isoform": "isoform_accession"})
    df = df.merge(val, on="isoform_accession", how="left")
    df["is_validated_mapping"] = df["is_validated_mapping"].fillna(False)
    if "resolved_is_validated" in df.columns:
        df["resolved_is_validated"] = df["resolved_is_validated"].fillna(False)

    # ── APPRIS vs UniProt canonical agreement ─────────────────────────────────
    df["is_uniprot_canonical"] = df["is_canonical"]
    df["canonical_appris_agree"] = np.where(
        df["appris_annotation"].isna(), np.nan,
        (df["is_uniprot_canonical"] & df["is_appris_principal"]) |
        (~df["is_uniprot_canonical"] & ~df["is_appris_principal"])
    )

    # ── per-isoform exon-boundary aggregates (validated only) ──────────────────
    idr = pd.read_parquet(C.INTERIM / "idr_exon_overlap_table.parquet")
    idr_agg = idr.groupby("isoform_accession").agg(
        n_idr_segments_mapped=("idr_segment_id", "nunique"),
        n_overlapping_idr_exons=("n_overlapping_coding_exons", "sum"),
        idr_contains_exon_boundary=("contains_exon_boundary", "any"),
        idr_boundary_density_per_kb=("idr_boundary_density_per_kb", "mean"),
    ).reset_index()

    nonidr = pd.read_parquet(C.INTERIM / "nonidr_exon_overlap_table.parquet")
    nonidr_agg = nonidr.groupby("isoform_accession").agg(
        nonidr_boundary_density_per_kb=("nonidr_boundary_density_per_kb", "mean"),
    ).reset_index()

    df = df.merge(idr_agg, on="isoform_accession", how="left")
    df = df.merge(nonidr_agg, on="isoform_accession", how="left")
    df["idr_vs_nonidr_boundary_density_delta"] = (
        df["idr_boundary_density_per_kb"] - df["nonidr_boundary_density_per_kb"])

    # ── true GENOMIC projection (s4b): scaffold validation + genomic boundaries ─
    scaf_path = C.INTERIM / "genomic_scaffold_validation.parquet"
    seg_path  = C.INTERIM / "genomic_projection_segments.parquet"
    if scaf_path.exists() and seg_path.exists():
        scaf = (pd.read_parquet(scaf_path)
                [["isoform_accession", "genomic_scaffold_validated", "cds_len_coherent"]]
                .drop_duplicates("isoform_accession"))
        df = df.merge(scaf, on="isoform_accession", how="left")
        df["genomic_scaffold_validated"] = df["genomic_scaffold_validated"].fillna(False)

        seg = pd.read_parquet(seg_path)
        gidr = (seg[seg["segment_type"] == "idr"].groupby("isoform_accession").agg(
                    idr_boundary_density_per_kb_genomic=("boundary_density_per_kb", "mean"),
                    idr_genomic_pieces_mean=("n_genomic_pieces", "mean"),
                    idr_projection_validated_frac=("projection_validated", "mean"))
                .reset_index())
        gnon = (seg[seg["segment_type"] == "nonidr"].groupby("isoform_accession").agg(
                    nonidr_boundary_density_per_kb_genomic=("boundary_density_per_kb", "mean"))
                .reset_index())
        df = df.merge(gidr, on="isoform_accession", how="left")
        df = df.merge(gnon, on="isoform_accession", how="left")

        # cross-check: aa-space (s4) vs genomic (s4b) IDR boundary density should agree
        df["boundary_density_aa_vs_genomic_delta"] = (
            df["idr_boundary_density_per_kb"] - df["idr_boundary_density_per_kb_genomic"])
        df["genomic_idr_vs_nonidr_boundary_density_delta"] = (
            df["idr_boundary_density_per_kb_genomic"]
            - df["nonidr_boundary_density_per_kb_genomic"])
    else:
        df["genomic_scaffold_validated"] = False

    # ── long-IDR flags ────────────────────────────────────────────────────────
    for thr in C.LONG_IDR_THRESHOLDS_AA:
        df[f"has_long_idr_{thr}aa"] = df["max_idr_len"].fillna(0) >= thr
    df["has_idr"] = df["n_idr_segments"].fillna(0) > 0

    # ── analysis inclusion flags ──────────────────────────────────────────────
    df["analysis_include_isoform_level"] = True    # all isoforms ok at protein level
    df["analysis_include_exon_level"]    = df["is_validated_mapping"]
    # genomic-coordinate claims require the self-verified DNA round-trip (s4b)
    df["analysis_include_genomic_projection"] = df["genomic_scaffold_validated"]
    df["exclude_reason"] = np.where(
        df["is_validated_mapping"], "",
        np.where(df["ensembl_transcript_id"].isna(),
                 "no_enst_" + df["mapping_level"].fillna("unmapped").astype(str),
                 "enst_not_sequence_validated"))

    df = C.stamp(df)
    C.PROCESSED.mkdir(exist_ok=True)
    df.to_parquet(OUT_PARQUET, index=False)
    df.to_csv(OUT_CSV, index=False)

    print(f"[out] master {len(df):,} isoforms x {df.shape[1]} cols")
    print(f"      -> {OUT_PARQUET.name} + .csv\n")
    print("exon-level analysis-eligible (validated):",
          f"{df['analysis_include_exon_level'].sum():,} "
          f"({df['analysis_include_exon_level'].mean()*100:.1f}%)")
    print("genomic-projection eligible (scaffold self-verified):",
          f"{df['analysis_include_genomic_projection'].sum():,} "
          f"({df['analysis_include_genomic_projection'].mean()*100:.1f}%)")
    if "boundary_density_aa_vs_genomic_delta" in df.columns:
        d = df["boundary_density_aa_vs_genomic_delta"].dropna()
        print(f"aa-space vs genomic IDR boundary density delta: "
              f"mean={d.mean():.4f}, median={d.median():.4f} (per kb; ~0 = methods agree)")
    print("APPRIS principal isoforms:", int(df['is_appris_principal'].sum()))
    ok = df["canonical_appris_agree"].dropna()
    print(f"canonical<->APPRIS agreement: {ok.mean()*100:.1f}% "
          f"(of {len(ok):,} isoforms with APPRIS)")
    print("\ncolumns:", list(df.columns))


if __name__ == "__main__":
    main()
