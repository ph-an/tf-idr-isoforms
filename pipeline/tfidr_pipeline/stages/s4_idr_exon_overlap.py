"""
Stage F — project IDR (and non-IDR control) segments onto coding exons.

For every sequence-validated isoform we take its MetaPredict IDR segments
(aa intervals) and its ENST's CDS->protein coordinate map (from s2), then count
how many coding exons each segment spans and how many internal exon boundaries
fall inside it. Non-IDR segments are the complement of the IDR union over the
protein length, giving a length-matched within-protein control.

Boundary density (identical formula for IDR and non-IDR):
    boundary_density_per_kb = ((n_exon_chunks - 1) / segment_length_bp) * 1000
    segment_length_bp       = segment_length_aa * 3

Restricted to is_validated_mapping == True so the aa coordinates (UniProt) align
with the Ensembl translation the exon map is built on.

Outputs:
    interim/idr_exon_overlap_table.parquet
    interim/nonidr_exon_overlap_table.parquet
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from tfidr_pipeline import config as C
from tfidr_pipeline import utils as U

IDR_OUT    = C.INTERIM / "idr_exon_overlap_table.parquet"
NONIDR_OUT = C.INTERIM / "nonidr_exon_overlap_table.parquet"


def count_exon_span(seg, exon_aa_ranges):
    """n coding exons a segment [s,e] (aa) overlaps, given sorted (a_start,a_end)."""
    s, e = seg
    return sum(1 for a0, a1 in exon_aa_ranges if a0 <= e and a1 >= s)


def build():
    idr_by_iso = U.idr_segments_by_isoform()

    val = pd.read_parquet(C.INTERIM / "mapping_validation_scaled.parquet")
    val = val[val["is_validated_mapping"]].copy()
    val["ENST_nover"] = val["ENST"].astype(str).str.split(".").str[0]

    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    exon_aa_by_enst = {
        enst: list(zip(g["aa_start"].astype(int), g["aa_end"].astype(int)))
        for enst, g in cds.groupby("ensembl_transcript_id", sort=False)
    }

    idr_rows, nonidr_rows = [], []
    for r in val.itertuples():
        iso, enst = r.uniprot_isoform, r.ENST_nover
        exon_ranges = exon_aa_by_enst.get(enst)
        if not exon_ranges:
            continue
        prot_len = r.uniprot_isoform_length_aa
        idrs = idr_by_iso.get(iso, [])

        for i, seg in enumerate(idrs, 1):
            n = count_exon_span(seg, exon_ranges)
            seg_len_aa = seg[1] - seg[0] + 1
            seg_len_bp = seg_len_aa * 3
            idr_rows.append({
                "isoform_accession": iso, "base_accession": r.base_accession,
                "ensembl_transcript_id": enst,
                "idr_segment_id": f"{iso}_idr{i}",
                "idr_aa_start": seg[0], "idr_aa_end": seg[1],
                "idr_length_aa": seg_len_aa,
                "n_overlapping_coding_exons": n,
                "contains_exon_boundary": n > 1,
                "n_internal_exon_boundaries": max(n - 1, 0),
                "idr_boundary_density_per_kb": (max(n - 1, 0) / seg_len_bp * 1000) if seg_len_bp else 0.0,
                "idr_encoded_by_single_exon": n == 1,
                "idr_encoded_by_multiple_exons": n > 1,
            })

        for j, seg in enumerate(U.complement_intervals(idrs, prot_len), 1):
            n = count_exon_span(seg, exon_ranges)
            seg_len_aa = seg[1] - seg[0] + 1
            seg_len_bp = seg_len_aa * 3
            nonidr_rows.append({
                "isoform_accession": iso, "base_accession": r.base_accession,
                "ensembl_transcript_id": enst,
                "nonidr_segment_id": f"{iso}_nonidr{j}",
                "nonidr_aa_start": seg[0], "nonidr_aa_end": seg[1],
                "nonidr_length_aa": seg_len_aa,
                "n_overlapping_coding_exons": n,
                "contains_exon_boundary": n > 1,
                "n_internal_exon_boundaries": max(n - 1, 0),
                "nonidr_boundary_density_per_kb": (max(n - 1, 0) / seg_len_bp * 1000) if seg_len_bp else 0.0,
            })

    idr_df    = pd.DataFrame(idr_rows)
    nonidr_df = pd.DataFrame(nonidr_rows)
    idr_df.to_parquet(IDR_OUT, index=False)
    nonidr_df.to_parquet(NONIDR_OUT, index=False)

    print(f"[out] IDR segments      {len(idr_df):,} rows "
          f"({idr_df['isoform_accession'].nunique():,} isoforms) -> {IDR_OUT.name}")
    print(f"[out] non-IDR segments  {len(nonidr_df):,} rows "
          f"({nonidr_df['isoform_accession'].nunique():,} isoforms) -> {NONIDR_OUT.name}")
    print("\nIDR: mean exons/segment = "
          f"{idr_df['n_overlapping_coding_exons'].mean():.2f}, "
          f"multi-exon = {idr_df['idr_encoded_by_multiple_exons'].mean()*100:.1f}%")
    print("IDR boundary density (per kb): median = "
          f"{idr_df['idr_boundary_density_per_kb'].median():.3f}")
    print("non-IDR boundary density (per kb): median = "
          f"{nonidr_df['nonidr_boundary_density_per_kb'].median():.3f}")
    print("\nsample IDR rows:")
    print(idr_df.head(4).to_string(index=False))


if __name__ == "__main__":
    build()
