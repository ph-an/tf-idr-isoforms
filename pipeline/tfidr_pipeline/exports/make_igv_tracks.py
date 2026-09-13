"""
Build IGV-ready BED12 tracks (GENCODE v50 / GRCh38 / hg38) from the s4b genomic
projection, so each IDR can be viewed as ONE feature — with its exon-split pieces
drawn as blocks — against the isoform's coding-exon structure.

Load in IGV with the **hg38** genome. Feature names are `GENE|ISOFORM|idrN`, so
IGV's search box finds e.g. "TP53". IDRs are colored by TF status.

Outputs (processed/igv/):
    idr_segments.hg38.bed        BED12 — one feature per IDR segment (exon-blocked)
    nonidr_segments.hg38.bed     BED12 — one feature per non-IDR segment
    cds_isoforms.hg38.bed        BED12 — one feature per isoform (coding-exon structure)

Usage:
    python make_igv_tracks.py                # all validated isoforms
    python make_igv_tracks.py P04637 TP53    # only these accessions / gene symbols
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
from tfidr_pipeline import config as C

OUTDIR = C.PROCESSED / "igv"
OUTDIR.mkdir(exist_ok=True)

TF_RGB, NONTF_RGB, CDS_RGB = "213,94,0", "0,114,178", "120,120,120"  # Okabe-Ito


def blocks_from_pieces(pieces):
    """pieces: list of (start0, end) 0-based half-open, disjoint. -> BED12 fields."""
    pieces = sorted(pieces)
    chrom_start = pieces[0][0]
    chrom_end   = pieces[-1][1]
    sizes  = [e - s for s, e in pieces]
    starts = [s - chrom_start for s, e in pieces]
    return (chrom_start, chrom_end, len(pieces),
            ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ",")


def write_bed12(rows, path, track_name, desc):
    rows = sorted(rows, key=lambda r: (r[0], r[1]))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f'track name="{track_name}" description="{desc}" '
                 f'itemRgb="On" visibility="pack"\n')
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return len(rows)


def seg_track(parquet, seg_id_col, label, rgb_map, gene_of, path, track_name):
    df = pd.read_parquet(parquet)
    df = df[df["projection_validated"]] if "projection_validated" in df else df
    rows = []
    for sid, g in df.groupby(seg_id_col):
        r0 = g.iloc[0]
        iso = r0["isoform_accession"]
        cs, ce, n, sizes, starts = blocks_from_pieces(
            list(zip(g["genomic_start_0based"].astype(int), g["genomic_end"].astype(int))))
        name = f"{gene_of.get(iso,'?')}|{iso}|{sid.split('_')[-1]}"
        rgb = rgb_map.get(iso, NONTF_RGB)
        rows.append((r0["chromosome"], cs, ce, name, 0, r0["strand"],
                     cs, ce, rgb, n, sizes, starts))
    return write_bed12(rows, path, track_name,
                       f"{label} segments, GENCODE v50/hg38, colored by TF status")


def main(filt=None):
    m = pd.read_parquet(C.PROCESSED / "tf_idr_isoform_master.parquet")
    m = m[m.genomic_scaffold_validated == True]
    gene_of = dict(zip(m.isoform_accession, m.gene_name))
    tf_of   = dict(zip(m.isoform_accession, m.is_tf))
    rgb_map = {i: (TF_RGB if tf_of.get(i) else NONTF_RGB) for i in m.isoform_accession}

    # optional subset by accession or gene symbol
    keep = None
    if filt:
        fset = set(filt)
        keep = set(m[m.isoform_accession.isin(fset) | m.base_accession.isin(fset)
                     | m.gene_name.isin(fset)].isoform_accession)
        print(f"[filter] {len(keep):,} isoforms match {filt}")

    def maybe_filter(pq):
        df = pd.read_parquet(pq)
        if keep is not None:
            df = df[df.isoform_accession.isin(keep)]
        return df

    # IDR + non-IDR segment tracks
    for pq, sidcol, lab, fname, tname in [
        (C.INTERIM / "idr_genomic_intervals.parquet", "idr_segment_id", "IDR",
         "idr_segments.hg38.bed", "IDR segments"),
        (C.INTERIM / "nonidr_genomic_intervals.parquet", "nonidr_segment_id", "non-IDR",
         "nonidr_segments.hg38.bed", "non-IDR segments")]:
        df = maybe_filter(pq)
        df = df[df["projection_validated"]] if "projection_validated" in df else df
        rows = []
        for sid, g in df.groupby(sidcol):
            r0 = g.iloc[0]; iso = r0["isoform_accession"]
            cs, ce, n, sizes, starts = blocks_from_pieces(
                list(zip(g["genomic_start_0based"].astype(int), g["genomic_end"].astype(int))))
            name = f"{gene_of.get(iso,'?')}|{iso}|{str(sid).split('_')[-1]}"
            rows.append((r0["chromosome"], cs, ce, name, 0, r0["strand"],
                         cs, ce, rgb_map.get(iso, NONTF_RGB), n, sizes, starts))
        nlab = write_bed12(rows, OUTDIR / fname, tname,
                           f"{lab} segments, GENCODE v50/hg38, colored by TF (red) / non-TF (blue)")
        print(f"[out] {fname:28} {nlab:,} features")

    # CDS-structure track (one feature per isoform = its coding exons as blocks)
    ex = pd.read_parquet(C.INTERIM / "ensembl_transcript_exon_coordinates.parquet")
    ex = ex[ex.is_coding_exon].copy()
    m2 = m.copy(); m2["ENST"] = m2.resolved_transcript.astype(str).str.split(".").str[0]
    iso_enst = m2[["isoform_accession", "ENST", "gene_name"]]
    if keep is not None:
        iso_enst = iso_enst[iso_enst.isoform_accession.isin(keep)]
    cds_by = {e: g for e, g in ex.groupby("ensembl_transcript_id")}
    rows = []
    for r in iso_enst.itertuples():
        g = cds_by.get(r.ENST)
        if g is None: continue
        pieces = list(zip((g.cds_genomic_start.astype(int) - 1), g.cds_genomic_end.astype(int)))
        cs, ce, n, sizes, starts = blocks_from_pieces(pieces)
        name = f"{r.gene_name}|{r.isoform_accession}|{r.ENST}"
        rows.append((g.chromosome.iloc[0], cs, ce, name, 0, g.strand.iloc[0],
                     cs, ce, CDS_RGB, n, sizes, starts))
    ncds = write_bed12(rows, OUTDIR / "cds_isoforms.hg38.bed",
                       "isoform CDS structure", "coding-exon structure per isoform, GENCODE v50/hg38")
    print(f"[out] cds_isoforms.hg38.bed        {ncds:,} features")
    print(f"\nLoad in IGV with genome = hg38 (GRCh38). Search by gene symbol, e.g. TP53.")
    print(f"Files -> {OUTDIR}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
