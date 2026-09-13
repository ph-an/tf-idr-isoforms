"""
Stage D.3 (s4b) — project IDR / non-IDR segments onto real GENOMIC coordinates.

This finishes the step the pipeline previously deferred (METHODOLOGY.md had the
aa->CDS->genomic projection "on hold"). It ports the approach from
Susie Song's 01_genomic_coords.ipynb and adapts it to our local-GENCODE, isoform-resolved
architecture, adding a DNA round-trip that SELF-VERIFIES every projection — the
key correctness guarantee Susie has and s3 (protein-only) does not.

For each sequence-validated isoform and its RESOLVED ENST (s3):

  1. Build a strand-aware flat array of genomic positions over the CDS, walking
     coding exons in ascending GENCODE exon_rank (5'->3'):
         + strand:  arange(cds_start, cds_end+1,  +1)
         - strand:  arange(cds_end,   cds_start-1, -1)
     concatenated in rank order -> one genomic position per CDS nucleotide.

  2. Pull the ENST's CDS nucleotides from GENCODE pc_transcripts (header CDS:a-b),
     aligned to the flat array. Scaffold check: translate(cds_nt) reproduces the
     protein exactly and len(cds_nt) == len(flat) -> genomic_scaffold_validated. Because
     the flat positions and cds_nt come from the same GENCODE release, cds_nt[i]
     is the nucleotide at genomic position flat[i].

  3. For every IDR segment (MetaPredict) and its non-IDR complement:
         nt slice = [(aa_start-1)*3 : aa_end*3]
         genomic intervals = array_to_intervals(flat[nt])     # collapses runs
         projection_validated = translate(cds_nt[nt]) == protein[aa_start-1:aa_end]
     The exon-piece count is the number of ACTUAL genomic interval breaks
     (n_pieces-1 boundaries), which resolves the split-codon +/-1 ambiguity of the
     aa-space estimate in s4.

Outputs:
    interim/idr_genomic_intervals.parquet       one row per IDR genomic piece (BED-like)
    interim/nonidr_genomic_intervals.parquet    one row per non-IDR genomic piece
    interim/genomic_projection_segments.parquet one row per segment (IDR & non-IDR)
    interim/genomic_scaffold_validation.parquet one row per isoform
    processed/idrs.bed  processed/nonidrs.bed  processed/exons.bed   BED6, 0-based, sorted
    qc/genomic_projection_failures.csv
"""
import sys, gzip, re
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from tfidr_pipeline import config as C
from tfidr_pipeline import utils as U

IDR_OUT      = C.INTERIM / "idr_genomic_intervals.parquet"
NONIDR_OUT   = C.INTERIM / "nonidr_genomic_intervals.parquet"
SEG_OUT      = C.INTERIM / "genomic_projection_segments.parquet"
SCAFFOLD_OUT = C.INTERIM / "genomic_scaffold_validation.parquet"
IDR_BED      = C.PROCESSED / "idrs.bed"
NONIDR_BED   = C.PROCESSED / "nonidrs.bed"
EXON_BED     = C.PROCESSED / "exons.bed"
FAIL_CSV     = C.QC / "genomic_projection_failures.csv"

# ── standard genetic code (NCBI table 1); no biopython dep in the s-stages ──────
_BASES = "TCAG"
_AAS   = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AAS[i]
               for i, (a, b, c) in enumerate(
                   (x, y, z) for x in _BASES for y in _BASES for z in _BASES)}


def translate(nt):
    """Translate a nucleotide string (multiple of 3 assumed) with table 1.
    Unknown/ambiguous codons -> 'X'; stop codons -> '*'."""
    nt = nt.upper()
    n = len(nt) - len(nt) % 3
    return "".join(CODON_TABLE.get(nt[i:i + 3], "X") for i in range(0, n, 3))


def array_to_intervals(arr):
    """Collapse a run of consecutive genomic positions (step +1 or -1) into
    (start, end) intervals. Ported from Susie Song's 01_genomic_coords.ipynb; a break in
    the +/-1 pattern (an exon junction) starts a new interval."""
    if len(arr) == 0:
        return []
    arr = np.asarray(arr)
    diffs = np.diff(arr)
    breaks = np.where((diffs != 1) & (diffs != -1))[0]
    intervals, start_idx = [], 0
    for b in breaks:
        intervals.append((int(arr[start_idx]), int(arr[b])))
        start_idx = b + 1
    intervals.append((int(arr[start_idx]), int(arr[-1])))
    return intervals


def load_pc_transcripts(wanted):
    """{ENST_unversioned: (transcript_seq, cds_start_1based, cds_end_1based)} for
    the transcripts we need. Header: >ENST..|..|CDS:start-end|.."""
    out, enst, cds, buf = {}, None, None, []
    cds_re = re.compile(r"CDS:(\d+)-(\d+)")
    with gzip.open(C.GENCODE_PC_TRANSCRIPTS_GZ, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if enst is not None and enst in wanted and cds is not None:
                    out[enst] = ("".join(buf), cds[0], cds[1])
                parts = line[1:].split("|")
                enst = next((p.split(".")[0] for p in parts if p.startswith("ENST")), None)
                m = next((cds_re.match(p) for p in parts if p.startswith("CDS:")), None)
                cds = (int(m.group(1)), int(m.group(2))) if m else None
                buf = []
            else:
                buf.append(line.strip())
        if enst is not None and enst in wanted and cds is not None:
            out[enst] = ("".join(buf), cds[0], cds[1])
    return out


def build_flat_positions(cds_rows):
    """Flat 1-based genomic position per CDS nucleotide, 5'->3', strand-aware.
    cds_rows: iterable of (cds_genomic_start, cds_genomic_end, strand) sorted by
    ascending exon_rank."""
    pieces = []
    for start, end, strand in cds_rows:
        start, end = int(start), int(end)
        if strand == "+":
            pieces.append(np.arange(start, end + 1, 1))
        else:                                   # '-': 5' end is the high coordinate
            pieces.append(np.arange(end, start - 1, -1))
    return np.concatenate(pieces) if pieces else np.array([], dtype=int)


def intervals_to_bed(intervals):
    """(start,end) 1-based inclusive genomic intervals -> list of (start0, end)
    BED half-open, min/max normalized (strand-agnostic coordinates)."""
    out = []
    for a, b in intervals:
        lo, hi = (a, b) if a <= b else (b, a)
        out.append((lo - 1, hi))               # 0-based start, exclusive end
    return out


def main():
    idr_by_iso = U.idr_segments_by_isoform()

    # sequence-validated isoforms + their resolved transcript (s3)
    val = pd.read_parquet(C.INTERIM / "mapping_validation_scaled.parquet")
    if "resolved_transcript" not in val.columns:
        sys.exit("mapping_validation_scaled.parquet lacks resolved_transcript — rerun s3.")
    val = val[val["resolved_is_validated"].fillna(False)].copy()
    val["ENST"] = val["resolved_transcript"].astype(str).str.split(".").str[0]

    # reference protein sequences (UniProt isoform) for the round-trip assert
    uni_seq = (pd.read_csv(C.ISOFORM_ANNOTATED)
               .set_index("isoform_accession")["Sequence"].to_dict())
    base_of = dict(zip(val["uniprot_isoform"], val["base_accession"]))

    # CDS genomic coords per transcript (from s2), grouped + rank-sorted
    cds = pd.read_parquet(C.INTERIM / "transcript_cds_to_protein_coordinate_map.parquet")
    cds = cds.sort_values(["ensembl_transcript_id", "exon_rank"])
    cds_by_enst = {
        enst: (g["chromosome"].iloc[0], g["strand"].iloc[0],
               list(zip(g["cds_genomic_start"], g["cds_genomic_end"], g["strand"])))
        for enst, g in cds.groupby("ensembl_transcript_id", sort=False)
    }

    wanted_enst = set(val["ENST"]) & set(cds_by_enst)
    print(f"[load] pc_transcripts for {len(wanted_enst):,} resolved transcripts ...")
    tx = load_pc_transcripts(wanted_enst)
    print(f"[load] {len(tx):,} transcript CDS sequences")

    idr_rows, nonidr_rows, seg_rows, scaffold_rows, exon_rows, fails = [], [], [], [], [], []

    for r in val.itertuples():
        iso, enst = r.uniprot_isoform, r.ENST
        meta = cds_by_enst.get(enst)
        prot = uni_seq.get(iso)
        if meta is None or not prot or enst not in tx:
            scaffold_rows.append({"isoform_accession": iso, "ensembl_transcript_id": enst,
                                  "genomic_scaffold_validated": False,
                                  "reason": "missing_cds_or_seq_or_transcript"})
            continue
        chrom, strand, cds_rows = meta
        flat = build_flat_positions(cds_rows)
        seq, cds_a, cds_b = tx[enst]
        cds_nt = seq[cds_a - 1:cds_b]                 # includes stop codon
        cds_nt = cds_nt[:len(flat)]                   # align to GTF CDS (drops stop)

        # ── scaffold self-check: does the CDS translate back to the protein? ──
        scaffold_ok = (len(cds_nt) == len(flat) and len(flat) % 3 == 0)
        prot_from_cds = translate(cds_nt) if scaffold_ok else ""
        # Conservative policy: the CDS translation must be fully identical.
        scaffold_ok = scaffold_ok and prot_from_cds == prot
        scaffold_rows.append({"isoform_accession": iso, "ensembl_transcript_id": enst,
                              "chromosome": chrom, "strand": strand,
                              "cds_len_bp": int(len(flat)),
                              "protein_len_aa": len(prot),
                              "cds_len_coherent": len(flat) == 3 * len(prot),
                              "genomic_scaffold_validated": bool(scaffold_ok)})
        if not scaffold_ok:
            fails.append({"isoform_accession": iso, "ensembl_transcript_id": enst,
                          "level": "scaffold", "detail": "cds_translation_mismatch"})
            continue

        # coding-exon (CDS) intervals for exons.bed
        for a, b in intervals_to_bed([(int(s), int(e)) for s, e, _ in cds_rows]):
            exon_rows.append((chrom, a, b, iso, 0, strand))

        base = base_of.get(iso, r.base_accession)
        idrs = idr_by_iso.get(iso, [])

        def project(seg, seg_type, seg_id, rows_sink):
            s, e = seg
            nt0, nt1 = (s - 1) * 3, e * 3
            if nt1 > len(flat):                       # segment beyond CDS
                fails.append({"isoform_accession": iso, "ensembl_transcript_id": enst,
                              "level": seg_type, "detail": f"segment_{s}_{e}_beyond_cds"})
                return None
            gpos = flat[nt0:nt1]
            intervals = array_to_intervals(gpos)
            aa_obs = translate(cds_nt[nt0:nt1])
            aa_ref = prot[s - 1:e]
            ok = (aa_obs == aa_ref)
            if not ok:
                fails.append({"isoform_accession": iso, "ensembl_transcript_id": enst,
                              "level": seg_type, "detail": f"segment_{s}_{e}_translation_mismatch"})
            n_pieces = len(intervals)
            seg_len_bp = (e - s + 1) * 3
            for k, bed in enumerate(intervals_to_bed(intervals), 1):
                rows_sink.append({
                    "isoform_accession": iso, "base_accession": base,
                    "ensembl_transcript_id": enst, "chromosome": chrom, "strand": strand,
                    f"{seg_type}_segment_id": seg_id,
                    "aa_start": s, "aa_end": e,
                    "genomic_start_0based": bed[0], "genomic_end": bed[1],
                    "piece_index": k, "n_pieces": n_pieces,
                    "projection_validated": ok,
                })
            seg_rows.append({
                "isoform_accession": iso, "base_accession": base,
                "ensembl_transcript_id": enst, "chromosome": chrom, "strand": strand,
                "segment_type": seg_type, "segment_id": seg_id,
                "aa_start": s, "aa_end": e, "length_aa": e - s + 1,
                "n_genomic_pieces": n_pieces,
                "n_internal_exon_boundaries": max(n_pieces - 1, 0),
                "contains_exon_boundary": n_pieces > 1,
                "boundary_density_per_kb": (max(n_pieces - 1, 0) / seg_len_bp * 1000) if seg_len_bp else 0.0,
                "projection_validated": ok,
            })
            return ok

        for i, seg in enumerate(idrs, 1):
            project(seg, "idr", f"{iso}_idr{i}", idr_rows)
        for j, seg in enumerate(U.complement_intervals(idrs, len(prot)), 1):
            project(seg, "nonidr", f"{iso}_nonidr{j}", nonidr_rows)

    # ── assemble + write ──────────────────────────────────────────────────────
    idr_df    = C.stamp(pd.DataFrame(idr_rows))
    nonidr_df = C.stamp(pd.DataFrame(nonidr_rows))
    seg_df    = C.stamp(pd.DataFrame(seg_rows))
    scaf_df   = C.stamp(pd.DataFrame(scaffold_rows))
    idr_df.to_parquet(IDR_OUT, index=False)
    nonidr_df.to_parquet(NONIDR_OUT, index=False)
    seg_df.to_parquet(SEG_OUT, index=False)
    scaf_df.to_parquet(SCAFFOLD_OUT, index=False)

    def write_bed(rows, path):
        b = pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "score", "strand"])
        b = b.sort_values(["chrom", "start", "end"])
        b.to_csv(path, sep="\t", header=False, index=False)
        return len(b)

    n_idr_bed = write_bed(
        [(d["chromosome"], d["genomic_start_0based"], d["genomic_end"],
          d["isoform_accession"], 0, d["strand"]) for d in idr_rows], IDR_BED)
    n_nonidr_bed = write_bed(
        [(d["chromosome"], d["genomic_start_0based"], d["genomic_end"],
          d["isoform_accession"], 0, d["strand"]) for d in nonidr_rows], NONIDR_BED)
    n_exon_bed = write_bed(exon_rows, EXON_BED)

    pd.DataFrame(fails).to_csv(FAIL_CSV, index=False)

    # ── report ────────────────────────────────────────────────────────────────
    n_iso = scaf_df["isoform_accession"].nunique()
    n_scaf_ok = int(scaf_df["genomic_scaffold_validated"].sum())
    print(f"\n[out] scaffold: {n_scaf_ok:,}/{n_iso:,} isoforms genomic-scaffold-validated "
          f"({n_scaf_ok/max(n_iso,1)*100:.1f}%)")
    if len(seg_df):
        pv = seg_df["projection_validated"]
        idr_seg = seg_df[seg_df["segment_type"] == "idr"]
        non_seg = seg_df[seg_df["segment_type"] == "nonidr"]
        print(f"[out] segments: {len(seg_df):,} projected, "
              f"{pv.mean()*100:.1f}% round-trip-validated")
        print(f"[out] IDR genomic pieces -> {IDR_OUT.name} ({len(idr_df):,} rows, {n_idr_bed:,} BED)")
        print(f"[out] non-IDR pieces     -> {NONIDR_OUT.name} ({len(nonidr_df):,} rows, {n_nonidr_bed:,} BED)")
        print(f"[out] exons.bed          -> {n_exon_bed:,} coding-exon rows")
        print(f"\nIDR median boundary density/kb (genomic):    {idr_seg['boundary_density_per_kb'].median():.3f}")
        print(f"non-IDR median boundary density/kb (genomic): {non_seg['boundary_density_per_kb'].median():.3f}")
        print(f"IDR mean genomic pieces/segment: {idr_seg['n_genomic_pieces'].mean():.2f}, "
              f"multi-exon = {(idr_seg['n_genomic_pieces'] > 1).mean()*100:.1f}%")
    print(f"[qc]  projection failures logged: {len(fails):,} -> {FAIL_CSV.name}")


if __name__ == "__main__":
    main()
