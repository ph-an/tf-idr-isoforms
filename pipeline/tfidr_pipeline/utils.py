"""Shared data-loading and interval helpers for the pipeline."""
import pandas as pd
from tfidr_pipeline import config as C


def load_idr_frame():
    """Load MetaPredict rows, normalize the UniProt accession to ``acc``, and keep IDRs
    of at least ``config.MIN_IDR_LEN_AA`` residues (the definition used for pct_idr)."""
    df = pd.read_csv(C.METAPREDICT_IDROME)
    df.columns = df.columns.str.strip()
    df["acc"] = df["FASTA header"].str.strip().str.extract(r"sp\|([^|]+)\|")[0]
    df = df.dropna(subset=["IDR start", "IDR end"])
    return df[df["IDR len"] >= C.MIN_IDR_LEN_AA]


def idr_segments_by_isoform(df=None):
    """Return ``{isoform: [(aa_start, aa_end), ...]}``."""
    df = load_idr_frame() if df is None else df
    return {
        acc: sorted(zip(g["IDR start"].astype(int), g["IDR end"].astype(int)))
        for acc, g in df.groupby("acc")
    }


def complement_intervals(intervals, length):
    """Inclusive intervals in ``[1, length]`` not covered by ``intervals``."""
    if length is None or length < 1:
        return []
    gaps, cursor = [], 1
    for start, end in sorted(intervals):
        if start > cursor:
            gaps.append((cursor, start - 1))
        cursor = max(cursor, end + 1)
    if cursor <= length:
        gaps.append((cursor, length))
    return gaps


def overlap_length(start, end, intervals):
    """Number of positions in inclusive ``[start, end]`` covered by intervals."""
    return sum(
        max(0, min(end, right) - max(start, left) + 1)
        for left, right in intervals
    )


def intervals_overlap(a, b):
    return a[0] <= b[1] and b[0] <= a[1]
