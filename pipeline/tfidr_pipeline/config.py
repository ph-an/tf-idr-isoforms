"""
Central configuration for the TF-IDR isoform genomic-coordinate pipeline.

Everything downstream imports paths, release versions, and binning schemes from
here so releases are configurable (never hard-coded) and reruns are reproducible.

Design decision (agreed with user): REUSE existing validated annotation rather
than rebuild from raw. data/annotated/ already carries TF family, gene age,
APPRIS, and IDR features; data/genomic/ carries the sequence-
validated UniProt->Ensembl isoform mapping. This pipeline fills the genuine
gaps: exon/CDS coordinates at scale, IDR<->exon overlap, and the assembled
master / gene / analysis / QC tables.
"""
from pathlib import Path

# ── Roots ─────────────────────────────────────────────────────────────────────
# Layout: pipeline/tfidr_pipeline/config.py under the repository root
PACKAGE   = Path(__file__).resolve().parent           # pipeline/tfidr_pipeline/
ROOT      = PACKAGE.parent                            # pipeline/
PROJECT   = ROOT.parent                               # repository root
DATA      = PROJECT / "data"
RAW       = DATA / "raw"                              # external source data (read-only)

# ── Pipeline I/O (all under data/genomic/) ────────────────────────────────────
OUT       = DATA / "genomic"
CACHE     = OUT / "cache"                             # parsed download caches
INTERIM   = OUT / "interim"
PROCESSED = OUT / "processed"
ANALYSIS  = OUT / "analysis"                       # subsets, viz tables, statistics
QC        = OUT / "qc"
ISOFORM_MAPPING = OUT / "isoform_mapping"
FIGURES   = PROJECT / "figures" / "genomic"
GENCODE_DIR = RAW / "gencode"                         # s1 downloads here if the files are absent
UNIPROT_XREF_TSV = RAW / "uniprot" / "uniprot_ensembl_xref_raw.tsv"   # fetch downloads it if absent
for _d in (CACHE, INTERIM, PROCESSED, ANALYSIS, QC, ISOFORM_MAPPING, FIGURES, GENCODE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ── Reused, already-validated inputs (external, under tphan/data/) ─────────────
# Isoform backbone + features (TF family, gene age, age bin, pct_idr, n segments)
ISOFORM_ANNOTATED   = DATA / "annotated" / "IDRisoforms_df_geneage.csv"
# Gene-level APPRIS (ENSG-level: appris_annotation, appris_tier, ...)
ISOFORM_APPRIS      = DATA / "annotated" / "IDRisoforms_df_appris.csv"
# Per-isoform, per-segment MetaPredict IDR intervals (aa coords + sequence features)
METAPREDICT_IDROME  = DATA / "raw" / "metapredict" / "IDRome_AllSwissProtHumansProIsos.csv"
# Sequence-validated UniProt isoform -> Ensembl transcript mapping (Stages A-C)
ISOFORM_MAP         = OUT / "isoform_transcript_map.csv"
ISOFORM_CANDIDATES  = OUT / "isoform_transcript_candidates_long.csv"

# ── GENCODE (configurable release) ────────────────────────────────────────────
GENCODE_RELEASE = "50"                    # GRCh38; bump freely
GENCODE_BASE = (f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
                f"release_{GENCODE_RELEASE}")
# Comprehensive annotation (not `.basic`): recovers sequence-valid ENSTs the basic
# set omits (the former `valid_no_cds` group). Extra low-support transcripts it adds
# (NMD/readthrough/partial) are filtered downstream by exact-match + s4b round-trip.
GENCODE_GTF_URL          = f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.annotation.gtf.gz"
GENCODE_TRANSLATIONS_URL = f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.pc_translations.fa.gz"
# pc_transcripts = spliced transcript nucleotides + CDS:start-end header offsets.
# Needed by s4b for the DNA round-trip that self-verifies each genomic projection.
GENCODE_PC_TRANSCRIPTS_URL = f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.pc_transcripts.fa.gz"
GENCODE_GTF_GZ           = GENCODE_DIR / f"gencode.v{GENCODE_RELEASE}.annotation.gtf.gz"
GENCODE_TRANSLATIONS_GZ  = GENCODE_DIR / f"gencode.v{GENCODE_RELEASE}.pc_translations.fa.gz"
GENCODE_PC_TRANSCRIPTS_GZ  = GENCODE_DIR / f"gencode.v{GENCODE_RELEASE}.pc_transcripts.fa.gz"
GENOME_BUILD = "GRCh38"

# ── GTEx (deferred module — configurable, not merged until mappings are clean) ─
GTEX_VERSION = "V10"                      # V10=GENCODE39, V11=GENCODE47; set when used
GTEX_ENABLED = False

# ── Binning schemes ───────────────────────────────────────────────────────────
AGE_BINS_MYA        = [0, 100, 500, 1000, float("inf")]
AGE_BIN_LABELS      = ["<100", "100-500", "500-1000", ">1000"]

ISOFORM_COUNT_BINS  = [3, 4, 5, 6, 7, 8, 9, 10]     # 10 => "10+"
SEGMENT_LEN_AA_BINS       = [0, 50, 100, 200, 300, 400, float("inf")]
SEGMENT_LEN_AA_LABELS     = ["0-50", "50-100", "100-200", "200-300", "300-400", "400+"]

# ── Analysis thresholds ───────────────────────────────────────────────────────
# Minimum MetaPredict IDR length (aa). Must match notebooks/01_cohort_and_idrs.ipynb, which
# drops shorter segments before computing pct_idr / idr_aaLen for the isoform master tables.
# Before 2026-09 this pipeline used every segment (no minimum), so protein-level and genomic
# IDR calls disagreed for 4,274 isoforms. Set TFIDR_MIN_IDR_LEN_AA=0 only to reproduce those
# legacy outputs (e.g. the 2026-07 poster numbers).
import os as _os
MIN_IDR_LEN_AA = int(_os.environ.get("TFIDR_MIN_IDR_LEN_AA", "20"))
MIN_ISOFORMS_FOR_GENE_DELTA = 3      # genes need >=3 isoforms for meaningful dIDR
LONG_IDR_THRESHOLDS_AA      = [30, 50, 100]

# Conservative validation policy. Other match classes are retained for QC but
# are not eligible for transcript-, exon-, or genomic-coordinate analyses.
VALID_SEQUENCE_MATCH_TYPES = ("exact",)

# ── Release provenance stamped into every output ──────────────────────────────
import datetime as _dt
RELEASE_META = {
    "data_release_uniprot":  "2025_reviewed_human",   # from step1 stream query
    "data_release_ensembl":  f"GENCODE_v{GENCODE_RELEASE}_{GENOME_BUILD}",
    "data_release_appris":   "reused_from_annotated",
    "data_release_gtex":     GTEX_VERSION if GTEX_ENABLED else "not_included",
    "pipeline_run_date":     _dt.date.today().isoformat(),
}


def stamp(df):
    """Attach release-provenance columns to any output frame."""
    for k, v in RELEASE_META.items():
        df[k] = v
    return df
