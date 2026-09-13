"""Run the TF-IDR genomic-coordinate pipeline in dependency order.

Order (each step reads only files written by earlier steps or by the notebooks in ../notebooks):

  fetch    UniProt -> Ensembl cross-references      (data/raw/uniprot; live UniProt API if absent)
  map      candidate isoform -> transcript mapping
  s1       GENCODE downloads                          (data/raw/gencode; downloaded if absent)
  s2..s9   coordinates, validation, IDR overlap, projection, master/gene tables, QC, splicing analysis
  resolve  strict transcript resolution
  events   splice-event classification
  exons    exon-level export table
  igv      IGV BED tracks
  s9b      alternative-exon IDR enrichment            -> figures/genomic/fig14
  s9c      continuous exon IDR fraction               -> fig16
  s9d      gene-level IDR fraction                    -> fig17
  s10b     IDR exon fate (conserved/removed/added)    -> fig15

Examples:
    python pipeline/run_pipeline.py                     # everything
    python pipeline/run_pipeline.py --from s4 --through s9
    python pipeline/run_pipeline.py --from resolve      # post-processing and downstream analyses only
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

STEPS = {
    "fetch": "tfidr_pipeline.mapping.fetch_uniprot_ensembl",
    "map": "tfidr_pipeline.mapping.build_mapping_table",
    "s1": "tfidr_pipeline.stages.s1_download_gencode",
    "s2": "tfidr_pipeline.stages.s2_exon_cds_coords",
    "s3": "tfidr_pipeline.stages.s3_validate_translations",
    "s4": "tfidr_pipeline.stages.s4_idr_exon_overlap",
    "s4b": "tfidr_pipeline.stages.s4b_genomic_projection",
    "s5": "tfidr_pipeline.stages.s5_build_master",
    "s6": "tfidr_pipeline.stages.s6_gene_summary",
    "s7": "tfidr_pipeline.stages.s7_viz_tables",
    "s8": "tfidr_pipeline.stages.s8_qc_and_subsets",
    "s9": "tfidr_pipeline.stages.s9_splicing_idr_analysis",
    "resolve": "tfidr_pipeline.mapping.resolve_isoform_transcripts",
    "events": "tfidr_pipeline.mapping.classify_splice_events",
    "exons": "tfidr_pipeline.exports.exon_level_table",
    "igv": "tfidr_pipeline.exports.make_igv_tracks",
    "s9b": "tfidr_pipeline.stages.s9b_idr_enrichment",
    "s9c": "tfidr_pipeline.stages.s9c_exon_idr_fraction",
    "s9d": "tfidr_pipeline.stages.s9d_gene_idr_fraction",
    "s10b": "tfidr_pipeline.stages.s10b_idr_gain_loss",
}


def select(start: str, end: str) -> list[tuple[str, str]]:
    names = list(STEPS)
    first, last = names.index(start), names.index(end)
    if first > last:
        raise SystemExit("--from must come before --through")
    return [(n, STEPS[n]) for n in names[first:last + 1]]


def run(steps: list[tuple[str, str]]) -> None:
    for name, module in steps:
        print(f"\n{'=' * 70}\n[{name}] {module}\n{'=' * 70}", flush=True)
        started = time.time()
        subprocess.run([sys.executable, "-m", module], cwd=HERE, check=True)
        print(f"[{name}] completed in {time.time() - started:.0f}s", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from", dest="start", choices=STEPS, default="fetch")
    parser.add_argument("--through", choices=STEPS, default="s10b")
    args = parser.parse_args()
    run(select(args.start, args.through))
    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
