"""Rebuild every table and figure from the raw inputs, in dependency order.

    python run.py                          # everything (~20 min on a laptop)
    python run.py --steps figures catalog  # one or more steps
    python run.py --list                   # show the steps

Steps
  inputs    check the SHA-256 of every input listed in data/manifest.tsv
  proteome  run notebooks/01-06 top to bottom in fresh kernels -> data/annotated/, figures/proteome*/,
            figures/exon_density/ (go_enrichment_CLUSTER_ONLY needs GO files on the AKEY cluster; not run)
  genomic   pipeline/run_pipeline.py (UniProt->Ensembl mapping, GENCODE coordinates, IDR projection,
            exon tables, splicing analyses) -> data/genomic/, figures/genomic/
  figures   scripts/figures/*.py (figures made outside the notebooks and the pipeline)
  catalog   render figure previews and FIGURES.md from figures/catalog.tsv

Logs go to results/logs/ and executed notebooks to results/executed_notebooks/ (neither tracked by git).
A complete run writes results/run_manifest.json (environment, parameters, timings, output checksums);
a partial run writes the same record to results/logs/partial_run_manifest.json.
PDFs carry a fixed creation date (SOURCE_DATE_EPOCH) so an unchanged figure is byte-identical after a rerun.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata as md
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent
RESULTS = REPO / "results"
LOGS = RESULTS / "logs"
MANIFEST = REPO / "data" / "manifest.tsv"

NOTEBOOKS = [
    "notebooks/01_cohort_and_idrs",               # cohort filters, IDR segments, TF labels -> isoform master tables
    "notebooks/02_proteome_analysis",             # gene age; figures/proteome (fig1-fig7)
    "notebooks/03_appris_annotation",             # APPRIS gene annotation (read by the genomic pipeline)
    "notebooks/04_tfiso_ppi_counts",              # TFIso1.0 isoform PPI counts
    "notebooks/05_proteome_analysis_extended",    # figures/proteome_extended
    "notebooks/06_exon_density_length_matched",   # Susie Song's RefSeq exon-density figures -> figures/exon_density
]
FIGURE_SCRIPTS = [
    # proteome level (read data/annotated)
    "fig2_canonical_vs_alternative_violins.py",
    "fig4A_gene_age_pct_idr_stars.py",
    "fig_family_disorder.py",
    "fig_idr_position_protein_exon.py",      # also reads data/genomic exon_level_table
    "fig4C_5Dalt_length_delta_proteome.py",  # reconstruction
    # genomic (read data/genomic)
    "fig03_idr_modularity.py",
    "fig09_gene_age.py",
    "fig13_gene_grouping_alt_fraction.py",   # reconstruction
    "fig15_idr_gain_loss_tf_vs_nontf.py",
    "fig15b_idr_fate_by_age.py",
    "fig16b_idrfrac_by_age.py",
    "fig_splice_idr_sticker_grammar.py",
    "genomic_preliminary_fig01_to_fig12.py", # reconstruction
    # RefSeq exon density
    "fig_length_matched_tf_vs_nontf_idr_only.py",  # reconstruction
]
OUTPUT_DIRS = ["data/annotated", "data/genomic/analysis", "data/genomic/processed", "data/genomic/isoform_mapping",
               "data/genomic/interim", "data/genomic/qc", "figures", "results/tables"]
PACKAGES = ["numpy", "pandas", "scipy", "matplotlib", "seaborn", "statsmodels", "pyarrow", "openpyxl",
            "biopython", "tqdm", "requests", "nbclient", "ipykernel", "pymupdf"]
ENV = {"PYTHONUTF8": "1", "PYTHONHASHSEED": "0", "MPLBACKEND": "Agg", "SOURCE_DATE_EPOCH": "0"}


def sha256(path: Path) -> str:
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def installed_version(package: str) -> str:
    try:
        return md.version(package)
    except md.PackageNotFoundError:
        return "not installed"


def run_cmd(name: str, cmd: list[str], cwd: Path) -> None:
    print(f"  -> {name}", flush=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    with (LOGS / f"{name}.log").open("w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=cwd, env=os.environ, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc:
        raise SystemExit(f"{name} failed (exit {rc}); see {LOGS / (name + '.log')}")


def step_inputs() -> dict:
    rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8"), delimiter="\t"))
    bad = []
    for r in rows:
        p = REPO / r["path"]
        if not p.exists():
            bad.append(f"MISSING  {r['path']}")
        elif sha256(p) != r["sha256"]:
            bad.append(f"CHANGED  {r['path']}")
    if bad:
        raise SystemExit("Input check failed (see data/README.md for how to obtain the inputs):\n  " + "\n  ".join(bad))
    print(f"  {len(rows)} inputs verified", flush=True)
    return {"verified_inputs": len(rows)}


def step_proteome() -> dict:
    import nbformat
    from nbclient import NotebookClient

    out_dir = RESULTS / "executed_notebooks"
    out_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    for rel in NOTEBOOKS:
        path = REPO / f"{rel}.ipynb"
        print(f"  -> {rel}.ipynb", flush=True)
        nb = nbformat.read(path, as_version=4)
        start = time.time()
        client = NotebookClient(nb, timeout=None, kernel_name="python3", resources={"metadata": {"path": str(path.parent)}})
        try:
            client.execute()
        finally:
            nbformat.write(nb, out_dir / f"{path.stem}.ipynb")
        timings[path.stem] = round(time.time() - start, 1)
    return {"notebook_seconds": timings}


def step_genomic() -> dict:
    start = time.time()
    run_cmd("genomic_pipeline", [sys.executable, "run_pipeline.py"], REPO / "pipeline")
    return {"genomic_seconds": round(time.time() - start, 1)}


def step_figures() -> dict:
    timings = {}
    for script in FIGURE_SCRIPTS:
        start = time.time()
        run_cmd(f"figure_{Path(script).stem}", [sys.executable, str(REPO / "scripts" / "figures" / script)], REPO)
        timings[script] = round(time.time() - start, 1)
    return {"figure_seconds": timings}


def step_catalog() -> dict:
    run_cmd("figure_catalog", [sys.executable, str(REPO / "scripts" / "build_figure_catalog.py")], REPO)
    return {}


STEPS = {"inputs": step_inputs, "proteome": step_proteome, "genomic": step_genomic, "figures": step_figures,
         "catalog": step_catalog}


def output_checksums() -> dict:
    sums = {}
    for d in OUTPUT_DIRS:
        for p in sorted((REPO / d).rglob("*")):
            if p.is_file():
                sums[p.relative_to(REPO).as_posix()] = sha256(p)
    return sums


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--steps", nargs="+", choices=list(STEPS), default=list(STEPS))
    parser.add_argument("--list", action="store_true", help="print the step descriptions and exit")
    args = parser.parse_args()
    if args.list:
        print(__doc__)
        return

    os.environ.update(ENV)   # inherited by notebook kernels and subprocesses
    min_idr = os.environ.get("TFIDR_MIN_IDR_LEN_AA", "20")
    record = {
        "started": dt.datetime.now().isoformat(timespec="seconds"),
        "steps": args.steps,
        "parameters": {"TFIDR_MIN_IDR_LEN_AA": min_idr},
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {p: installed_version(p) for p in PACKAGES},
    }
    if min_idr != "20":
        print(f"WARNING: TFIDR_MIN_IDR_LEN_AA={min_idr} (not the paper's IDR definition)", flush=True)
    for name in args.steps:
        print(f"[{name}]", flush=True)
        record[name] = STEPS[name]()
    record["finished"] = dt.datetime.now().isoformat(timespec="seconds")
    record["outputs_sha256"] = output_checksums()
    # the tracked manifest describes a complete rebuild; partial runs are logged next to their step logs
    target = RESULTS / "run_manifest.json" if args.steps == list(STEPS) else LOGS / "partial_run_manifest.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    print(f"done -> {target}")


if __name__ == "__main__":
    main()
