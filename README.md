# Alternative splicing and intrinsic disorder in human transcription-factor isoforms

Tommy Phan, Susie Song, Joshua Akey — Lewis-Sigler Institute for Integrative Genomics, Princeton University

**Question.** How are intrinsically disordered regions (IDRs) distributed across protein isoforms and coding
exons? Does alternative splicing add or remove them? And is that different in transcription factors (TFs)?

**Approach.** Reviewed human UniProt isoforms (genes with ≥3 isoforms) are annotated with MetaPredict IDRs,
HumanTFs TF status, GenOrigin gene age and APPRIS. Each isoform is then matched to a sequence-identical
GENCODE v50 transcript. That lets its IDRs be placed on coding exons and genome coordinates, and
alternative and constitutive exons be compared.

> **Status:** pre-publication. The repository holds every figure produced so far; the ones that go into the
> paper are still being chosen. Start with **[FIGURES.md](FIGURES.md)**.

---

## Choosing figures for the paper

* **[FIGURES.md](FIGURES.md)** is a gallery of all 75 figures, grouped by question. Each has a preview, what it
  shows, its statistics, the code that makes it, and caveats (duplicates, reconstructed code, known issues).
* **[figures/catalog.tsv](figures/catalog.tsv)** is the same list as a spreadsheet. Fill in the `paper` column
  (e.g. `Fig 2B`, `Supp 3`, `drop`), then run `python run.py --steps catalog` to update the gallery.
* PDFs are in `figures/<topic>/`; `figures/previews/` holds PNG renders for the gallery.

## Quick start

```bash
# 1. environment (Python 3.13)
python -m venv .venv
.venv\Scripts\activate                 # macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt

# 2. inputs: put the files listed in data/README.md under data/raw/ (checked against data/manifest.tsv)

# 3. rebuild every table and figure from the raw inputs (~25 min on a laptop)
python run.py
```

`run.py` stops at the first failing step and names its log in `results/logs/`. Steps can be run separately:

```bash
python run.py --steps inputs     # check input checksums
python run.py --steps proteome   # notebooks 01-06: cohort, IDRs, gene age, APPRIS, PPI, proteome figures
python run.py --steps genomic    # genomic-coordinate pipeline and splicing analyses
python run.py --steps figures    # figures drawn by standalone scripts
python run.py --steps catalog    # figure previews and FIGURES.md
```

Reruns are deterministic: seeds are fixed and PDFs carry a fixed creation date, so an unchanged figure stays
byte-identical in git.

## How the analysis runs

```
data/raw  UniProt FASTA + TSV, MetaPredict IDRome, GenOrigin, APPRIS, TFIso1.0 Y2H, GENCODE v50,
  │       UniProt->Ensembl cross-references, RefSeq exon-density files
  │  + data/external/hpc_snapshot   TF labels and gene ages computed on the AKEY cluster
  ▼
notebooks/01_cohort_and_idrs             filters (Swiss-Prot, named gene, ≥20 aa, PE<5, ≥3 isoforms/gene),
  │                                      TF labels, IDR segments ≥20 aa  -> data/annotated/
  ▼
notebooks/02_proteome_analysis           gene age, isoform counts, %IDR range   -> figures/proteome/
  ├─ notebooks/03_appris_annotation      APPRIS gene annotation (read by the pipeline)
  ├─ notebooks/04_tfiso_ppi_counts       TFIso1.0 isoform PPI counts
  ├─ notebooks/05_proteome_analysis_extended                                  -> figures/proteome_extended/
  └─ notebooks/06_exon_density_length_matched  (Susie Song, RefSeq inputs)    -> figures/exon_density/
  ▼
pipeline/run_pipeline.py                 UniProt->Ensembl candidates; GENCODE exon/CDS coordinates; exact
  │                                      translation match; IDR->exon overlap; genomic projection with a DNA
  │                                      round-trip check; master and gene tables; QC; splicing x IDR analyses;
  │                                      transcript resolution; splice events; exon table  -> data/genomic/,
  │                                      figures/genomic/
  ▼
scripts/figures/*.py                     figures drawn from the tables above   -> figures/*, results/tables/
scripts/build_figure_catalog.py          previews + FIGURES.md
```

## Definitions

| Term | Definition |
|---|---|
| Cohort | 21,801 reviewed human UniProt isoforms from 5,143 genes with ≥3 isoforms (5,143 canonical, 16,658 alternative) |
| TF | Gene labelled "Is TF? = Yes" in HumanTFs v1.01 (Lambert et al. 2018): 416 genes, 1,888 isoforms |
| IDR | MetaPredict disordered segment of **≥20 residues** (36,862 segments in the cohort) |
| %IDR (`pct_idr`) | 100 × residues in IDRs / protein length |
| %IDR range (`pct_change_idr`) | Maximum minus minimum %IDR across a gene's isoforms |
| Canonical / alternative isoform | UniProt's canonical sequence / every other isoform of the gene |
| Gene age | GenOrigin age in Ma, binned <100, 100–500, 500–1000, >1000 |
| MAIN set | 15,888 isoforms with one resolved GENCODE transcript (unique or differing only in UTRs) |
| IDR exon | Coding exon with ≥50% of its residues in IDRs (`idr_frac ≥ 0.5`) |
| **Alternative exon** | Exon present in **some but not all** of its gene's mapped transcripts (genes with ≥2 mapped transcripts); **constitutive** = present in all of them. Column `exon_usage` of `data/genomic/processed/exon_level_table` (`usage` in `data/genomic/analysis/exon_idr_annotation.csv`). The table also keeps `is_constitutive`, computed over *all* GENCODE coding transcripts, for reference only; no figure uses it |
| Removed / added IDR exon | IDR exon present only in the canonical / only in the alternative isoform of a pair (relative to the canonical, not evolutionary loss) |

## Repository layout

```
README.md  FIGURES.md  requirements.txt  run.py
notebooks/                 proteome-level analysis (run in number order by run.py)
  01_cohort_and_idrs  02_proteome_analysis  03_appris_annotation  04_tfiso_ppi_counts
  05_proteome_analysis_extended  06_exon_density_length_matched
  go_enrichment_CLUSTER_ONLY   needs GO annotation files on the AKEY cluster; not run by run.py
pipeline/                  genomic-coordinate pipeline (package tfidr_pipeline); README.md, METHODOLOGY.md
scripts/
  figures/                 standalone figure scripts
  build_figure_catalog.py  FIGURES.md + previews from figures/catalog.tsv
data/
  README.md, manifest.tsv  inputs: where to get them, SHA-256
  raw/                     third-party inputs, read-only                          [not in git]
  external/hpc_snapshot/   TF labels and gene ages from the cluster (see its README)
  annotated/               tables written by the notebooks                        [not in git]
  genomic/                 tables written by the pipeline                         [not in git]
figures/
  catalog.tsv              one row per figure (edit the `paper` column)
  proteome/  proteome_extended/  genomic/  exon_density/  go_enrichment/  previews/
results/
  tables/                  tables behind standalone figures
  run_manifest.json        environment, parameters, timings and output checksums of the last full run
docs/PROVENANCE.md         where the code came from, reconstructed figures, decisions, open issues
```

## What is not rebuilt here

* **MetaPredict IDR prediction** was run once in a Colab notebook (Holehouse lab: https://colab.research.google.com/github/holehouse-lab/ALBATROSS-colab/blob/main/idrome_constructor/idrome_constructor.ipynb). Its
  output, `IDRome_AllSwissProtHumansProIsos.csv`, was inputted.
* **TF labels and gene ages** need an Ensembl-gene→UniProt table that exists only on the AKEY cluster. Their
  per-gene results are frozen in `data/external/hpc_snapshot/`; the notebooks use the cluster file when it is present.
* **GO enrichment** (`notebooks/go_enrichment_CLUSTER_ONLY.ipynb`, `figures/go_enrichment/`) needs GO files on the cluster.
* **Susie Song's RefSeq IDR/exon BED files** (`data/raw/exon_density/`) come from her cluster pipeline.

## Known issues

Open issues are currently outlines in [docs/PROVENANCE.md](docs/PROVENANCE.md#open-issues) and in each figure's notes in [FIGURES.md](FIGURES.md).

## License and citation

To be decided before public release.
