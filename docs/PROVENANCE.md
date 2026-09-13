# Provenance

Where the code and figures in this repository came from, how far each can be trusted, the analysis decisions
behind them, and the issues still open.

## Origin

This repository was exported on 2026-09-13 from the project's working ("lab") repository, `tphan/`, at commit
`c6e44d5`, by `reproducibility/export_release.py` in that repository. The export copies the analysis code and
reorganises paths; it does not change analyses. Each path rewrite is asserted in the export script. The lab
repository stays as the record of the audit and cleanup: the pre-cleanup snapshot, the comparison harness, run
records, the recovered-script transcripts, archived originals of superseded or hand-edited figures, the
poster, and side projects that are not part of this paper.

| Here | Lab repository |
|---|---|
| `notebooks/01_cohort_and_idrs.ipynb` | `notebooks/00_data_processing.ipynb` |
| `notebooks/02_proteome_analysis.ipynb` | `notebooks/01_analysis.ipynb` |
| `notebooks/03_appris_annotation.ipynb` | `notebooks/appris_annot.ipynb` |
| `notebooks/04_tfiso_ppi_counts.ipynb` | `notebooks/03_ppi_counts.ipynb` |
| `notebooks/05_proteome_analysis_extended.ipynb` | `notebooks/analysis_2.ipynb` |
| `notebooks/06_exon_density_length_matched.ipynb` | `susie/2024_06_10_exon_enrichement.ipynb` (Susie Song) |
| `notebooks/go_enrichment_CLUSTER_ONLY.ipynb` | `notebooks/02_go_enrichment.ipynb` |
| `pipeline/tfidr_pipeline/`, `pipeline/run_pipeline.py` | `genomiccoords/src/tfidr_pipeline/`, `genomiccoords/src/run_all.py` |
| `scripts/figures/` | `scripts/figures/` |
| `run.py` | `scripts/reproduce_all.py` |
| `figures/proteome/` | `figures/main/` (+ `figures/miscellaneous/tf_canonical_vs_isoform_length_delta.pdf`) |
| `figures/proteome_extended/` | `figures/analysis_2/` (earlier `figures/claude_figs/`) |
| `figures/genomic/` | `genomiccoords/figures/`, its `preliminary_analysis/` and `isoform_transcript_mapping/` |
| `figures/exon_density/` | `susie/figures/` |
| `figures/go_enrichment/` | `figures/miscellaneous/Q4_*` |
| `data/genomic/` | `genomiccoords/data/` |
| `data/raw/gencode/`, `data/raw/uniprot/uniprot_ensembl_xref_raw.tsv` | `genomiccoords/data/cache/` |
| `data/raw/exon_density/` | `susie/` (FASTA, IDRome, xlsx) and `susie/metapredict_results/` |

Not exported: the GENCODE-REST pilot module `mapping/ensembl_exons_pilot.py` (a one-off, not part of the run),
the pipeline's narrative notebooks, and the notebooks' saved outputs (`run.py` writes executed copies to
`results/executed_notebooks/`).

## Verification of this repository

On 2026-09-13 this repository was rebuilt from its raw inputs alone. The run used a new virtual environment
installed from `requirements.txt` (Python 3.13.14, Windows 11), passed the checksums of all 21 inputs, and
`run.py` finished every step with exit code 0 in 23 minutes (notebooks 4 min, pipeline 15 min, figure scripts
3 min). Each of the 177 outputs was then compared with the same file in the lab repository
(`reproducibility/compare_release.py` there):

| Result | Outputs |
|---|---:|
| byte-identical | 71 |
| tables with identical values (only the run-date stamp differs) | 34 |
| PDFs with identical rendering and text (bytes differ only in the embedded creation date) | 71 |
| `data/genomic/qc/qc_report.md`: differs only in its generation date | 1 |

No output differs in content. Lab files with no counterpart here are the GENCODE and UniProt downloads (kept
under `data/raw/`), older unused GENCODE v46 and basic-annotation downloads, two hand-written README notes, the
pilot module's outputs, and legacy tables that nothing in the workflow reads (`gene_df.csv` and the pre-filter
`*UNFILTERED*` tables).

The committed repository was then cloned fresh, given only `data/raw/`, and rebuilt with `python run.py` (all five
steps, exit code 0, 23 minutes). Every tracked file (figures, previews, `FIGURES.md`, `results/tables/`) came out
byte-identical, and all 175 outputs of the first run, data tables included, had the same checksums. The only
changed file was `results/run_manifest.json`, which records that run and is the version committed here.

## How much to trust each figure

The 2026-09 audit reran the whole project from raw inputs and compared every table and figure with the
pre-cleanup files (tables by value, PDFs by rendered pixels and by every printed number). Figures fall into
four groups; each figure's notes in [FIGURES.md](../FIGURES.md) say which group it is in.

1. **Original code** (notebooks and pipeline stages): reruns reproduce the saved outputs.
2. **Recovered scripts** — 11 figures from 10 scripts that existed only in a Claude Code session transcript: B1,
   B2 (fig2 canonical vs alternative), D8 (family disorder), E1 (fig4A stars), E14 (fig09), E15 (fig15b), E16
   (fig16b), F1 (IDR position), F3 (fig03), G12 (fig15 TF vs non-TF), G13 (sticker grammar). The scripts were
   replayed byte-for-byte from the transcript. On the July inputs their outputs were pixel-identical to the
   archived figures, except F1 (0.3% of pixels, because it now reads the current IDR table). The ≥20 aa and
   alternative-exon decisions below have since changed some of their numbers.
3. **Reconstructions** — 16 figures whose code was lost: A1, A3, B5, B8, C11, C13, D5, D7, E4, F2, F4, F7, G2, G5,
   G6, G7. For each, candidate definitions were tested until every number printed on the archived original was
   reproduced from the July tables. Layout differs, and a definition that matches every printed number could
   still differ from the original in something the figure does not print.
4. **Cluster only** — the 4 GO figures (H1–H4) were copied, not regenerated.

In the lab repository, numbers differed for one set of figures: the June 2026 renders of fig1A, fig1B, fig2A/2B
(all isoforms), fig4A (p-value box) and fig5A–D. They were drawn on the cluster from a slightly different table
(1,885 instead of 1,888 TF isoforms). Every rebuild gives 1,888, and no significance call changes. The figures
here are the rebuilt versions.

## Analysis decisions

| Decision | Detail |
|---|---|
| IDR = MetaPredict segment ≥20 aa, everywhere | The proteome notebooks always used ≥20 aa. The genomic pipeline used every segment until 2026-09-12 and now reads `MIN_IDR_LEN_AA` (default 20; `TFIDR_MIN_IDR_LEN_AA=0` reproduces the legacy results). This moved several genomic numbers (e.g. fig14 all genes 39.7% vs 25.7% → 38.3% vs 25.5%; fig03 caption ~4× → ~3×; TF frame-symmetry p 0.035 → 0.068). **Poster numbers for genomic panels predate this.** |
| Alternative exon = present in some but not all of the gene's mapped transcripts | Decided 2026-09-13. The GENCODE-wide definition (absent from any GENCODE coding transcript) labels ~92% of coding exons as alternative and would erase the fig13 result (p=0.86 instead of 1.4×10⁻¹⁰). Five figures that had used it were switched: fig06 all isoforms 31% vs 23% → 38% vs 26%, TF 56% vs 40% → 61% vs 50%; fig07; fig09 (the 100–500 Ma bin goes from **** to ns, and the youngest bin becomes * with non-TFs higher); fig10 (35 → 25 qualifying families, different top 12); sticker-grammar baseline n |
| `IDRome_AllSwissProtHumansProIsos.csv` is ground truth | The author's MetaPredict Colab run; not regenerated |
| TF labels and gene ages from the cluster snapshot | The cluster's Ensembl-gene→UniProt table gives 416 TF genes; local mappings give 414 and 45 gene-age differences (see `data/external/hpc_snapshot/README.md`) |
| fig5 families: deterministic tie-break | The 7th family ties on gene count; broken by gene count, isoform count, then name (Rel over HMG/Sox), which matches the archived figure |
| pandas 3 | The notebooks were written for pandas 2; string-dtype assignment, `groupby.apply` key columns and `sum(numeric_only=…)` were updated, and outputs were verified identical |

## Open issues

Scientific — flagged for the paper, not changed:

1. **Multiple testing is uncorrected** in B1/B2 (4 tests), B7 (6), D1–D5 (~10 per figure), E1/E2 (4 bins), E14–E16 (per bin), F5–F7 (per bin), G13 (4).
2. **F1 "Permutation P < 10⁻³⁰"** is a normal-approximation p from a z-score over 200 circular shifts; an empirical permutation p cannot go below ~1/201. Report it as z-based or as p ≤ 0.005.
3. **F3 (fig03) is not length-matched.** IDR segments (median 69 aa) are much shorter than ordered segments (median 191 aa), and shorter segments span fewer exons by construction.
4. **G9 panel C / s9 gene model**: the unadjusted TF coefficient on %IDR range is not significant (β=0.97, t=0.80). It moves to −0.21 with controls and −0.76 after adding IDR-splicing (β=2.74, t=7.2). The title's "TF effect vanishes … IDR-splicing drives it" reads as mediation; the model shows association only.
5. **G5/G6 (fig13) unit**: isoform × exon rows, so a constitutive exon is counted once per isoform. With unique exons the TF difference holds (0.69 vs 0.56, p=3.7×10⁻¹¹).
6. **"Removed"/"added" IDR exons (E15, G11, G12)** are relative to the UniProt canonical isoform, not evolutionary loss or gain.
7. **Selective display**: D7 (fig10) shows the top 12 of 25 qualifying families; D5 picks 3 of 5 families tied at median 5 isoforms/gene by sort order; A1's bars are not nested.
8. **Duplicates** (pick one of each pair): B1/B3, B2/B4, C10/C11, E4/E5, G1a/G2, G11/G12, G5/G6; C5 overlaps C4, C7 overlaps C6, E6 overlaps E3.
9. **GO (H2–H4)**: no plotted term passes FDR 0.05.
10. **Cohort selection**: genes with ≥3 isoforms are a selected subset of the proteome; APPRIS flags from `03_appris_annotation` are gene-level.
11. **Poster numbers** for genomic panels use the legacy IDR definition and the GENCODE-wide exon definition where applicable.

Presentational:

12. **D6 (figG)**: the "n genes" side bars are drawn in reverse row order (the heatmap is correct).
13. **A5 (fig7A)**: the legend overlaps the non-TF p-value box.
14. **E4 (fig4C)**: the file name says "pctchangeidr" but the figure shows %IDR.

Computational:

15. The cluster-only `ensg_to_upkb.parquet` (and `DatabaseExtract_v_1.01.csv`) would make TF labels and gene ages rebuildable if copied off the cluster.
16. `uniprot_ensembl_xref_raw.tsv` came from the live UniProt API; a new download reflects the current release.
17. The pipeline stamps `data_release_uniprot = "2025_reviewed_human"`, a label rather than a release identifier.
