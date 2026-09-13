# Genomic-coordinate pipeline (`tfidr_pipeline`)

Matches each UniProt isoform in the cohort to a sequence-identical GENCODE v50 transcript. It then places the
isoform's IDRs on coding exons and GRCh38 coordinates, checks every projection by translating it back, and runs
the splicing × IDR analyses. Method details, with counts at every step, are in [METHODOLOGY.md](METHODOLOGY.md).

UniProt supplies the protein and isoform backbone and GENCODE the transcript and genome structure. TF family,
gene age, APPRIS and summary IDR features are reused from the notebooks' tables in `data/annotated/`, not
rebuilt. Paths, the GENCODE release and thresholds are set in [`tfidr_pipeline/config.py`](tfidr_pipeline/config.py).

## Running

Normally run as part of `python run.py` from the repository root, after the notebooks have written
`data/annotated/IDRisoforms_df_geneage.csv` and `IDRisoforms_df_appris.csv`. On its own:

```bash
python pipeline/run_pipeline.py                         # every step, in order
python pipeline/run_pipeline.py --from s4 --through s9  # a range of steps
python pipeline/run_pipeline.py --from resolve          # post-processing and downstream analyses only
```

| Step | Module | Does |
|---|---|---|
| `fetch` | `mapping/fetch_uniprot_ensembl` | UniProt → Ensembl cross-references (`data/raw/uniprot/`; downloaded if absent) |
| `map` | `mapping/build_mapping_table` | candidate transcripts per isoform |
| `s1` | `stages/s1_download_gencode` | GENCODE v50 GTF and FASTAs (`data/raw/gencode/`; downloaded if absent) |
| `s2` | `stages/s2_exon_cds_coords` | exon and CDS coordinates, CDS ↔ protein map |
| `s3` | `stages/s3_validate_translations` | exact translation match against UniProt |
| `s4` | `stages/s4_idr_exon_overlap` | IDR / non-IDR segments vs exons in protein space |
| `s4b` | `stages/s4b_genomic_projection` | genomic projection with DNA round-trip check; BED files |
| `s5`–`s8` | `stages/s5_build_master` … `s8_qc_and_subsets` | isoform master, gene summary, figure-ready tables, QC report |
| `s9` | `stages/s9_splicing_idr_analysis` | per-exon IDR table with alternative/constitutive usage; splicing × IDR statistics |
| `resolve` | `mapping/resolve_isoform_transcripts` | one transcript per isoform by CDS structure (MAIN set) |
| `events` | `mapping/classify_splice_events` | splice events between canonical and alternative isoforms |
| `exons` | `exports/exon_level_table` | exon-level table + data dictionary |
| `igv` | `exports/make_igv_tracks` | BED12 tracks for IGV |
| `s9b` | `stages/s9b_idr_enrichment` | alternative-exon IDR enrichment → `figures/genomic/fig14` |
| `s9c` | `stages/s9c_exon_idr_fraction` | continuous exon IDR fraction → `fig16` |
| `s9d` | `stages/s9d_gene_idr_fraction` | per-gene paired IDR fraction → `fig17` |
| `s10b` | `stages/s10b_idr_gain_loss` | IDR exons conserved / removed / added → `fig15` |

Single modules can be run for debugging from `pipeline/`, e.g. `python -m tfidr_pipeline.exports.exon_level_table`.

## Settings

* `MIN_IDR_LEN_AA` (default **20**, the notebooks' definition): shorter MetaPredict segments are ignored.
  `TFIDR_MIN_IDR_LEN_AA=0` reproduces the pre-2026-09 genomic results, which used every segment.
* `GENCODE_RELEASE = "50"` (GRCh38). Changing it re-runs everything against another annotation, and the input checksums will no longer match.
* Only **exact** translation matches are used for transcript, exon or genomic analyses.

## Outputs (`data/genomic/`)

| Location | Main files |
|---|---|
| `processed/` | `tf_idr_isoform_master` (one row per isoform: mapping, validation, TF, age, APPRIS, IDR features, analysis flags); `tf_idr_gene_summary`; `exon_level_table` (+ `_DATA_DICTIONARY.md`); `idrs.bed`, `nonidrs.bed`, `exons.bed`; `igv/` |
| `analysis/` | `exon_idr_annotation.csv` (per mapped isoform × coding exon, with `usage`); `splicing_idr_*.csv` (s9 statistics); `idr_gain_loss_*.csv`; figure-ready `viz_*.csv`; `analysis_*.parquet` subsets |
| `isoform_mapping/` | resolved mapping, `splice_events`, exon region classes |
| `interim/` | exon/CDS coordinates, translation validation, overlap tables, projection tables |
| `qc/` | `qc_report.md`, mapping-status and missingness tables, projection failures |

Figures go to `figures/genomic/`.

## Coverage (GENCODE v50 / GRCh38, IDRs ≥20 aa)

* 21,801 isoforms, 5,143 genes (the notebooks' cohort)
* 16,249 isoforms (74.5%) with ≥1 candidate transcript; 15,997 of those match the GENCODE translation exactly
* **15,888 isoforms in the MAIN set** (11,528 unique + 4,360 UTR-collapsed); 15,889 genomic-scaffold validated; 15,822 in both
* 51,955 segments projected to the genome (26,618 IDR, 25,337 non-IDR), all round-trip validated
* exon-level table: 219,377 exon rows (206,445 coding) across 15,888 isoforms, 5,102 Ensembl genes

## Alternative exons

`exon_usage` in the exon table (`usage` in `analysis/exon_idr_annotation.csv`) is the definition every analysis uses.
Among a gene's mapped transcripts (genes with ≥2), an exon in all of them is **constitutive** and one in some but
not all is **alternative**; exons of genes with one mapped transcript are `single_tx_gene`. `is_constitutive` and
`transcript_inclusion_fraction` count *all* GENCODE protein-coding transcripts of the gene and are kept for reference only.

## Caveats

* Restrict genomic and exon-coordinate claims to validated isoforms (`genomic_scaffold_validated`, the MAIN set).
  Protein-level measures (%IDR etc.) are fine for all isoforms.
* APPRIS fields in the isoform master come from the notebooks' gene-level (Ensembl gene) annotation, so they are
  not per-isoform calls. The exon table's `is_appris_principal` / `is_mane_select` are the chosen transcript's own GENCODE tags.
* Exons and events are nested within genes; the s9 tests are gene-paired where possible.
* GTEx expression is not included (`GTEX_ENABLED = False`).
