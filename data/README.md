# Data

`python run.py` first checks every input against [`manifest.tsv`](manifest.tsv) (path, SHA-256, size, source,
which step reads it). Inputs are not in git (about 620 MB). Get them as below, or copy the `data/raw/` folder from
an existing copy of the project, and put them at these paths.

## Inputs (`data/raw/`)

| Path | What | How to obtain |
|---|---|---|
| `uniprot/uniprotkb_allproteins.fasta` | UniProtKB/Swiss-Prot human canonical and isoform sequences (42,562 entries) | uniprot.org, query `reviewed:true AND model_organism:9606`, FASTA including isoforms; this copy was downloaded 2025-06-02 |
| `uniprot/uniprotkb_annotated` | The same proteome as TSV with all return fields | uniprot.org, same query, TSV, all columns |
| `uniprot/idmapping/upkb_to_ensemble_humanisos.tsv` | UniProt isoform → Ensembl gene | UniProt ID mapping (UniProtKB AC/ID → Ensembl) |
| `uniprot/uniprot_ensembl_xref_raw.tsv` | UniProt → Ensembl transcript/protein cross-references | Downloaded by the pipeline's `fetch` step from the live UniProt REST API if absent. A new download reflects the current UniProt release and will not match the checksum, so keep this copy for exact reproduction |
| `metapredict/IDRome_AllSwissProtHumansProIsos.csv` | MetaPredict v3.1 IDR segments and ALBATROSS properties for the FASTA above (72,333 segments) | Made by Tommy Phan with the Holehouse-lab `idrome_constructor` Colab notebook. **Ground truth; not regenerated** |
| `genorigin/Homo_sapiens.csv` | GenOrigin gene ages | chenzxlab.hzau.edu.cn/GenOrigin → download → Homo sapiens |
| `appris/appris_data.appris.txt` | APPRIS per-transcript annotation | appris.bioinfo.cnio.es → downloads → Homo sapiens (GENCODE) |
| `ppi/Table_S1.tsv`, `ppi/SuppTable_PairwiseY2HResults.txt` | TFIso1.0 isoform clone library and pairwise Y2H results | Supplementary tables of the TFIso1.0 TF-isoform interaction study (as in the lab repository's `data/raw/ppi/`) |
| `gencode/gencode.v50.{annotation.gtf,pc_translations.fa,pc_transcripts.fa}.gz` | GENCODE release 50, GRCh38 | Downloaded by pipeline stage `s1` from `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/` if absent |
| `exon_density/TF_mRNA_hg38_CDS.fa`, `IDRome_all.csv`, `bed_all_tfs/`, `bed_all_transcripts/` | Susie Song's RefSeq-based TF/proteome IDR and exon files | AKEY cluster (`/home/ss4521/akeylab/`, `/scratch/gpfs/ss4521/akeylab/metapredict_results/`); used only by the exon-density figures |
| `exon_density/1-s2.0-S0092867420304815-mmc7.xlsx` | Lambert et al. 2018 (Cell) Table S7, TF IDR classification | Journal supplement |

## Frozen cluster results (`data/external/hpc_snapshot/`, in git)

TF labels (416 genes) and gene ages (5,108 proteins) that the notebooks derive from a cluster-only
Ensembl-gene→UniProt table. See [`external/hpc_snapshot/README.md`](external/hpc_snapshot/README.md) for why a
local mapping cannot replace them and how to remove the snapshot.

## Derived data (rebuilt by `run.py`, not in git)

| Location | Written by | Main tables |
|---|---|---|
| `data/annotated/` | notebooks 01–04 | `IDRisoforms_df.csv` (isoform master: sequence, TF status, %IDR, IDR properties), `idr_df.csv` (one row per IDR segment), `*_geneage.csv` (+ gene age, isoform counts, %IDR range), `*_appris.csv`, TFIso PPI tables |
| `data/genomic/` | `pipeline/run_pipeline.py` | `processed/tf_idr_isoform_master`, `processed/exon_level_table` (+ data dictionary), `analysis/exon_idr_annotation.csv`, `analysis/splicing_idr_*.csv`, `isoform_mapping/splice_events`, `qc/qc_report.md`, BED/IGV tracks |
| `results/tables/` | `scripts/figures/` | tables behind standalone figures |

## Cohort

42,562 Swiss-Prot human entries → reviewed, named gene, length ≥20 aa, protein existence <5 (isoforms inherit
the canonical's PE) → genes with ≥3 isoforms: **21,801 isoforms, 5,143 genes**. TFs: **416 genes, 1,888 isoforms**.
IDRs: MetaPredict segments **≥20 aa** (36,862 in the cohort). Genomic analyses use isoforms whose GENCODE v50
translation matches the UniProt sequence exactly (MAIN set, 15,888 isoforms).
