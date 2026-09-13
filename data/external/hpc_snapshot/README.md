# HPC snapshot: TF labels and gene ages

Two per-gene annotations were built on the AKEY cluster from an Ensembl-gene → UniProt table that
is not in this repository:

    /scratch/gpfs/AKEY/ssong/tf_idr_paper/data/idmapping/ensg_to_upkb.parquet

| File | Built by | Rows | Content |
|---|---|---:|---|
| `tf_labels_from_hpc_mapping.csv` | `notebooks/01_cohort_and_idrs.ipynb`, "Annotate transcription factors" cell | 416 | UniProt base accession and HumanTFs v1.01 DBD family of every TF gene in the ≥3-isoform cohort |
| `gene_age_from_hpc_mapping.csv` | `notebooks/02_proteome_analysis.ipynb`, GenOrigin cell | 5,108 | GenOrigin gene age (Ma; `>4290` → 4300) and age bin (`[0,100)`, `[100,500)`, `[500,1000)`, `[1000,∞)`) per UniProt accession |

Both are the exact per-gene results of those cells, extracted from the tables the cluster run saved
(the lab repository's `reproducibility/freeze_hpc_annotations.py`). The notebooks use the original derivation when the
parquet is reachable and these files otherwise.

**Why not rebuild from local files?** No local mapping reproduces the cluster table:

* the lab repository's `data/raw/uniprot/idmapping/ensg_to_upkb.csv` gives 414 TF genes instead of 416 (it holds versioned
  ENSG IDs and stray `ENSGT` gene-tree IDs, and maps CUX1's gene to Q13948 instead of P39880) and 45
  gene-age differences — which changes TF counts and p-values in most figures.
* The unversioned `Ensembl` rows of the lab repository's `data/raw/uniprot/idmapping/idmapping.csv` recover all 416 TF genes
  only with an arbitrary sort order and still differ for 6 gene ages (CUX1/Q13948, EPM2A, TMPO, GNAS,
  NRXN3, IPCEF1).

To remove this snapshot: copy `ensg_to_upkb.parquet` (and `DatabaseExtract_v_1.01.csv`) from the cluster,
point the two HPC paths in those notebook cells at the copies, and rerun `python run.py`.
