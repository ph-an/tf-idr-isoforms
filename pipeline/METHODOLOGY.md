# Methodology

Goal: place transcription-factor (TF) isoform IDRs onto exact GRCh38 coordinates and test whether alternative splicing preferentially rewires IDRs. All coordinates are GENCODE v50 / GRCh38.

**Design decision (reuse, not rebuild).** Validated annotation already existed and was reused rather than regenerated: TF family, gene age, APPRIS, and summarized IDR features come from `data/annotated/`; the sequence-validated UniProt→Ensembl mapping backbone is reused. This pipeline fills the genuine gaps — exon/CDS coordinates at scale, IDR↔exon overlap, genomic projection, and the assembled master/gene/analysis tables. Release, paths, and thresholds are centralized in `config.py` (never hard-coded) and stamped onto every output.

**Primary sources.** UniProt reviewed human isoforms (`IDRisoforms_df_geneage.csv`); gene-level APPRIS (`IDRisoforms_df_appris.csv`); per-segment MetaPredict IDR intervals (`IDRome_AllSwissProtHumansProIsos.csv`); GENCODE v50 comprehensive GTF + `pc_translations` + `pc_transcripts`.

**Reading the percentages.** Unless stated otherwise, all percentages are of the 21,801 starting isoforms. Nothing is physically deleted — a "removed" isoform is one that fails a filter's eligibility flag and is held out of that analysis tier (see the attrition funnel before _Which population to use_). All counts below are from the pipeline run with the comprehensive GTF (GENCODE v50 / GRCh38), IDR segments ≥20 aa, and the mapped-transcript definition of alternative exons (2026-09).

## Data flow & reproducibility

The numbered steps below are a **narrative** order. The code executes as stages `s1`–`s9`, four post-modules (`resolve`, `events`, `exons`, `igv`) and the analyses `s9b`, `s9c`, `s9d`, `s10b`, driven by `pipeline/run_pipeline.py`. Every module falls into one of four roles:

| Role                           | Modules (methodology step)                                                                                                    | What it does                                                                                                             |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| **Data acquisition** (network) | `fetch_uniprot_ensembl` (1) · `s1_download_gencode` (feeds 3/5)                                                               | Streams UniProt xrefs; downloads the GENCODE GTF + FASTAs. _Plus the reused `data/annotated/` + MetaPredict inputs (0)._ |
| **Coordinate extraction**      | `s2_exon_cds_coords` (3) · `s4b_genomic_projection` (5)                                                                       | Parses the GTF into exon/CDS tables; builds the genomic position arrays and writes the `.bed` coordinates.               |
| **Validation / filtering**     | `s3_validate_translations` (2) · `resolve_isoform_transcripts` (8) · `s8_qc_and_subsets`                                      | Assigns match types, resolves transcripts by CDS structure, writes QC + subset tables.                                   |
| **Analysis / assembly**        | `s4` (4) · `s5` (6) · `s6` (7) · `exon_level_table` (9) · `classify_splice_events` (10) · `s9` (11) · `s7`, `make_igv_tracks` | Joins, group-bys, statistics, and exports — pure computation on already-loaded data.                                     |

**Reproducibility crux:** only **two modules touch the network** (`fetch_uniprot_ensembl`, `s1_download_gencode`); pin those inputs (UniProt release + GENCODE v50) and every downstream stage is deterministic. The basic→comprehensive GTF switch reaches only the acquisition/extraction layer (`s1`, `s2`; `s4b` uses the unchanged `pc_transcripts`) — everything from Step 6 onward simply recomputes on the new coordinates with no logic change.

---

## 0 — Starting dataset

21,801 reviewed human UniProt isoforms (100%) from 5,143 genes (`IDRisoforms_df_geneage.csv`): 5,143 canonical (23.6%) + 16,658 alternative (76.4%); 416 TF genes / 1,888 TF isoforms (8.7%). Carries protein sequence, length, TF family, gene age, and summarized %IDR; MetaPredict supplies the individual IDR amino-acid intervals and biophysics.
_Filter:_ reviewed human proteins from genes with ≥3 isoforms. No isoforms are ever deleted — later stages add per-analysis eligibility flags instead.

## 1 — Candidate Ensembl transcripts

Parsed UniProt's isoform-tagged `xref_ensembl_full` field (e.g. `ENST…; ENSP…; ENSG…. [P04637-1]`) to link each isoform to candidate ENST/ENSP/ENSG IDs. _Assumption:_ project-canonical accessions (bare, `P04637`) are reconciled to UniProt's `-1` tag as the same isoform; alternatives must match their exact `-N` tag. Version suffixes kept in raw columns, stripped in normalized columns for cross-release joins.
_Coverage:_ 1 candidate ENST 11,688 (53.6%); several 4,561 (20.9%); gene-only ENSG 5,461 (25.0%); no mapping 91 (0.4%) → **≥1 ENST for 16,249 (74.5%)**; the 5,552 with no isoform-specific ENST (25.5%) are removed from transcript/exon/genomic analyses here.
_Out:_ `isoform_transcript_map.csv` (per isoform), `isoform_transcript_candidates_long.csv` (per candidate).

## 2 — Sequence-validate each ENST

Compared each candidate ENST's local GENCODE translation against the UniProt isoform sequence (no REST calls). Classes: **exact** (identical) · near_exact (≥99% identity or clean prefix) · length_match_only · mismatch · no_translation.
_Assumption (conservative):_ **valid = exact only**; all weaker classes are retained for QC but excluded from transcript/exon/genomic analyses. For isoforms with several candidates, all are scored and the best-ranking (exact > near_exact > … ) becomes `resolved_transcript`; the original primary ENST is preserved for reproducibility.
_Result:_ 15,997 / 16,249 primary ENSTs validate exactly (98.4% of those with an ENST; 73.4% of the 21,801 start); scoring the alternate candidates too raises the sequence-valid resolved set to 16,007. The 252 non-exact primaries (1.6% of those with an ENST) are removed from downstream coordinate work.
_Out:_ `mapping_validation_scaled.parquet`.

## 3 — Transcript exon & coding-coordinate maps

From the GENCODE v50 comprehensive GTF, pulled genomic exon and CDS coordinates for mapped ENSTs. Coding exons ordered 5′→3′ with cumulative CDS offsets, linking amino-acid position ↔ CDS nucleotide ↔ genomic coding exon. CDS phase and strand retained so split codons and minus-strand transcripts project correctly.
_Assumption/caveat:_ uses the **comprehensive** annotation so every sequence-valid ENST can receive coordinates (basic omitted ~712 of them; `valid_no_cds` fell 763→51 after the switch). The comprehensive set also adds low-support transcripts (NMD, readthrough, partial), but these only survive if they pass exact sequence match (Step 2) and the genomic round-trip (Step 5), so they are filtered by design.
_Out:_ `ensembl_transcript_exon_coordinates.parquet`, `transcript_cds_to_protein_coordinate_map.parquet`.

## 4 — IDR/exon overlap in protein space (approximate)

Overlaid MetaPredict IDR aa-intervals onto coding-exon aa-ranges; per IDR and non-IDR complement segment computed length, overlapping-exon count, boundary crossing, internal boundary count, and boundary density/kb.
_Caveat:_ aa-space only — approximate because codons can split across exons (±1 boundary ambiguity). Step 5 replaces this with the exact genomic projection.
_Out:_ `idr_exon_overlap_table.parquet`, `nonidr_exon_overlap_table.parquet`.

## 5 — Project IDRs to GRCh38 and self-verify

For each validated isoform's resolved ENST, built a strand-aware flat array of one genomic position per CDS nucleotide (walking exons in rank order). An IDR at `[start,end]` maps to CDS nucleotides `[(start−1)×3 : end×3]`; slicing the array and collapsing runs yields genomic (BED) intervals, so exon-piece counts come from actual genomic breaks (resolving Step 4's split-codon ambiguity).

_Two round-trip checks (the correctness guarantee):_ (a) scaffold — the full CDS must translate back to the intended UniProt protein; (b) segment — every projected slice must translate back to its protein subsequence. Failures logged and excluded.
_Result:_ 15,889 isoforms pass scaffold validation (72.9% of the 21,801 start; 99.3% of the 16,007 sequence-valid resolved transcripts — the other 118, 0.7%, fail the round-trip and are excluded); 51,955 segments projected (26,618 IDR segments ≥20 aa and 25,337 non-IDR segments), 100% round-trip valid.
_Out:_ `genomic_scaffold_validation.parquet`, `genomic_projection_segments.parquet`, `idr/nonidr_genomic_intervals.parquet`, `idrs.bed`/`nonidrs.bed`/`exons.bed`, `genomic_projection_failures.csv`.

## 6 — Isoform master table

Left-joined mapping, validation, projection, TF annotation, gene age, APPRIS, and IDR measures onto all 21,801 isoforms. Nothing deleted; eligibility flags added: `analysis_include_isoform_level` (all), `_exon_level` (sequence-valid), `_genomic_projection` (scaffold-valid), plus `exclude_reason`.
_Caveat:_ APPRIS here is reused gene-level annotation; transcript-level principal flags come directly from GENCODE in the exon table (Step 9).
_Out:_ `tf_idr_isoform_master.{csv,parquet}`.

## 7 — Gene-level summary

One row per base UniProt accession. Main within-gene variable `pct_change_idr_gene = max isoform %IDR − min isoform %IDR`, plus isoform counts, mapping coverage, canonical-vs-alternative IDR differences, lengths, and boundary measures.
_Out:_ `tf_idr_gene_summary.{csv,parquet}` — 5,143 genes.

## 8 — Resolve multiple valid transcripts by CDS structure

Never `.iloc[0]`. For each isoform with ≥1 valid ENST, built a coding-structure signature (chromosome, strand, ordered CDS-exon intervals) and classified: resolved_unique 11,528 (52.9%); resolved_collapsed_equivalent 4,360 (20.0%, ENSTs differing only in UTRs collapsed); ambiguous_multiple_cds 68 (0.3%, flagged, never auto-picked); valid_no_cds 51 (0.2%); unresolved_no_valid 242 (1.1%); gene_only 5,461 (25.0%); no_mapping 91 (0.4%).
_Key result:_ the 20.9% "multiple candidate" ambiguity is almost all UTR-only redundancy — only 68 isoforms (0.3%) are genuinely ambiguous by coding structure. **MAIN set** = resolved_unique ∪ collapsed_equivalent = 15,888 isoforms (72.9%, ~equal TF/non-TF retention: 73.9% TF vs 72.8% non-TF); the remaining 5,913 (27.1%) are held out of chosen-transcript analyses.
_Out:_ `isoform_transcript_mapping_resolved.{csv,parquet}`.

## 9 — Comprehensive exon table

One row per MAIN isoform × exon of its chosen transcript: ENSG/ENST/ENSP/exon IDs; genomic/CDS/aa coordinates; CDS phase; per-exon IDR fraction and region class; canonical/scaffold flags; MANE, Ensembl-canonical, transcript-level APPRIS; alternative vs constitutive status among the gene's mapped transcripts (`exon_usage`, from S9 — the definition used by all analyses); GENCODE-wide inclusion fraction and `is_constitutive` (reference only).
_Result:_ 219,377 rows across 15,888 isoforms / 5,102 genes (99.2% of the 5,143 genes; the 41 missing genes, 0.8%, retain no MAIN isoform). Of 206,445 coding exon rows, 62,560 are IDR-encoding (`idr_frac ≥ 0.5`). By `exon_usage` (Step 11 definition) 76,500 are alternative, 120,595 constitutive, 8,658 belong to genes with one mapped transcript and 692 to isoforms outside the Step 11 table; under the GENCODE-wide reference column `is_constitutive` only 16,658 (8.1%) would be constitutive.
_Out:_ `exon_level_table.{csv,parquet}` (+ data dictionary).

## 10 — Classify canonical-vs-alternative splice events

For genes where canonical and ≥1 alternative isoform are both MAIN-resolved with CDS coordinates, compared coding-exon structures and classified each differential exon: exon skip · alt 5′ ss · alt 3′ ss · alt first/last exon · intron retention · complex. _Assumption:_ an exon with IDR fraction ≥ 0.5 is IDR-encoding, so events are labeled adding/removing/resizing an IDR.
_Result:_ 4,204 genes (81.7% of 5,143 — the rest lack two MAIN-resolved isoforms to compare), 10,151 canonical–alternative pairs, 39,422 events (IDR removal 10,290 outweighs addition 2,492).
_Out:_ `splice_events.{csv,parquet}`, `splice_event_idr_summary.csv`.

## 11 — Does splicing preferentially change IDRs?

Restricted to scaffold-valid isoforms. Exon classes: IDR (frac ≥ 0.5) · ordered (≤ 0.1) · mixed (between) · symmetric (CDS length divisible by 3). Tested: IDR vs non-IDR boundary density; enrichment of frame-symmetric IDR exons; IDR content of alternative vs constitutive exons (an exon is **alternative** if present in some but not all of the gene's mapped transcripts, among genes with ≥2; **constitutive** if present in all); whether IDR-targeted splicing explains gene-level IDR change; biophysics of spliced IDRs; TF vs non-TF exon IDR fraction. Gene-paired tests reduce pseudoreplication; the gene-level regression controls for isoform count, protein length, and gene age.
_Main result:_ alternative exons encode an IDR more often than constitutive exons (38.3% vs 25.5%, OR 1.82, gene-paired p = 7.7e-172; TFs 61.0% vs 49.8%, OR 1.58). In the gene-level model the TF coefficient on `pct_change_idr` is not significant before adjustment (β = 0.97, t = 0.80); it becomes −0.21 with controls (isoform count, length, gene age) and −0.76 after adding IDR-targeted splicing, whose own coefficient is β = 2.74 (t = 7.2). This is an association; the model does not establish mediation.
_Out:_ `exon_idr_annotation.csv`, `splicing_idr_*.csv`.

---

## Cohort attrition (isoforms)

Each tier is a superset filter on the 21,801 start; "removed" = held out of that tier, not deleted.

| Tier                                | Isoforms | % of start | Removed vs start |
| ----------------------------------- | -------: | ---------: | ---------------: |
| Starting dataset                    |   21,801 |       100% |                — |
| ≥1 candidate ENST (Step 1)          |   16,249 |      74.5% |    5,552 (25.5%) |
| Exact sequence-valid (Step 2)       |   15,997 |      73.4% |    5,804 (26.6%) |
| Genomic-scaffold validated (Step 5) |   15,889 |      72.9% |    5,912 (27.1%) |
| MAIN structurally resolved (Step 8) |   15,888 |      72.9% |    5,913 (27.1%) |

MAIN (structure-resolved) and scaffold-validated (round-trip) are now near-identical (15,888 vs 15,889) though defined by _different_ criteria; the most conservative coordinate set is their intersection. The comprehensive GTF lifted both coordinate tiers by ~760 isoforms over the previous basic-GTF run (which sat at ~15,150).

## Which population to use for a claim

| Claim                                       | Population                        |
| ------------------------------------------- | --------------------------------- |
| Protein sequence, %IDR, TF family, gene age | All 21,801 isoforms               |
| Transcript / exon analysis                  | Sequence-valid mapping            |
| Exact genomic-coordinate claim              | `genomic_scaffold_validated=True` |
| Chosen-transcript splice comparison         | MAIN resolved set (15,888)        |
| Most conservative coordinate analysis       | MAIN ∩ scaffold-validated         |

## Key caveats

- **Match policy:** only `exact` sequence matches are eligible; `near_exact` is a QC label, never used in analysis.
- **Comprehensive GTF:** the pipeline uses the comprehensive v50 annotation (set in `config.py`), which recovered the sequence-valid transcripts basic omitted — `valid_no_cds` fell from 763 to 51 (712 recovered), lifting the MAIN/scaffold coordinate set to ~15,888. Extra low-support transcripts it introduces (NMD, readthrough, partial) are removed by the exact-match + round-trip gates.
- **Boundary density** is sensitive to per-segment vs pooled aggregation and its direction is not robust — mechanistic conclusions rest on exon _modularity_, not density.
- **Pseudoreplication:** exons/events are nested within genes; naive row-level tests overstate significance, so gene-paired tests are used.
- **Coordinates:** GTF-derived tables are 1-based inclusive; BED files are 0-based half-open.
- **Deferred (not built here):** GTEx PSI, ClinVar/gnomAD, conservation, splice-site/RBP motifs, exon age, SLiM/phosphosite overlap.
