# Gene annotation and paired RCTD comparison

## Inputs

The original Visium folder contains `reference_rctd_heart.rds`, the Ensembl-based `heart_mergeobj_preProject_full_sketchHarmony.rds`, a biomaRt-renamed counterpart, all four `RCTD_fullref_cutoff30_<sample>.rds` results, and coordinate files. BPCells matrices may reference directories on the original computer. An RDS filename alone does not guarantee portable access to its on-disk counts.

## Execution

Copy `config.example.R` to `config.local.R` and set the original computer's paths. From this directory run in order:

```sh
Rscript 01_preflight.R config.local.R
Rscript 02_compare_mappings.R config.local.R
Rscript 03_rerun_paired_rctd.R config.local.R
Rscript 04_summarize_results.R config.local.R
```

Alternatively, run all steps with `bash run_comparison.sh config.local.R`.

Step 4 combines per-sample agreement, retention, and primary-label composition tables.

Step 1 reads the runtime and original RCTD configurations without inference. Step 2 compares the same Ensembl universe, records one-to-many mappings, and freezes both tables. Step 3 performs full paired doublet-mode RCTD runs for four samples and measures label agreement, changed classifications, and filter-related loss of bins. It also compares the new org.Mm.eg.db arm with the saved original output to expose environment or stochastic rerun differences.

Prefer the original BioMart query result if recoverable. The original `useMart()` call did not pin an Ensembl release. A newly retrieved table is explicitly labeled as a new retrieval, and must not be described as the exact historical mapping. Display names modified by `make.unique()` are not interchangeable with gene symbols and are not directly used as RCTD input here.

Both arms use the original rule of retaining the first nonempty mapping per Ensembl ID and summing counts for multiple Ensembl IDs assigned to one symbol. All ambiguous mappings are separately listed. The initial filter is raw `nCount_RNA >= 30`; each sample's saved `UMI_min` is used after gene mapping, with retained-bin differences reported.

The output tables report mapping categories, retained spatial units, paired RCTD label agreement, and affected cell types. Use these tables to summarize the annotation comparison and its impact on deconvolution.

References: [spacexr original workflow](https://github.com/dmcable/spacexr/blob/master/vignettes/spatial-transcriptomics.Rmd), [biomaRt archived datasets](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html).
