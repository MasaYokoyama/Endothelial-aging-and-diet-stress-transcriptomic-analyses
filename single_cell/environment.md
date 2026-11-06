# Analysis environment

The study Methods report the following package versions:

| Package | Version |
|---|---|
| Seurat | 5.0.3 |
| SCENIC | 1.3.1 |
| RcisTarget | 1.29.0 |
| UCell | 2.3.0 |
| ggplot2 | 3.5.1 |
| igraph | 1.5.0 |
| CellChat | 1.6.1 |

Additional libraries are loaded explicitly in the corresponding scripts. Each analysis entry point writes the active R `sessionInfo()` to its output directory.

SCENIC uses the following cisTarget resources:

- `mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather`
- `mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather`
- RcisTarget annotation object `motifAnnotations_mgi_v9`

Pass the directory containing these database files to the SCENIC runner. The source workflow uses SCT depth-corrected counts, nTrees = 500 and seed = 123. Each organ/comparison is run in a separate output directory.
