# Single cell transcriptomic analyses

This folder supplies reusable analysis code, organ/condition settings, recovered gene sets, and source provenance for the single-cell portions of the Supplementary Code. It was assembled from the supplied `HFD_Aging_reanalysis` notebooks and saved outputs. The source correspondence is provided in `provenance/code_map.tsv`.

## Contents and execution

Run the commands below from this `single_cell` directory. Install the study's R environment before executing them; the scripts do not install or upgrade packages automatically. Large Seurat objects and cisTarget databases are external inputs, not included in this code package.

| Analysis | Entry point | Settings and inputs |
|---|---|---|
| SCENIC and RSS | `scripts/run_scenic.R` | `config/scenic_execution.tsv`; 20 organ/comparison runs |
| Aging signatures | `scripts/run_aging_signatures.R` | SCT normalized data and the four gene sets in `data/gene_sets.tsv` |
| IFN score | `scripts/run_ifn_scores.R` | RNA normalized data, normal-diet cells, tissue-specific control matching, IFN65 gene set |
| Regulon–pathway associations | `scripts/run_pathway_correlations.R` | Seurat object with extended regulon AUC and Reactome scores; pathway-name CSV |
| CellChat | `scripts/run_cellchat.R` | `config/cellchat_runs.tsv`; all 4 organs × 4 age/diet conditions |
| Fig. 4j regulator bubble plot | `scripts/plot_fig4j_regulators.R` | Included author-supplied `data/regulator_table.csv` |
| Junb high/low gene DEGs | `scripts/run_junb_gene_degs.R` | Explicit RNA or SCT gene-expression assay argument |
| Gene-name exclusion audit | `R/gene_filter_audit.R` | Gene universe from the exact marker/correlation input object |

Examples (replace paths with actual inputs):

```sh
Rscript scripts/run_scenic.R SCENIC_Aging_Adipose /path/to/HFD_Aging_reanalysis /path/to/new_output/adipose_aging /path/to/cisTarget_databases
Rscript scripts/run_aging_signatures.R /path/to/CapillarySinusoidal_AgingComparison_subset_reclustering_SeuratV5.rds /path/to/new_output/aging_signatures
Rscript scripts/run_ifn_scores.R /path/to/Endothelial_subset_reclustering.rds /path/to/new_output/ifn
Rscript scripts/run_pathway_correlations.R /path/to/Capillary_Adipose_SCENIC_regulonAUC_extended_Reactome.rds /path/to/Reactome_annotations.csv /path/to/new_output/pathways
Rscript scripts/run_cellchat.R Adipose_Normal_Young /path/to/Adipose_Merged_CellChat_prep.rds /path/to/new_output/cellchat_adipose_young
Rscript scripts/plot_fig4j_regulators.R /path/to/new_output/regulator_bubble_Adipose_Junb.pdf
Rscript scripts/run_junb_gene_degs.R /path/to/Adipose_SCENIC_regulonAUC_extended.rds RNA /path/to/new_output/junb Type5 "1 year"
```

SCENIC is computationally expensive. Choose an unused output directory for each run. The manifest records the original core counts, nTrees = 500, and seed = 123. Four-group organ-specific runs explicitly subset capillary ECs, or liver sinusoidal ECs, after reading the organ endothelial object.

## Biological interpretation

The code preserves the authors' distinction between positive Pearson r² ranking for individual regulon–pathway associations and Spearman coefficients for the multi-TF network. These correlations and the cell-level IFN Welch tests are descriptive analyses. Cells are not independent biological replicates, r² is not statistical significance, and none of these functions estimates between-animal reproducibility.

SCENIC uses SCT depth-corrected counts, not Pearson residuals or Harmony-corrected expression. The correlation step uses filtered corrected counts, GENIE3 uses log2(filtered counts + 1), and AUCell uses log2(full corrected counts + 1). The code retains both nonduplicated regulon and extended-only outputs because the source notebooks used both for distinct downstream displays. Existing extended regulon assay names should be checked before joining old and newly generated outputs.

The historical name-pattern filter is preserved. Its terms include ribosomal, mitochondrial, predicted-name and explicitly selected genes, but a name match is not a validated noise or gene-biotype label. The audit exports each excluded symbol and its matching term. It should be applied only where the original marker/correlation code used it.

## Source data recovered during code review

`data/gene_sets.tsv` contains Aging_1 (31 genes), Aging_2 (14), Aging_3 (34), Aging_4 (6), and IFN65 (65). These were extracted from the original notebooks' saved vector outputs, with notebook/cell references attached to every row. They are historical inputs, not newly computed results. The Aging gene sets came from specified ranges of the coexpression heatmap row order; the derivation is retained in `source_excerpts/aging_signature_derivation.R`.

`data/regulator_table.csv` is copied from the author's existing table. The Fig. 4j script recreates its final saved bubble-plot code. It does not rerun Metascape/TRRUST enrichment or establish the original enrichment background.

## Input and source correspondence

`R/` contains generalized functions and `scripts/` contains entry points. `input_contracts.md` specifies the required assays and metadata. `source_excerpts/` retains selected detailed analysis blocks, including the Aging signature derivation, the Reactome network visualization, and the four-group RSS plot. Each excerpt identifies its upstream objects. `provenance/code_map.tsv` and `provenance/source_excerpt_cells.tsv` map the functions and excerpts to one-based source notebook cell positions.

CellChat retains the source-specific minimum cell-group sizes: 10 for adipose and 1 for heart, liver and skeletal muscle. The manifest also records the outgoing/incoming communication pattern counts used for each condition. These parameters remain distinct from sample-level biological replication.

## Environment and software checks

See `environment.md` for the versions reported in the Methods and the cisTarget resources. Entry points write `sessionInfo.txt` when executed. The code includes cell-identifier alignment checks and explicit input-assay selection. `checks/check_alignment.R` verifies positive Pearson ranking, Spearman signs and cell-label alignment using small artificial matrices.
