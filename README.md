# Supplementary code for transcriptomic analyses

Generalized analysis code and executable workflows accompanying the Nature Aging
manuscript. Methods follow the author-designated latest manuscript dated 4 November 2026.

| Analysis | Folder | Main entry points |
|---|---|---|
| Single-cell RNA-seq and regulons | `single_cell/` | SCENIC, signature scores, pathway correlations, CellChat, Junb and Figure 4j |
| Visium HD | `visium_hd/` | Spatial neighborhoods, LIANA, preprocessing and annotation |
| Xenium | `xenium/` | Bst2–Pira2 proximity for Heart, Muscle and Adipose |
| Gene annotation comparison | `validation/mapping_rctd/` | Mapping tables, paired RCTD runs and agreement summaries |

`Supplementary_Code_Transcriptomic_Analyses.md` presents the generalized workflow.
Each folder's README provides its commands and input specifications. Repeated analyses
use shared functions with organ/sample settings in configuration tables.

Set input paths to the analysis computer's data locations and use a separate output
directory. For Seurat objects backed by BPCells, keep the associated matrix directories
accessible. Software versions described in Methods are listed in
`docs/software_versions.tsv`; execution scripts record the actual session used.

Small gene-set, regulator-summary and spatial annotation tables are included. Large
Seurat objects, raw expression matrices and cisTarget databases are supplied through
input paths. Source notebook and cell references accompany the analysis modules.
