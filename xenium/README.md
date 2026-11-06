# Xenium Bst2–Pira2 proximity analysis

The workflow calculates distances from monocyte/macrophage-lineage cells to the
nearest Bst2-high endothelial-containing segment within each organ and age specimen.
It exports cumulative proximity curves and Pira2-positive minus Pira2-negative differences.

## Run

From the repository root:

```sh
Rscript xenium/run_xenium.R /path/to/R_objects /path/to/output
```

To reproduce the original cell-resampling bands:

```sh
Rscript xenium/run_xenium.R /path/to/R_objects /path/to/output_with_bands --cell-bootstrap
```

The runner processes Heart, Muscle and Adipose, records settings/session information,
and writes `COMPLETED.txt` after all three organs finish. Use a new output directory.

## Inputs

Supply the final annotated objects:

- `heart_qc_annotated_refined_updated.rds`
- `muscle_qc_annotated_refined.rds`
- `adipose_qc_annotated_refined_updated.rds`

Each object requires the normalized Xenium assay `data` layer with Bst2/Pira2,
matching assay/metadata cell identifiers, Young/Aged organ FOVs with centroid
coordinates in micrometers, and the organ-specific annotation column listed in
`sample_manifest.tsv`. Required classes are `Endothelial-containing segment` and
`Monocyte/macrophage lineage`.

## Settings

Heart and Adipose use the all-cell pooled Young/Aged Bst2 90th percentile within each
organ, followed by endothelial-segment selection. Muscle uses Bst2 > 0 within
endothelial-containing segments. Pira2-positive myeloid cells have Pira2 > 0.
Distances are calculated within each age FOV and summarized from 0–200 µm in 5-µm steps.

Each organ–age condition contains one mouse. Curves summarize within-specimen spatial
organization. The optional independent-cell bootstrap uses 1,000 resamples and seed
123 to describe within-specimen resampling variation; it is not animal-level inference.

`proximity_functions.R` contains reusable functions. `preprocessing_settings.tsv`,
`annotation_cluster_maps.tsv` and `annotation_rules.md` document the upstream processing.
`provenance.tsv` links the analysis to original Notebook cells. `notebook_observations.tsv`
records the original saved outputs. R packages are Seurat, SeuratObject, dplyr, ggplot2 and RANN.
