# Visium HD analysis code

Heart 8-µm-bin and adipose segmented-cell workflows for eight samples.

## Run spatial neighborhoods and LIANA

From this folder, set the input data directory and separate output directories:

```sh
Rscript scripts/run_neighborhood_export.R /path/to/Visium_data /path/to/exports config/samples.tsv
python scripts/run_liana.py /path/to/exports /path/to/liana_results config/samples.tsv
```

The R script reads saved `Aging_1_Score1` scores, constructs core, neighbor and distal
groups, and exports expression matrices with aligned cell metadata. The Python script
runs LIANA and produces per-sample estimates and two-sample age-group summaries.

## Inputs and settings

- `config/samples.tsv`: all eight samples, object paths, assays and spatial thresholds.
- `config/input_objects.tsv`: required object classes and metadata.
- `annotations/`: four left-ventricular ROI and four valve-selection tables.
- `scripts/*_example.R`: preprocessing, reference, mapping and annotation steps.
- `curated_notebooks/`: selected original steps in analysis order.
- `provenance/selected_cells.tsv`: original notebook and cell references.

The core cutoffs are 0.3 for heart and 0.05 for adipose. The neighborhood radius is
100 µm. Heart module scores are calculated before selecting LV regions. Adipose
coordinate conversion follows the original centroid-span/6,500-µm calculation.
The additional count filter of 50 applies to Heart_Aged_1 and Heart_Young_1.

Set the raw Space Ranger paths and BPCells matrix directory in preprocessing examples.
Supply upstream objects as listed in the input table and run steps in their numbered
order. Keep on-disk matrix directories accessible to the RDS objects that use them.

R packages include Seurat, Matrix, RANN and dplyr; Python packages include LIANA,
Scanpy, AnnData, pandas, NumPy and SciPy. Each runner records its package versions.
Gene-ID mapping comparison and paired RCTD execution are in `../validation/mapping_rctd/`.
