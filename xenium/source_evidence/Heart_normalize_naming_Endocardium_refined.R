# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Heart_normalize_naming_Endocardium_refined.ipynb CELL 5 (one-based)
heart_endocardium_dir <- "/work/Xenium/Heart_endocardium"

export_files <- list.files(
  heart_endocardium_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

artery_files <- export_files[
  grepl("artery", basename(export_files), ignore.case = TRUE)
]

vein_files <- export_files[
  grepl("vein", basename(export_files), ignore.case = TRUE)
]


# SOURCE Heart_normalize_naming_Endocardium_refined.ipynb CELL 6 (one-based)
read_selected_cells <- function(file) {

  x <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  ids <- unique(trimws(x[["Cell ID"]]))
  ids <- ids[!is.na(ids) & ids != ""]

  if (grepl("^young_", basename(file), ignore.case = TRUE)) {
    ids <- paste0("Young_Heart_", ids)
  }

  if (grepl("^aged_", basename(file), ignore.case = TRUE)) {
    ids <- paste0("Aged_Heart_", ids)
  }

  ids
}


# SOURCE Heart_normalize_naming_Endocardium_refined.ipynb CELL 7 (one-based)
artery_selected <- unique(
  unlist(
    lapply(artery_files, read_selected_cells),
    use.names = FALSE
  )
)

vein_selected <- unique(
  unlist(
    lapply(vein_files, read_selected_cells),
    use.names = FALSE
  )
)

c(
  artery_selected = length(artery_selected),
  vein_selected = length(vein_selected)
)


# SOURCE Heart_normalize_naming_Endocardium_refined.ipynb CELL 8 (one-based)
old_label <- "Venous-like large-vessel EC-rich"

original_target_cells <- colnames(heart_qc)[
  as.character(heart_qc$heart_annotation_refined) == old_label
]

artery_cells <- intersect(
  artery_selected,
  original_target_cells
)

vein_cells <- intersect(
  vein_selected,
  original_target_cells
)

endocardial_cells <- setdiff(
  original_target_cells,
  union(artery_cells, vein_cells)
)

c(
  original_target_total = length(original_target_cells),
  artery_matched = length(artery_cells),
  vein_matched = length(vein_cells),
  endocardial_remaining = length(endocardial_cells),
  artery_vein_overlap = length(intersect(artery_cells, vein_cells))
)


# SOURCE Heart_normalize_naming_Endocardium_refined.ipynb CELL 9 (one-based)
heart_qc$heart_annotation_refined <- as.character(
  heart_qc$heart_annotation_refined
)

heart_qc$heart_annotation_refined[artery_cells] <-
  "Artery-like large-vessel EC-rich"

heart_qc$heart_annotation_refined[vein_cells] <-
  "Venous-like large-vessel EC-rich"

heart_qc$heart_annotation_refined[endocardial_cells] <-
  "Endocardial cell-rich"

heart_qc$heart_annotation_refined <- factor(
  heart_qc$heart_annotation_refined
)
