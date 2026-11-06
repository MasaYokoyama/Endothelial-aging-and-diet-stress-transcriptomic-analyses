# Run from the repository root or any directory:
# Rscript xenium/run_xenium.R /path/to/R_objects /path/to/new_output
# Append --cell-bootstrap to reproduce the legacy cell-resampling bands.
# Inputs are read-only. Existing output files are never overwritten.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L ||
    (length(args) == 3L && args[3L] != "--cell-bootstrap")) {
  stop("Usage: Rscript run_xenium.R INPUT_RDS_DIR NEW_OUTPUT_DIR [--cell-bootstrap]")
}
input_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
include_bands <- length(args) == 3L
script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, mustWork = TRUE))
packages <- c("Seurat", "SeuratObject", "dplyr", "ggplot2", "RANN")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(dplyr); library(ggplot2); library(RANN)
})
source(file.path(script_dir, "proximity_functions.R"))
manifest <- read.delim(file.path(script_dir, "sample_manifest.tsv"), stringsAsFactors = FALSE)
manifest <- manifest[manifest$run_default, , drop = FALSE]
input_paths <- unique(file.path(input_dir, manifest$object_file))
if (!all(file.exists(input_paths))) stop("Missing final annotated RDS input(s): ",
                                        paste(input_paths[!file.exists(input_paths)], collapse = ", "))
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; choose a new directory to protect prior results.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
write.table(manifest, file.path(output_dir, "sample_manifest_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("Execution UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Cell bootstrap requested:", include_bands),
  "Biological replication: one mouse per organ-age condition, per response letter.",
  "Descriptive spatial analysis; cell-resampling bands are not animal-level confidence intervals.",
  "Spatial dependence is not addressed by the legacy independent-cell bootstrap."
), file.path(output_dir, "analysis_scope.txt"))
input_audit <- data.frame(input_file = input_paths,
                         size_bytes = file.info(input_paths)$size,
                         modified = as.character(file.info(input_paths)$mtime),
                         stringsAsFactors = FALSE)
write.table(input_audit, file.path(output_dir, "input_file_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summary_list <- list()
cell_counts <- list()
for (organ in unique(manifest$organ)) {
  settings <- manifest[manifest$organ == organ, , drop = FALSE]
  stopifnot(nrow(settings) == 2L, setequal(settings$age, c("Young", "Aged")))
  settings <- settings[match(c("Young", "Aged"), settings$age), , drop = FALSE]
  message("Processing ", organ)
  obj <- readRDS(file.path(input_dir, settings$object_file[1L]))
  stopifnot(inherits(obj, "Seurat"), "Xenium" %in% Assays(obj))
  stopifnot(all(settings$fov %in% Images(obj)))
  md <- obj[[]]
  if (!settings$annotation_col[1L] %in% names(md)) stop("Missing annotation column: ", settings$annotation_col[1L])
  spatial_df <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
    make_xenium_spatial_table(obj, settings$fov[i], settings$age[i], settings$annotation_col[i], assay = "Xenium")
  }))
  stopifnot(!anyDuplicated(spatial_df$cell))
  definition <- define_bst2_high_ec(spatial_df, organ)
  result <- run_bst2_pira2_proximity(spatial_df, definition, organ,
                                   thresholds_um = seq(0, 200, 5),
                                   n_boot = 1000, seed = 123,
                                   include_cell_bootstrap = include_bands)
  result$analysis_scope <- "Within-specimen description; one mouse per age; not between-animal inference"
  result$cell_bootstrap_included <- include_bands
  saveRDS(result, file.path(output_dir, paste0(organ, "_Bst2_Pira2_proximity_analysis.rds")))
  write.csv(result$macrophage_distance, file.path(output_dir, paste0(organ, "_macrophage_distances.csv")), row.names = FALSE)
  write.csv(result$cumulative_proximity, file.path(output_dir, paste0(organ, "_cumulative_proximity.csv")), row.names = FALSE)
  write.csv(result$bootstrap_difference, file.path(output_dir, paste0(organ, "_proximity_difference.csv")), row.names = FALSE)
  ggsave(file.path(output_dir, paste0(organ, "_proximity_difference.pdf")), result$proximity_plot, width = 9, height = 6)
  summary_list[[organ]] <- definition$summary
  cell_counts[[organ]] <- spatial_df %>% count(Age, annotation, name = "Cell_n") %>% mutate(Organ = organ)
  rm(obj, spatial_df, result, definition)
  gc()
}
write.csv(bind_rows(summary_list), file.path(output_dir, "Bst2_high_definition_summary.csv"), row.names = FALSE)
write.csv(bind_rows(cell_counts), file.path(output_dir, "annotation_counts_by_age.csv"), row.names = FALSE)
writeLines("Execution completed for Heart, Muscle and Adipose.", file.path(output_dir, "COMPLETED.txt"))
