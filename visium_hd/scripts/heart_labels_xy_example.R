# Generalized methods example; requires inputs and upstream objects documented in README.
# Merge RCTD predictions, treat rejected spots as Unknown, and attach to Seurat metadata.
paths <- c(
  Heart_Young_1 = "RCTD_fullref_cutoff30_Heart_Young_1.rds",
  Heart_Young_2 = "RCTD_fullref_cutoff30_Heart_Young_2.rds",
  Heart_Aged_1  = "RCTD_fullref_cutoff30_Heart_Aged_1.rds",
  Heart_Aged_2  = "RCTD_fullref_cutoff30_Heart_Aged_2.rds"
)
 
read_one_rctd <- function(sample, fp) {
  r <- readRDS(fp)
  cellid <- colnames(r@spatialRNA@counts)
  W <- r@results$weights
  results_df <- r@results$results_df
  stopifnot(all(cellid %in% rownames(W)), all(cellid %in% rownames(results_df)))
  W <- W[cellid, , drop = FALSE]
  pred_type <- colnames(W)[max.col(W, ties.method = "first")]
  spot_class <- as.character(results_df[cellid, "spot_class"])
  pred_type2 <- pred_type
  pred_type2[spot_class == "reject"] <- "Unknown"
  data.frame(cellid = cellid, sample = sample, spot_class = spot_class,
             rctd_type = pred_type, rctd_type_unknown = pred_type2,
             stringsAsFactors = FALSE)
}
 
rctd_merged <- do.call(rbind, Map(read_one_rctd, names(paths), paths))
stopifnot(!anyDuplicated(rctd_merged$cellid))
saveRDS(rctd_merged, "RCTD_merged_cutoff30_rejectUnknown.rds")
 
obj <- readRDS("heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming.rds")
spotclass_vec <- setNames(rctd_merged$spot_class, rctd_merged$cellid)
type_vec <- setNames(rctd_merged$rctd_type_unknown, rctd_merged$cellid)
obj$rctd_spot_class <- unname(spotclass_vec[colnames(obj)])
obj$rctd_type <- unname(type_vec[colnames(obj)])
 
# Retain units with an RCTD assignment, including the explicit Unknown class.
# Prefer the marker/cluster annotation; otherwise use the RCTD label.
use <- !is.na(obj$rctd_type)
labels <- as.character(obj$cluster_naming)
fallback <- is.na(labels) | labels == ""
labels[fallback] <- as.character(obj$rctd_type)[fallback]
map_new <- c(B_cell_lineage="B cell lineage", Cardiac_muscle_cell="Ventricular cardiac muscle",
 Cardiac_neuron="Cardiac neuron", Endocardial_cell="Endothelial cell",
 Endothelial_cell="Endothelial cell", Erythroid_cell="Red blood cell", Fibroblast="Fibroblast",
 Lymphatic_endothelial="Lymphatic endothelial cell", Macrophage___Monocyte___DC="Macrophage / Monocyte / DC",
 Mast_cell="Mast cell", Pericyte="VSMC / Pericyte", T_cell_lineage="T cell lineage",
 Unknown="Unknown", VSMC="VSMC / Pericyte")
i <- labels %in% names(map_new)
labels[i] <- unname(map_new[labels[i]])
# Each supplied Valve CSV contains manually selected base barcodes in its first column.
# Only fibroblast-labelled units in these selections are reassigned.
valve_files <- list.files("annotations", pattern = "^Heart_.*_Valve\\.csv$", full.names=TRUE)
stopifnot(length(valve_files) == 4)
valve_ids <- unlist(lapply(valve_files, function(f) {
 s <- sub("_Valve\\.csv$", "", basename(f))
 ids <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)[[1]]
 paste0(s, "_", gsub("^'+|'+$", "", ids))
}))
labels[which(colnames(obj) %in% valve_ids & labels == "Fibroblast")] <- "Cardiac valve"
labels[!use] <- NA_character_
obj$Merged_cellname <- labels
obj_use <- subset(obj, subset = !is.na(Merged_cellname))
obj_use$Merged_cellname <- factor(obj_use$Merged_cellname)
saveRDS(obj_use, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming_RCTDmerged_MergedCellname_ASSIGNEDonly.rds")
 
# Per-sample Visium HD 8-um objects were filtered and XY coordinates were added from tissue_positions.parquet.
add_xy_heart <- function(obj, pos_path) {
  pos <- read_parquet(pos_path)
  coords_raw <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
                           row.names = pos$barcode)
  cells <- colnames(obj)
  cells_base <- sub("^.*_(s_008um_)", "\\1", cells)
  idx <- match(cells_base, rownames(coords_raw))
  stopifnot(all(!is.na(idx)))
  obj$x <- coords_raw$x[idx]
  obj$y <- coords_raw$y[idx]
  obj
}
 
# The additional spatial-count filter was applied to Heart_Aged_1 and Heart_Young_1.
# The saved source does not apply this 50 threshold to all four samples.
obj_hd_aged1_50 <- subset(obj_hd_aged1_lab, subset = nCount_Spatial.008um >= 50)
obj_hd_aged1_50 <- add_xy_heart(obj_hd_aged1_50,
  "Heart_Aged_1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet")
saveRDS(obj_hd_aged1_50, "Heart_Aged_1_visiumHD8um_ASSIGNEDonly_withXY.rds")
