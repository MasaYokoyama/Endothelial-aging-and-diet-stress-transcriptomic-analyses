merge_heart_labels <- function(obj, valve_files) {
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

stopifnot(length(valve_files) == 4)
valve_ids <- unlist(lapply(valve_files, function(f) {
 s <- sub("_Valve\\.csv$", "", basename(f))
 ids <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)[[1]]
 paste0(s, "_", gsub("^'+|'+$", "", ids))
}))
labels[which(colnames(obj) %in% valve_ids & labels == "Fibroblast")] <- "Cardiac valve"
labels[!use] <- NA_character_
obj$Merged_cellname <- labels

obj
}
