# Common exporters, source-matched to original Heart/Adipose CoreNeighborDistal notebooks.
# Heart LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_heart <- function(obj, outdir) {
  DefaultAssay(obj) <- "Spatial.008um"
  obj <- JoinLayers(obj, assay = "Spatial.008um", layers = "data")
  keep_ct <- c("Ventricular cardiac muscle", "B cell lineage", "Cardiac neuron", "Endocardium",
               "Endothelial cell", "Epicardium", "Fibroblast", "Lymphatic endothelial cell",
               "Macrophage / Monocyte / DC", "Mast cell", "Platelet", "T cell lineage",
               "VSMC / Pericyte")
  obj_use <- subset(obj, subset = Merged_cellname %in% keep_ct)
  obj_use$Aging1_role <- dplyr::case_when(
    obj_use$Aging1_nbhd100um == "seed" ~ "Core",
    obj_use$Aging1_nbhd100um == "neighbor" ~ "Neighbor",
    obj_use$Aging1_nbhd100um == "other" ~ "Distal",
    TRUE ~ NA_character_
  )
  obj_use$CCI_group <- dplyr::case_when(
    obj_use$Aging1_role == "Core" ~ "Core",
    obj_use$Aging1_role == "Neighbor" ~ paste0("Neighbor_", obj_use$Merged_cellname),
    obj_use$Aging1_role == "Distal" ~ paste0("Distal_", obj_use$Merged_cellname),
    TRUE ~ NA_character_
  )
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = "Spatial.008um", layer = "data")
  stopifnot(identical(colnames(X), colnames(obj_use)))
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- data.frame(cell_id = colnames(obj_use), CCI_group = obj_use$CCI_group,
                    stringsAsFactors = FALSE)
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
 
# Adipose LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_adipose <- function(obj, outdir) {
  assay_use <- "SCT"
  DefaultAssay(obj) <- assay_use
  obj <- tryCatch(JoinLayers(obj, assay = assay_use, layers = "data"), error = function(e) obj)
  exclude_ct <- c("Erythroid cell", "Unknown")
  obj$CCI_group <- NA_character_
  obj$CCI_group[obj$Aging1_nbhd100um == "seed"] <- "Core"
  obj$CCI_group[obj$Aging1_nbhd100um == "neighbor"] <- paste0("Neighbor_", obj$Merged_cellname[obj$Aging1_nbhd100um == "neighbor"])
  obj$CCI_group[obj$Aging1_nbhd100um == "other"] <- paste0("Distal_", obj$Merged_cellname[obj$Aging1_nbhd100um == "other"])
  obj_use <- subset(obj, subset = !is.na(CCI_group) &
                    (CCI_group == "Core" | !(Merged_cellname %in% exclude_ct)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = assay_use, layer = "data")
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- obj_use@meta.data[, c("CCI_group"), drop = FALSE]
  obs$cell_id <- rownames(obs)
  obs <- obs[, c("cell_id", "CCI_group")]
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
