# Generalized methods example. See README and config/input_objects.tsv for input specifications.
# Label transfer from adipose scRNA-seq reference
adipose <- readRDS("adipose_mergeobj_preProject_sketchHarmony.fixed_onDiskPath_joinCounts.rds")
ref_obj <- readRDS("adipose_merged_obj_reference.rds")
DefaultAssay(ref_obj) <- "RNA"
DefaultAssay(adipose) <- "RNA"
 
ref_obj <- FindVariableFeatures(ref_obj, verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
features_use <- intersect(VariableFeatures(ref_obj), VariableFeatures(adipose))
features_use <- head(features_use, 500)
 
anchors <- FindTransferAnchors(reference = ref_obj, query = adipose,
                               features = features_use, dims = 1:30)
pred <- TransferData(anchorset = anchors, refdata = ref_obj$Cellname, dims = 1:30)
adipose <- AddMetaData(adipose, metadata = pred)
adipose$Cellname_final <- ifelse(adipose$prediction.score.max < 0.4,
                                 "Unknown", adipose$predicted.id)
saveRDS(adipose, "adipose_mergeobj_preProject_sketchHarmony_labelTransferred.rds")
 
# Generate the SCT assay used for marker review and downstream module scoring.
# Source: Adipose_Combination_5, cell 17 (one-based).
DefaultAssay(adipose) <- "RNA"
adipose <- SCTransform(adipose, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
# Ambiguous labels were reclassified using Harmony SNN reclustering.
DefaultAssay(adipose) <- "SCT"
adipose <- FindNeighbors(adipose, reduction = "harmony", dims = 1:30,
                         k.param = 50, graph.name = "harmony_snn_k50")
adipose <- FindClusters(adipose, graph.name = "harmony_snn_k50", resolution = 0.2,
                        random.seed = 1)
adipose$clusters_k50_r020 <- adipose[["harmony_snn_k50_res.0.2"]][, 1]
 
# Assign broad cluster identities used to resolve mixed/unknown labels.
cl <- as.character(adipose$clusters_k50_r020)
new_name <- rep("Unassigned", length(cl))
new_name[cl == "0"] <- "Adipocyte (mature; lipid-signal mixed)"
new_name[cl == "1"] <- "Endothelial (capillary~general)"
new_name[cl == "2"] <- "Adipose stromal / Fibroblast"
new_name[cl == "3"] <- "Myeloid (Macrophage/DC; MHCII+; incl M2-like)"
new_name[cl == "4"] <- "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)"
new_name[cl == "5"] <- "B cell"
new_name[cl == "8"] <- "Erythroid"
new_name[cl %in% as.character(10:29)] <- "Adipocyte-like / lipid-signal dominated"
adipose$cluster_based_naming <- new_name
 
# Cluster-based curation of mixed/unknown categories.
adipose$Cellname_partial <- as.character(adipose$Cellname_final)
adipose$Cellname_partial[adipose$Cellname_partial == "Adjacent parenchymal cells"] <- "Unknown"
adipose$Cellname_partial[adipose$Cellname_partial == "Endothelial cell"] <- "Endothelial_mixed"
adipose$Cellname_partial[adipose$Cellname_partial == "Macrophage/Monocyte/DC"] <- "Myeloid_mixed"
 
# Final merged categories used downstream.
adipose$Merged_cellname <- as.character(adipose$Cellname_partial)
map_cb_to_merged <- c(
  "Adipocyte (mature; lipid-signal mixed)" = "Adipocyte",
  "Adipocyte-like / lipid-signal dominated" = "Adipocyte",
  "Adipose stromal / Fibroblast" = "Adipose stromal cell",
  "B cell" = "B cell",
  "Endothelial (capillary~general)" = "Endothelial cell",
  "Erythroid" = "Erythroid cell",
  "Unassigned" = "Unknown",
  "Myeloid (Macrophage/DC; MHCII+; incl M2-like)" = "Macrophage / Monocyte / DC",
  "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)" = "VSMC/Pericyte"
)
idx_change <- as.character(adipose$Cellname_partial) %in% c("Endothelial_mixed", "Myeloid_mixed", "Unknown")
cb <- as.character(adipose$cluster_based_naming[idx_change])
new_label <- unname(map_cb_to_merged[cb])
new_label[is.na(new_label) | new_label == ""] <- "Unknown"
adipose$Merged_cellname[idx_change] <- new_label
saveRDS(adipose, "adipose_mergeobj_sketchHarmony_SCT_labelTransferred_k50_r020_Naming_260220.rds")
