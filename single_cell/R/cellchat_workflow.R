# Generalized from CellChat_EachCondition/CellChat_reanalysis_*.ipynb.
# This function constructs a condition-specific network. Plot selections and
# NMF pattern numbers remain in the per-condition manifest/source extracts.
run_cellchat_condition <- function(object, condition, output_file, min_cells) {
  stopifnot(requireNamespace("CellChat", quietly = TRUE))
  md <- object[[]]
  stopifnot(all(c("Celltype_CellChat", "Condition") %in% names(md)))
  mat <- Seurat::GetAssayData(object, assay = "SCT", slot = "counts")
  # As in the source: normalize SCT corrected counts using CellChat::normalizeData.
  mat <- CellChat::normalizeData(mat)
  cells <- rownames(md)[md$Condition == condition & !is.na(md$Condition)]
  stopifnot(length(cells) > 0L)
  meta <- data.frame(labels = md[cells, "Celltype_CellChat"],
                      sample = md[cells, "Condition"], row.names = cells)
  chat <- CellChat::createCellChat(object = mat[, cells, drop = FALSE],
                                  meta = meta, group.by = "labels")
  chat@DB <- CellChat::CellChatDB.mouse
  chat <- CellChat::subsetData(chat)
  chat <- CellChat::identifyOverExpressedGenes(chat)
  chat <- CellChat::identifyOverExpressedInteractions(chat)
  chat <- CellChat::projectData(chat, CellChat::PPI.mouse)
  chat <- CellChat::computeCommunProb(chat)
  # Original minimum group size differs: adipose 10, other organs 1.
  # A threshold of 1 retains nonempty groups; do not silently replace it with 10.
  chat <- CellChat::filterCommunication(chat, min.cells = min_cells)
  chat <- CellChat::computeCommunProbPathway(chat)
  chat <- CellChat::aggregateNet(chat)
  chat <- CellChat::netAnalysis_computeCentrality(chat, slot.name = "netP")
  saveRDS(chat, output_file)
  invisible(chat)
}
