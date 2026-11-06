# Run from single_cell/: Rscript scripts/run_pathway_correlations.R INPUT_RDS REACTOME_CSV OUTPUT_DIR
# REACTOME_CSV must contain column x with the pathway score metadata column names.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages(library(Seurat))
source("R/regulon_pathway_associations.R")
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(args[1])
paths <- read.csv(args[2], stringsAsFactors = FALSE)$x
stopifnot(length(paths) > 0L, !anyDuplicated(paths), all(paths %in% names(obj[[]])))
auc <- GetAssayData(obj, assay = "regulonAUC_extended", slot = "data")
scores <- obj[[]][, paths, drop = FALSE]
pearson <- rank_positive_pearson(auc, scores, tf = "Irf7")
spearman <- spearman_network_edges(auc, scores)
write.csv(pearson, file.path(args[3], "Irf7_positive_Pearson_r_squared.csv"), row.names = FALSE)
write.csv(spearman, file.path(args[3], "multi_TF_Spearman_coefficients.csv"), row.names = FALSE)
top20 <- do.call(rbind, lapply(split(spearman, spearman$TF), function(x) {
  x <- x[!is.na(x$Coefficient), , drop = FALSE]
  head(x[order(x$Coefficient, decreasing = TRUE), , drop = FALSE], 20)
}))
write.csv(top20, file.path(args[3], "multi_TF_top20_network_edges.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[3], "sessionInfo.txt"))
