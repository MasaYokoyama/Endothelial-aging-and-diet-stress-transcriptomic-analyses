# Run from single_cell/:
# Rscript scripts/run_junb_gene_degs.R INPUT_RDS RNA_OR_SCT OUTPUT_DIR [AGE_COLUMN] [AGED_LABEL]
# Defaults match the age labels used in the Supplementary Code: Type5, "1 year".
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L, length(args) <= 5L)
suppressPackageStartupMessages(library(Seurat))
source("R/junb_gene_degs.R")
gene_assay <- match.arg(args[2], c("RNA", "SCT"))
age_column <- if (length(args) >= 4L) args[4] else "Type5"
aged_label <- if (length(args) >= 5L) args[5] else "1 year"
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(args[1])
degs <- junb_gene_degs(obj, gene_assay = gene_assay, age_column = age_column,
                       aged_label = aged_label, cutoff = 0.15)
degs$gene <- rownames(degs)
write.csv(degs, file.path(args[3], "Junb_high_vs_low_gene_DEGs.csv"), row.names = FALSE)
ranked <- degs[!is.na(degs$avg_log2FC) & degs$avg_log2FC > 0, , drop = FALSE]
ranked <- ranked[order(ranked$avg_log2FC, decreasing = TRUE), , drop = FALSE]
write.csv(head(ranked, 100), file.path(args[3], "Junb_high_top100_by_log2FC.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[3], "sessionInfo.txt"))
