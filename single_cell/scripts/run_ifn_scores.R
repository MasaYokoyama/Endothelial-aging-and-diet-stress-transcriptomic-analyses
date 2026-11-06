# Usage from single_cell/: Rscript scripts/run_ifn_scores.R ENDOTHELIAL_RDS OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages(library(Seurat))
source("R/signature_scores.R")
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
genes <- read.delim("data/gene_sets.tsv", stringsAsFactors = FALSE)
genes <- genes$gene[genes$signature == "IFN65"]
object <- readRDS(args[1])
scores <- score_ifn_by_tissue(object, genes)
write.csv(scores, file.path(args[2], "IFN_scores_by_cell.csv"), row.names = FALSE)
tests <- descriptive_cell_welch_tests(scores)
write.csv(tests, file.path(args[2], "IFN_descriptive_cell_tests.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[2], "sessionInfo.txt"))
