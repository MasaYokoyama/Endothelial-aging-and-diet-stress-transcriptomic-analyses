# Run from single_cell/: Rscript scripts/run_cellchat.R RUN_ID INPUT_RDS OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages({library(Seurat); library(CellChat)})
source("R/cellchat_workflow.R")
cfg <- read.delim("config/cellchat_runs.tsv", stringsAsFactors = FALSE)
cfg <- cfg[cfg$run_id == args[1], , drop = FALSE]
stopifnot(nrow(cfg) == 1L)
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(args[2])
chat <- run_cellchat_condition(obj, cfg$condition,
                               file.path(args[3], "cellchat_aggregated.rds"),
                               min_cells = cfg$min_cells)
# Source notebooks select K from diagnostics, then fit outgoing/incoming NMF.
# Stored K values are not automatically optimized or changed during publication.
ks <- as.integer(strsplit(cfg$outgoing_then_incoming_NMF_patterns, ";", fixed = TRUE)[[1]])
if (length(ks) == 2L && all(is.finite(ks))) {
  chat <- CellChat::identifyCommunicationPatterns(chat, pattern = "outgoing", k = ks[1])
  chat <- CellChat::identifyCommunicationPatterns(chat, pattern = "incoming", k = ks[2])
  saveRDS(chat, file.path(args[3], "cellchat_with_patterns.rds"))
}
writeLines(capture.output(sessionInfo()), file.path(args[3], "sessionInfo.txt"))
