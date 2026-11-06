# Usage from single_cell/: Rscript scripts/run_scenic.R RUN_ID SOURCE_ROOT OUTPUT_DIR DB_DIR
# SOURCE_ROOT is the original HFD_Aging_reanalysis directory (or its relocated copy).
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
suppressPackageStartupMessages(library(Seurat))
source("R/scenic_workflow.R")
cfg <- read.delim("config/scenic_execution.tsv", stringsAsFactors = FALSE)
cfg <- cfg[cfg$run_id == args[1], , drop = FALSE]
if (nrow(cfg) != 1L) stop("Select exactly one run_id from config/scenic_execution.tsv")
obj <- readRDS(file.path(args[2], cfg$input_relative_path))
if (nzchar(cfg$vessel_subset)) {
  obj <- obj[, rownames(obj[[]])[obj$Vessel_Type == cfg$vessel_subset]]
}
run_scenic_workflow(obj, args[3], args[4], cfg$run_id, cfg$group_by,
                    tissue = cfg$tissue_subset, ncores = cfg$ncores_original,
                    seed = cfg$seed_documented, ntrees = cfg$ntrees)
