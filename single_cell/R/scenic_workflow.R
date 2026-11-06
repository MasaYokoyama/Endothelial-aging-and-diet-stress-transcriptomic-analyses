# Generalized from the source SCENIC notebooks; see provenance/code_map.tsv.
# Requires the historical SCENIC/RcisTarget/AUCell environment and v9 motif files.
run_scenic_workflow <- function(object, output_dir, db_dir, dataset_title,
                                group_by, tissue = NULL, ncores = 4L,
                                seed = 123L, ntrees = 500L) {
  stopifnot(requireNamespace("Seurat", quietly = TRUE),
            requireNamespace("SCENIC", quietly = TRUE),
            requireNamespace("RcisTarget", quietly = TRUE),
            requireNamespace("AUCell", quietly = TRUE))
  if (!is.null(tissue) && tissue != "ALL") {
    stopifnot("Tissue" %in% colnames(object[[]]))
    object <- object[, rownames(object[[]])[object$Tissue == tissue]]
  }
  stopifnot(ncol(object) > 0L, group_by %in% colnames(object[[]]))
  # SCT counts are depth-corrected UMI counts, not Pearson residuals.
  # Use the SCT count matrix used in the original study.
  expr <- as.matrix(Seurat::GetAssayData(object, assay = "SCT", slot = "counts"))
  stopifnot(identical(colnames(expr), rownames(object[[]])))
  dbs <- c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
           "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
  db_dir <- normalizePath(db_dir, mustWork = TRUE)
  stopifnot(all(file.exists(file.path(db_dir, dbs))))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE); setwd(output_dir)
  dir.create("int", showWarnings = FALSE); dir.create("output", showWarnings = FALSE)
  saveRDS(object[[]], "int/cellInfo.Rds")
  # Historical RcisTarget v9 annotations used by the source notebooks.
  data("motifAnnotations_mgi_v9", package = "RcisTarget", envir = .GlobalEnv)
  scenic <- SCENIC::initializeScenic(org = "mgi", dbDir = db_dir, dbs = dbs,
                                    datasetTitle = dataset_title, nCores = ncores)
  scenic@settings$seed <- seed
  set.seed(seed)
  # Original expression filter: counts >= 3 * 1% of cells and detected in 1%.
  genes <- SCENIC::geneFiltering(expr, scenicOptions = scenic,
                                 minCountsPerGene = 3 * .01 * ncol(expr),
                                 minSamples = .01 * ncol(expr))
  filtered <- expr[genes, , drop = FALSE]
  # Correlation on corrected counts precedes the log transform in the notebooks.
  SCENIC::runCorrelation(filtered, scenic)
  SCENIC::runGenie3(log2(filtered + 1), scenic, nTrees = ntrees)
  scenic <- SCENIC::runSCENIC_1_coexNetwork2modules(scenic)
  scenic <- SCENIC::runSCENIC_2_createRegulons(scenic)
  # AUCell uses the full gene matrix, not only the genes retained for GENIE3.
  scenic <- SCENIC::runSCENIC_3_scoreCells(scenic, log2(expr + 1))
  auc_all <- AUCell::getAUC(SCENIC::loadInt(scenic, "aucell_regulonAUC"))
  auc <- auc_all[SCENIC::onlyNonDuplicatedExtended(rownames(auc_all)), , drop = FALSE]
  groups <- object[[]][colnames(auc), group_by]
  names(groups) <- colnames(auc)
  rss <- SCENIC::calcRSS(AUC = auc, cellAnnotation = groups)
  # Preserve the separately reported extended-only assay/RSS, as in original
  # notebook cells 46–56; do not confuse it with onlyNonDuplicatedExtended.
  extended <- auc_all[grep("_extended", rownames(auc_all)), , drop = FALSE]
  rownames(extended) <- gsub(" \\(\\d+g\\)", "", rownames(extended))
  rownames(extended) <- gsub("_extended", "", rownames(extended))
  if (nrow(extended)) {
    stopifnot(!anyDuplicated(rownames(extended)))
    object[["regulonAUC_extended"]] <- Seurat::CreateAssayObject(counts = extended)
    rss_extended <- SCENIC::calcRSS(AUC = extended, cellAnnotation = groups)
    write.csv(rss_extended, "regulon_rss_extended.csv")
    write.csv(extended, "regulon_auc_extended.csv")
  }
  saveRDS(scenic, "scenicOptions_final.rds")
  saveRDS(auc, "regulon_auc.rds"); saveRDS(rss, "regulon_rss.rds")
  saveRDS(auc_all, "regulon_auc_all.rds")
  saveRDS(object, "seurat_with_regulon_auc_extended.rds")
  write.csv(rss, "regulon_rss.csv")
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  invisible(list(options = scenic, auc = auc, rss = rss))
}
