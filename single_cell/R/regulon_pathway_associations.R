# Preserve the two distinct analyses described in Reviewer 4 response 3.
# Inputs: TF-by-cell AUC matrix and cell-by-pathway score matrix.
regulon_pathway_associations <- function(auc, pathways, tfs,
                                         method = c("pearson", "spearman")) {
  method <- match.arg(method)
  stopifnot(!anyDuplicated(colnames(auc)), !anyDuplicated(rownames(pathways)),
            all(tfs %in% rownames(auc)),
            setequal(colnames(auc), rownames(pathways)))
  pathways <- as.matrix(pathways[colnames(auc), , drop = FALSE])
  out <- lapply(tfs, function(tf) {
    r <- vapply(seq_len(ncol(pathways)), function(j) {
      # Match the original complete-vector correlation: NA does not trigger
      # silent pairwise deletion; the association is retained as NA.
      stats::cor(as.numeric(auc[tf, ]), pathways[, j], method = method)
    }, numeric(1))
    data.frame(TF = tf, Pathway = colnames(pathways), Coefficient = r,
               method = method, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

rank_positive_pearson <- function(auc, pathways, tf = "Irf7") {
  x <- regulon_pathway_associations(auc, pathways, tf, "pearson")
  x <- x[!is.na(x$Coefficient) & x$Coefficient > 0, , drop = FALSE]
  x$R_squared <- x$Coefficient^2
  x[order(x$R_squared, decreasing = TRUE), , drop = FALSE]
}

spearman_network_edges <- function(auc, pathways,
                                   tfs = c("Irf7", "Irf9", "Stat1", "Stat2",
                                           "Junb", "Pparg", "Cebpd", "Jund", "Jun", "Klf2")) {
  regulon_pathway_associations(auc, pathways, tfs, "spearman")
}
# Correlations and positive Pearson r^2 rankings are exploratory cell-level
# associations. r^2 is not a significance test or an estimate of causation.
