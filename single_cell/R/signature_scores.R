# Source: IFNscore_AllECs_redgeplot.ipynb, cell 15 (one-based notebook cells).
signature_controls <- function(mat, genes, controls_per_gene = 10L) {
  stopifnot(!anyDuplicated(rownames(mat)), all(genes %in% rownames(mat)))
  stats <- data.frame(gene = rownames(mat),
                      mean = Matrix::rowMeans(mat),
                      zero = Matrix::rowMeans(mat == 0))
  stats$mean <- as.numeric(scale(stats$mean))
  stats$zero <- as.numeric(scale(stats$zero))
  if (anyNA(stats)) stop("Control-matching coordinates have zero variance or missing values")
  signature <- stats[stats$gene %in% genes, , drop = FALSE]
  background <- stats[!stats$gene %in% genes, , drop = FALSE]
  controls <- character()
  for (i in seq_len(nrow(signature))) {
    distance <- sqrt((background$mean - signature$mean[i])^2 +
                     (background$zero - signature$zero[i])^2)
    ordered <- background$gene[order(distance)]
    controls <- unique(c(controls, head(ordered[!ordered %in% controls], controls_per_gene)))
  }
  controls
}

calculate_signature_score <- function(mat, genes, controls_per_gene = 10L) {
  # Missing genes cause an explicit error instead of silently changing a signature.
  genes <- unique(genes)
  stopifnot(length(genes) > 0L, all(genes %in% rownames(mat)))
  controls <- signature_controls(mat, genes, controls_per_gene)
  stopifnot(length(controls) > 0L)
  score <- Matrix::colMeans(mat[genes, , drop = FALSE]) -
           Matrix::colMeans(mat[controls, , drop = FALSE])
  names(score) <- colnames(mat)
  attr(score, "control_genes") <- controls
  score
}

score_ifn_by_tissue <- function(object, ifn_genes) {
  stopifnot(all(c("Tissue", "Diet", "Condition", "Vessel_Type") %in% names(object[[]])))
  md <- object[[]]
  object <- object[, rownames(md)[md$Diet == "Normal"]]
  output <- lapply(c("Adipose", "Heart", "Muscle", "Liver"), function(tissue) {
    cells <- rownames(object[[]])[object$Tissue == tissue]
    if (!length(cells)) return(NULL)
    obj <- object[, cells]
    # RNA normalized data, as in the four original tissue-specific blocks.
    mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
    score <- calculate_signature_score(mat, ifn_genes)
    d <- obj[[]]
    d$genes.score <- score[rownames(d)]
    stopifnot(!anyNA(d$genes.score))
    d$cell_id <- rownames(d)
    d
  })
  do.call(rbind, output)
}

descriptive_cell_welch_tests <- function(score_metadata) {
  # Cell-level descriptive comparisons only. Cells are not independent animals.
  # Do not interpret these p-values as evidence of between-animal reproducibility.
  pieces <- split(score_metadata, list(score_metadata$Tissue,
                                       score_metadata$Vessel_Type), drop = TRUE)
  do.call(rbind, lapply(pieces, function(d) {
    d <- d[is.finite(d$genes.score) & !is.na(d$Condition), , drop = FALSE]
    n <- table(d$Condition)
    if (length(n) != 2L || any(n < 2L)) return(NULL)
    z <- t.test(genes.score ~ Condition, data = d)
    data.frame(Tissue = d$Tissue[1], Vessel_Type = d$Vessel_Type[1],
               p_value = z$p.value, n_cells = nrow(d),
               inference = "descriptive cell-level only")
  }))
}
