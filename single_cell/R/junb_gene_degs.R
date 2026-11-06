# Classify cells by Junb regulon activity and test gene expression in the
# explicitly selected expression assay.
junb_gene_degs <- function(object, gene_assay, age_column = "Type5",
                           aged_label = "1 year", cutoff = 0.15) {
  if (missing(gene_assay)) stop("Specify the gene expression assay: RNA or SCT")
  stopifnot(gene_assay %in% c("RNA", "SCT"), gene_assay %in% names(object@assays))
  stopifnot(age_column %in% names(object[[]]))
  auc <- Seurat::GetAssayData(object, assay = "regulonAUC_extended", slot = "data")
  stopifnot("Junb" %in% rownames(auc), setequal(colnames(auc), colnames(object)))
  # 0.15 was selected as the midpoint of the activity range in Extended Data
  # Fig. 4j (author response), to separate relatively high and low Junb activity.
  object$Regulon_Junb <- ifelse(auc["Junb", colnames(object)] >= cutoff,
                               "Junb high", "Junb low")
  age_values <- object[[]][[age_column]]
  cells <- rownames(object[[]])[!is.na(age_values) & age_values == aged_label]
  stopifnot(length(cells) > 0L)
  aged <- object[, cells]
  # Explicitly use a GENE expression assay. Never infer this from DefaultAssay
  # after a regulon visualization; regulonAUC_extended contains TF activities.
  Seurat::FindMarkers(aged, assay = gene_assay, group.by = "Regulon_Junb",
                       ident.1 = "Junb high", ident.2 = "Junb low",
                       logfc.threshold = 0.1)
}
