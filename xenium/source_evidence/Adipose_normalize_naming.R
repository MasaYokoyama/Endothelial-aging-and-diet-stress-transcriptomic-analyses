# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Adipose_normalize_naming.ipynb CELL 8 (one-based)
DefaultAssay(adipose_qc) <- "Xenium"

adipose_cluster_annotation <- c(
  "0"  = "Endothelial-containing segment",
  "1"  = "Adipocyte-rich",
  "2"  = "Fibroblast-rich",
  "3"  = "Mesothelial cell-rich",
  "4"  = "Endothelial-containing segment",
  "5"  = "Adipocyte-rich",
  "6"  = "Pericyte-rich",
  "7"  = "Hematopoietic cell",
  "8"  = "Hematopoietic cell",
  "9"  = "Fibroblast-rich",
  "10" = "Endothelial-containing segment",
  "11" = "Adipocyte-rich",
  "12" = "Large-vessel EC-rich",
  "13" = "Hematopoietic cell",
  "14" = "VSMC-rich",
  "15" = "Fibroblast-rich",
  "16" = "Hematopoietic cell",
  "17" = "Mixed segment",
  "18" = "Endothelial-containing segment",
  "19" = "Hematopoietic cell",
  "20" = "Fibroblast-rich",
  "21" = "Lymphatic endothelial-rich",
  "22" = "Adipocyte-rich",
  "23" = "Adipocyte-rich"
)

adipose_qc$adipose_annotation <- unname(
  adipose_cluster_annotation[
    as.character(adipose_qc$adipose_cluster)
  ]
)

adipose_qc$adipose_annotation <- factor(
  adipose_qc$adipose_annotation
)

table(
  adipose_qc$adipose_cluster,
  adipose_qc$adipose_annotation
)


# SOURCE Adipose_normalize_naming.ipynb CELL 24 (one-based)
adipose_hema_annotation <- c(
  "0"  = "Monocyte/macrophage lineage",
  "1"  = "Monocyte/macrophage lineage",
  "2"  = "Macrophage-adipocyte mixed",
  "3"  = "T cell lineage",
  "4"  = "pDC-like",
  "5"  = "DC-associated",
  "6"  = "B cell lineage",
  "7"  = "Macrophage-mural/vascular mixed",
  "8"  = "Macrophage-endothelial mixed",
  "9"  = "Monocyte/macrophage lineage",
  "10" = "Macrophage-fibroblast mixed",
  "11" = "Macrophage-adipocyte mixed",
  "12" = "Macrophage-mural/vascular mixed",
  "13" = "B cell lineage",
  "14" = "Monocyte/macrophage lineage",
  "15" = "Macrophage-mesothelial mixed"
)

adipose_hematopoietic$hema_annotation <- unname(
  adipose_hema_annotation[
    as.character(adipose_hematopoietic$hema_cluster)
  ]
)

adipose_hematopoietic$hema_annotation <- factor(
  adipose_hematopoietic$hema_annotation
)

table(
  adipose_hematopoietic$hema_cluster,
  adipose_hematopoietic$hema_annotation
)


# SOURCE Adipose_normalize_naming.ipynb CELL 26 (one-based)
# 元のannotationをコピー
adipose_qc$adipose_annotation_refined <- as.character(
  adipose_qc$adipose_annotation
)

# hematopoietic subsetのcell名
hema_cells <- colnames(adipose_hematopoietic)

# Hematopoietic cellだった部分だけ置換
adipose_qc$adipose_annotation_refined[hema_cells] <- as.character(
  adipose_hematopoietic$hema_annotation
)

adipose_qc$adipose_annotation_refined <- factor(
  adipose_qc$adipose_annotation_refined
)

# 確認
table(
  before = adipose_qc$adipose_annotation,
  after = adipose_qc$adipose_annotation_refined
)
