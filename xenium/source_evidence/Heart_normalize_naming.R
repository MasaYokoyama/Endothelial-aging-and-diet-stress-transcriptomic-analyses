# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Heart_normalize_naming.ipynb CELL 17 (one-based)
DefaultAssay(heart_qc) <- "Xenium"

cluster_annotation <- c(
  "0"  = "Endothelial-containing segment",
  "1"  = "Ventricular cardiomyocyte-rich",
  "2"  = "Endothelial-containing segment",
  "3"  = "Ventricular cardiomyocyte-rich",
  "4"  = "Endothelial-containing segment",
  "5"  = "Pericyte-rich",
  "6"  = "Hematopoietic cell",
  "7"  = "Fibroblast-rich",
  "8"  = "Macrophage-EC mixed",
  "9"  = "Endothelial-containing segment",
  "10" = "Ventricular cardiomyocyte-rich",
  "11" = "Atrial cardiomyocyte-rich",
  "12" = "Stromal-EC mixed",
  "13" = "VSMC-rich",
  "14" = "Fibroblast-rich",
  "15" = "Macrophage-stromal mixed",
  "16" = "Venous-like large-vessel EC-rich",
  "17" = "Endothelial-containing segment",
  "18" = "Arterial-associated mural/VSMC-rich",
  "19" = "Venous-like large-vessel EC-rich",
  "20" = "Adipocyte-rich",
  "21" = "Hematopoietic cell",
  "22" = "Valve-associated fibroblast-rich",
  "23" = "Cardiomyocyte-glial mixed",
  "24" = "Epicardial/mesothelial-rich",
  "25" = "Lymphatic endothelial-rich",
  "26" = "Ventricular cardiomyocyte-rich",
  "27" = "Pericyte-rich",
  "28" = "Schwann/glial-rich",
  "29" = "Endothelial-containing segment"
)

heart_qc$heart_annotation <- unname(
  cluster_annotation[as.character(heart_qc$heart_cluster)]
)

ptprc_expression <- FetchData(
  heart_qc,
  vars = "Ptprc",
  layer = "data"
)[, 1]

hematopoietic_flag <-
  as.character(heart_qc$heart_cluster) %in% c("6", "21") |
  ptprc_expression >= 4

heart_qc$heart_annotation[hematopoietic_flag] <- "Hematopoietic cell"
heart_qc$heart_annotation <- factor(heart_qc$heart_annotation)

# cluster番号との対応を確認
table(
  cluster = heart_qc$heart_cluster,
  annotation = heart_qc$heart_annotation
)

# annotation後のUMAP
options(
  repr.plot.width = 14,
  repr.plot.height = 9,
  repr.plot.res = 120
)

DimPlot(
  heart_qc,
  reduction = "umap_d20",
  group.by = "heart_annotation",
  label = TRUE,
  repel = TRUE,
  label.size = 3
)


# SOURCE Heart_normalize_naming.ipynb CELL 27 (one-based)
DefaultAssay(heart_hematopoietic) <- "Xenium"

hema_expr <- FetchData(
  heart_hematopoietic,
  vars = c(
    "Cd3d", "Ms4a1",
    "Ppbp", "Csf3r",
    "Siglech", "Cd209a", "Itgax"
  ),
  layer = "data"
)

hema_cluster <- as.character(heart_hematopoietic$hema_cluster)

# 基本annotation
hema_annotation <- rep(
  "Monocyte/macrophage lineage",
  ncol(heart_hematopoietic)
)

names(hema_annotation) <- colnames(heart_hematopoietic)

# T/B lineage：cluster 3内だけ
is_cluster3 <- hema_cluster == "3"

hema_annotation[
  is_cluster3 &
  hema_expr$Ms4a1 >= 5.0
] <- "B cell lineage"

hema_annotation[
  is_cluster3 &
  hema_expr$Cd3d >= 4.5 &
  hema_expr$Ms4a1 < 5.0
] <- "T cell lineage"

# Neutrophil / Platelet：cluster 10・13内だけ
is_granulo_pool <- hema_cluster %in% c("10", "13")

hema_annotation[
  is_granulo_pool &
  hema_expr$Ppbp >= 6.0
] <- "Platelet-rich"

hema_annotation[
  is_granulo_pool &
  hema_expr$Ppbp < 6.0 &
  hema_expr$Csf3r >= 5.0
] <- "Neutrophil"

# DC-associated：cluster 5・12内だけ
is_dc_pool <- hema_cluster %in% c("5", "12")

hema_annotation[
  is_dc_pool &
  hema_expr$Itgax >= 4.5
] <- "DC-associated"

hema_annotation[
  is_dc_pool &
  hema_expr$Cd209a >= 5.0
] <- "moDC-like"

hema_annotation[
  is_dc_pool &
  hema_expr$Siglech >= 4.5 &
  hema_expr$Cd209a < 5.0
] <- "pDC-associated"

heart_hematopoietic$hema_annotation <- factor(hema_annotation)


# SOURCE Heart_normalize_naming.ipynb CELL 31 (one-based)
# 元のannotationをコピー
heart_qc$heart_annotation_refined <- as.character(
  heart_qc$heart_annotation
)

# hematopoietic subsetのcell名と新しいannotation
hema_cells <- colnames(heart_hematopoietic)

heart_qc$heart_annotation_refined[hema_cells] <- as.character(
  heart_hematopoietic$hema_annotation
)

heart_qc$heart_annotation_refined <- factor(
  heart_qc$heart_annotation_refined
)

# 置換結果を確認
table(
  before = heart_qc$heart_annotation,
  after = heart_qc$heart_annotation_refined
)

options(
  repr.plot.width = 14,
  repr.plot.height = 9,
  repr.plot.res = 120
)

DimPlot(
  heart_qc,
  reduction = "umap_d20",
  group.by = "heart_annotation_refined",
  label = TRUE,
  repel = TRUE,
  label.size = 3
)
