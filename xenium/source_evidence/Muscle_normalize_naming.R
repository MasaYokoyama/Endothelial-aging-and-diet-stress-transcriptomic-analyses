# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Muscle_normalize_naming.ipynb CELL 14 (one-based)
muscle_cluster_annotation <- c(
  "0"  = "Skeletal muscle-rich",
  "1"  = "Endothelial-containing segment",
  "2"  = "Adventitial fibroblast-rich",
  "3"  = "Endothelial-containing segment",
  "4"  = "Skeletal muscle-rich",
  "5"  = "Skeletal muscle-rich",
  "6"  = "Skeletal muscle-rich",
  "7"  = "Skeletal muscle-rich",
  "8"  = "Hematopoietic cell",
  "9"  = "Pericyte-rich",
  "10" = "ECM-rich fibroblast",
  "11" = "Skeletal muscle-rich",
  "12" = "Skeletal muscle-rich",
  "13" = "Adventitial fibroblast-rich",
  "14" = "Endothelial-containing segment",
  "15" = "Large-vessel EC-rich",
  "16" = "Hematopoietic cell",
  "17" = "Hematopoietic cell",
  "18" = "Schwann/glial-rich",
  "19" = "Skeletal muscle-epithelial mixed",
  "20" = "VSMC-rich",
  "21" = "ECM-rich fibroblast",
  "22" = "Muscle spindle-associated muscle-rich",
  "23" = "Mesothelial/epithelial-rich",
  "24" = "Adipocyte-rich",
  "25" = "Skeletal muscle-rich",
  "26" = "Lymphatic endothelial-rich",
  "27" = "Skeletal muscle-rich"
)

muscle_qc$muscle_annotation <- unname(
  muscle_cluster_annotation[
    as.character(muscle_qc$muscle_cluster)
  ]
)

muscle_qc$muscle_annotation <- factor(
  muscle_qc$muscle_annotation
)

table(
  muscle_qc$muscle_cluster,
  muscle_qc$muscle_annotation
)

options(
  repr.plot.width = 8,
  repr.plot.height = 6,
  repr.plot.res = 140
)

DimPlot(
  muscle_qc,
  reduction = "umap_d20",
  group.by = "muscle_annotation",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + NoLegend()


# SOURCE Muscle_normalize_naming.ipynb CELL 22 (one-based)
muscle_hema_annotation <- c(
  "0"  = "Monocyte/macrophage lineage",
  "1"  = "Macrophage-fibroblast mixed",
  "2"  = "Macrophage-muscle mixed",
  "3"  = "Macrophage-endothelial mixed",
  "4"  = "Monocyte/macrophage lineage",
  "5"  = "Macrophage-muscle mixed",
  "6"  = "Macrophage-muscle mixed",
  "7"  = "Trem2-enriched macrophage",
  "8"  = "DC-associated",
  "9"  = "Macrophage-muscle mixed",
  "10" = "Hematopoietic-associated mixed",
  "11" = "Macrophage-vascular/mural mixed",
  "12" = "Macrophage-adipocyte mixed"
)

muscle_hematopoietic$hema_annotation <- unname(
  muscle_hema_annotation[
    as.character(muscle_hematopoietic$hema_res0.8)
  ]
)

muscle_hematopoietic$hema_annotation <- factor(
  muscle_hematopoietic$hema_annotation
)

table(
  muscle_hematopoietic$hema_res0.8,
  muscle_hematopoietic$hema_annotation
)


# SOURCE Muscle_normalize_naming.ipynb CELL 23 (one-based)
tb_expr <- FetchData(
  muscle_hematopoietic,
  vars = c("Cd3d", "Ms4a1"),
  layer = "data"
)

summary(tb_expr$Cd3d)
summary(tb_expr$Ms4a1)


# SOURCE Muscle_normalize_naming.ipynb CELL 24 (one-based)
muscle_hematopoietic$hema_annotation <-
  as.character(muscle_hematopoietic$hema_annotation)

muscle_hematopoietic$hema_annotation[
  tb_expr$Cd3d >= 5.0
] <- "T cell lineage"

muscle_hematopoietic$hema_annotation[
  tb_expr$Ms4a1 >= 5.5
] <- "B cell lineage"

muscle_hematopoietic$hema_annotation <- factor(
  muscle_hematopoietic$hema_annotation
)


# SOURCE Muscle_normalize_naming.ipynb CELL 25 (one-based)
muscle_hematopoietic$hema_annotation <- as.character(
  muscle_hematopoietic$hema_annotation
)

muscle_hematopoietic$hema_annotation[
  muscle_hematopoietic$hema_annotation == "Hematopoietic-associated mixed"
] <- "B cell-associated mixed"

muscle_hematopoietic$hema_annotation <- factor(
  muscle_hematopoietic$hema_annotation
)


# SOURCE Muscle_normalize_naming.ipynb CELL 32 (one-based)
muscle_qc$muscle_annotation_refined <- as.character(
  muscle_qc$muscle_annotation
)

muscle_qc$muscle_annotation_refined[
  colnames(muscle_hematopoietic)
] <- as.character(
  muscle_hematopoietic$hema_annotation
)

present_labels <- unique(
  muscle_qc$muscle_annotation_refined
)

final_order <- c(
  intersect(naming1_order, present_labels),
  setdiff(present_labels, naming1_order)
)

muscle_qc$muscle_annotation_refined <- factor(
  muscle_qc$muscle_annotation_refined,
  levels = final_order
)

DimPlot(
  muscle_qc,
  reduction = umap_name,
  group.by = "muscle_annotation_refined",
  label = FALSE,
  raster = TRUE
)
