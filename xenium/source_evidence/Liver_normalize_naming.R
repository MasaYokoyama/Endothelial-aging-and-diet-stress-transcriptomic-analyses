# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Liver_normalize_naming.ipynb CELL 16 (one-based)
liver_cluster_annotation <- c(
  "0"  = "Hepatocyte-rich",
  "1"  = "Hepatocyte-rich",
  "2"  = "Hepatocyte-rich",
  "3"  = "Hepatocyte-rich",
  "4"  = "Hepatocyte-rich",
  "5"  = "LSEC-rich",
  "6"  = "Kupffer cell",
  "7"  = "Hepatocyte-rich",
  "8"  = "LSEC-rich",
  "9"  = "Hepatic stellate cell-rich",
  "10" = "Hepatocyte-rich",
  "11" = "Kupffer-LSEC mixed",
  "12" = "Hepatocyte-rich",
  "13" = "Hepatocyte-rich",
  "14" = "Hepatocyte-rich",
  "15" = "Cholangiocyte",
  "16" = "Hematopoietic cell",
  "17" = "Hepatic stellate cell-rich",
  "18" = "Large-vessel EC-rich",
  "19" = "Hematopoietic cell",
  "20" = "Endothelial-containing mixed",
  "21" = "Mesothelial-hepatocyte mixed",
  "22" = "VSMC-rich",
  "23" = "IFN-responsive hepatocyte-rich",
  "24" = "Mesothelial cell"
)

liver_qc$liver_annotation <- unname(
  liver_cluster_annotation[
    as.character(liver_qc$liver_cluster)
  ]
)

liver_qc$liver_annotation <- factor(
  liver_qc$liver_annotation
)


# SOURCE Liver_normalize_naming.ipynb CELL 35 (one-based)
liver_hema_annotation <- c(
  "0"  = "Hepatic stellate cell-rich",
  "1"  = "T cell lineage",
  "2"  = "T cell lineage",
  "3"  = "Macrophage lineage",
  "4"  = "B cell lineage",
  "5"  = "DC-associated",
  "6"  = "DC-associated",
  "7"  = "T cell lineage",
  "8"  = "Macrophage lineage",
  "9"  = "Macrophage lineage",
  "10" = "pDC-like",
  "11" = "Neutrophil",
  "12" = "Lymphatic endothelial-rich",
  "13" = "Macrophage lineage",
  "14" = "Hepatocyte-rich",
  "15" = "Schwann/glial-rich",
  "16" = "Portal/adventitial fibroblast-rich"
)

liver_hematopoietic$hema_annotation <- unname(
  liver_hema_annotation[
    as.character(liver_hematopoietic$hema_res0.8)
  ]
)

liver_hematopoietic$hema_annotation <- as.character(
  liver_hematopoietic$hema_annotation
)

t_cells <- rownames(tb_expr)[tb_expr$Cd3d >= 4.0]
b_cells <- rownames(tb_expr)[tb_expr$Ms4a1 >= 4.5]

liver_hematopoietic$hema_annotation <- as.character(
  liver_hematopoietic$hema_annotation
)

liver_hematopoietic$hema_annotation[
  colnames(liver_hematopoietic) %in% t_cells
] <- "T cell lineage"

liver_hematopoietic$hema_annotation[
  colnames(liver_hematopoietic) %in% b_cells
] <- "B cell lineage"

liver_hematopoietic$hema_annotation <- factor(
  liver_hematopoietic$hema_annotation
)


# SOURCE Liver_normalize_naming.ipynb CELL 41 (one-based)
t_cells <- rownames(tb_expr)[tb_expr$Cd3d >= 3.0]
b_cells <- rownames(tb_expr)[tb_expr$Ms4a1 >= 4.0]

liver_hematopoietic$hema_annotation <- as.character(
  liver_hematopoietic$hema_annotation
)

liver_hematopoietic$hema_annotation[
  colnames(liver_hematopoietic) %in% t_cells
] <- "T cell lineage"

liver_hematopoietic$hema_annotation[
  colnames(liver_hematopoietic) %in% b_cells
] <- "B cell lineage"

liver_hematopoietic$hema_annotation <- factor(
  liver_hematopoietic$hema_annotation
)


# SOURCE Liver_normalize_naming.ipynb CELL 46 (one-based)
liver_qc$liver_annotation_refined <- as.character(
  liver_qc$liver_annotation
)

liver_qc$liver_annotation_refined[
  colnames(liver_hematopoietic)
] <- as.character(
  liver_hematopoietic$hema_annotation
)

liver_qc$liver_annotation_refined <- factor(
  liver_qc$liver_annotation_refined
)


# SOURCE Liver_normalize_naming.ipynb CELL 55 (one-based)
liver_hematopoietic$hema_annotation <- as.character(
  liver_hematopoietic$hema_annotation
)

liver_hematopoietic$hema_annotation[
  liver_hematopoietic$hema_cluster %in% c("0", "16")
] <- "Fibroblast-rich"

liver_hematopoietic$hema_annotation <- factor(
  liver_hematopoietic$hema_annotation
)


# SOURCE Liver_normalize_naming.ipynb CELL 56 (one-based)
liver_qc$liver_annotation_refined <- as.character(
  liver_qc$liver_annotation_refined
)

common_cells <- intersect(
  colnames(liver_qc),
  colnames(liver_hematopoietic)
)

liver_qc$liver_annotation_refined[common_cells] <-
  as.character(liver_hematopoietic$hema_annotation)[
    match(common_cells, colnames(liver_hematopoietic))
  ]

liver_qc$liver_annotation_refined <- factor(
  liver_qc$liver_annotation_refined
)
