# Source evidence; depends on the loaded/merged organ object.

# SOURCE Heart_normalize.ipynb CELL 14 (one-based)
heart_qc <- subset(
  heart,
  subset = nCount_Xenium > 5 & nFeature_Xenium > 3
)

table(heart$condition)
table(heart_qc$condition)

# SOURCE Heart_normalize.ipynb CELL 17 (one-based)
DefaultAssay(heart_qc) <- "Xenium"

heart_qc <- NormalizeData(heart_qc)

heart_qc <- ScaleData(
  heart_qc,
  features = rownames(heart_qc)
)

heart_qc <- RunPCA(
  heart_qc,
  features = rownames(heart_qc),
  npcs = 30
)

ElbowPlot(heart_qc, ndims = 30)

# SOURCE Heart_normalize.ipynb CELL 30 (one-based)
heart_qc$heart_cluster <- heart_qc$cluster_d20_res1
Idents(heart_qc) <- "heart_cluster"

DimPlot(
  heart_qc,
  reduction = "umap_d20",
  group.by = "heart_cluster",
  label = TRUE,
  repel = TRUE
) + NoLegend()
