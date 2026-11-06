# Source evidence; depends on the loaded/merged organ object.

# SOURCE Liver_normalize.ipynb CELL 13 (one-based)
liver_qc <- subset(
  liver,
  subset = nCount_Xenium > 5 & nFeature_Xenium > 3
)

table(liver$condition)
table(liver_qc$condition)

# SOURCE Liver_normalize.ipynb CELL 14 (one-based)
DefaultAssay(liver_qc) <- "Xenium"

liver_qc <- NormalizeData(liver_qc)

liver_qc <- ScaleData(
  liver_qc,
  features = rownames(liver_qc)
)

liver_qc <- RunPCA(
  liver_qc,
  features = rownames(liver_qc),
  npcs = 30
)

ElbowPlot(liver_qc, ndims = 30)

# SOURCE Liver_normalize.ipynb CELL 26 (one-based)
liver_qc$liver_cluster <- liver_qc$cluster_d15_res1.2
Idents(liver_qc) <- "liver_cluster"

DimPlot(
  liver_qc,
  reduction = "umap_d15",
  group.by = "liver_cluster",
  label = TRUE,
  repel = TRUE
) + NoLegend()
