# Source evidence; depends on the loaded/merged organ object.

# SOURCE Adipose_normalize.ipynb CELL 14 (one-based)
adipose_qc <- subset(
  adipose,
  subset = nCount_Xenium > 5 & nFeature_Xenium > 3
)

table(adipose$condition)
table(adipose_qc$condition)

# SOURCE Adipose_normalize.ipynb CELL 16 (one-based)
DefaultAssay(adipose_qc) <- "Xenium"

adipose_qc <- NormalizeData(adipose_qc)

adipose_qc <- ScaleData(
  adipose_qc,
  features = rownames(adipose_qc)
)

adipose_qc <- RunPCA(
  adipose_qc,
  features = rownames(adipose_qc),
  npcs = 30
)

ElbowPlot(adipose_qc, ndims = 30)

# SOURCE Adipose_normalize.ipynb CELL 27 (one-based)
adipose_qc$adipose_cluster <- adipose_qc$cluster_d15_res1.2
Idents(adipose_qc) <- "adipose_cluster"

DimPlot(
  adipose_qc,
  reduction = "umap_d15",
  group.by = "adipose_cluster",
  label = TRUE,
  repel = TRUE
) + NoLegend()
