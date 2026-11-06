# Source evidence; depends on the loaded/merged organ object.

# SOURCE Muscle_normalize.ipynb CELL 13 (one-based)
muscle_qc <- subset(
  muscle,
  subset = nCount_Xenium > 5 & nFeature_Xenium > 3
)

table(muscle$condition)
table(muscle_qc$condition)

# SOURCE Muscle_normalize.ipynb CELL 16 (one-based)
DefaultAssay(muscle_qc) <- "Xenium"

muscle_qc <- NormalizeData(muscle_qc)

muscle_qc <- ScaleData(
  muscle_qc,
  features = rownames(muscle_qc)
)

muscle_qc <- RunPCA(
  muscle_qc,
  features = rownames(muscle_qc),
  npcs = 30
)

ElbowPlot(muscle_qc, ndims = 30)

# SOURCE Muscle_normalize.ipynb CELL 28 (one-based)
muscle_qc$muscle_cluster <- muscle_qc$cluster_d20_res1.2
Idents(muscle_qc) <- "muscle_cluster"
