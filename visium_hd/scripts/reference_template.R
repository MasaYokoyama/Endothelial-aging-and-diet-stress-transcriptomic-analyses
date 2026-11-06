# Template retained from submitted methods; original adipose source notebook still required.
# Tissue-specific references used for spatial cell-type annotation.
# Public mouse aging atlas datasets were obtained from CZ CELLxGENE Discover
# as h5ad files, including Tabula Muris Senis heart, adipose tissue,
# and gonadal fat pad datasets. Cell category columns were converted to strings
# before conversion/merging into Seurat objects.

# Heart reference:
#   Tabula Muris Senis heart dataset + present-study heart single-cell RNA-seq dataset
#   -> normalization/integration
#   -> RCTD reference object used by spacexr
reference_rctd_heart <- readRDS("reference_rctd_heart.rds")

# Adipose reference inputs:
#   obj_wt      : previously analyzed adipose single-nucleus RNA-seq dataset
#                 deposited under BioProject accession PRJDB17809
#   seu_ageing  : Tabula Muris Senis adipose tissue reference
#   seu_gonadal : Tabula Muris Senis gonadal fat pad reference
#
# The public h5ad-derived RNA count matrices were converted to integer-like
# count matrices before Seurat reference construction.
counts_ageing <- round(GetAssayData(seu_ageing, assay = "RNA", layer = "counts"))
seu_ageing[["RNA"]] <- CreateAssayObject(counts = counts_ageing)

counts_gonadal <- round(GetAssayData(seu_gonadal, assay = "RNA", layer = "counts"))
seu_gonadal[["RNA"]] <- CreateAssayObject(counts = counts_gonadal)

# Harmonize source and cell-type metadata before merging.
obj_wt$sample      <- "snRNAseq"
seu_ageing$sample  <- "ageing10x"
seu_gonadal$sample <- "gonadal10x"

obj_wt$ref_celltype      <- as.character(obj_wt$celltype)
seu_ageing$ref_celltype  <- as.character(seu_ageing$cell_type)
seu_gonadal$ref_celltype <- as.character(seu_gonadal$cell_type)

DefaultAssay(obj_wt)      <- "RNA"
DefaultAssay(seu_ageing)  <- "RNA"
DefaultAssay(seu_gonadal) <- "RNA"

# Merge the three adipose reference datasets, normalize with SCTransform,
# and correct source-associated variation with Harmony.
merged_obj <- merge(
  x = obj_wt,
  y = list(seu_ageing, seu_gonadal),
  add.cell.ids = c("snRNAseq", "ageing10x", "gonadal10x"),
  project = "Adipose_ref_merged"
)

merged_obj <- SCTransform(
  merged_obj,
  method = "glmGamPoi",
  vars.to.regress = NULL,
  verbose = TRUE
)
merged_obj <- RunPCA(merged_obj, verbose = FALSE)
merged_obj <- RunHarmony(
  merged_obj,
  group.by.vars = "sample",
  verbose = TRUE
)
merged_obj <- RunUMAP(
  merged_obj,
  reduction = "harmony",
  dims = 1:30
)
merged_obj <- FindNeighbors(
  merged_obj,
  reduction = "harmony",
  dims = 1:30
)
merged_obj <- FindClusters(
  merged_obj,
  resolution = 0.5
)

# For final reference-label curation, the Harmony reduction was also examined
# using dimensions 1-20, with a dedicated neighbor graph and clustering.
dims_use <- 1:20
merged_obj <- RunUMAP(
  merged_obj,
  reduction = "harmony",
  dims = dims_use,
  reduction.name = "umap_harmony20"
)
merged_obj <- FindNeighbors(
  merged_obj,
  reduction = "harmony",
  dims = dims_use,
  graph.name = "harmony20_snn"
)
merged_obj <- FindClusters(
  merged_obj,
  graph.name = "harmony20_snn",
  resolution = 0.5,
  algorithm = 1
)

# Reference labels were reviewed using canonical marker expression and
# harmonized into the Cellname categories used for spatial label transfer.
merged_obj$Cellname <- merged_obj$ref_celltype
merged_obj$Cellname[merged_obj$seurat_clusters == "14"] <- "Neutrophil"

celltype_map <- c(
  "Adipocyte" = "Adipocyte",
  "ASC" = "Adipose stromal cell",
  "B cell" = "B cell",
  "CD4-positive, alpha-beta T cell" = "T cell",
  "CD8-positive, alpha-beta T cell" = "T cell",
  "endothelial cell" = "Endothelial cell",
  "Endothelial cell" = "Endothelial cell",
  "Epididymal cell" = "Adjacent parenchymal cells",
  "epithelial cell" = "Epithelial cell",
  "erythroblast" = "Erythroid cell",
  "Hematopoietic cell" = "Macrophage/Monocyte/DC",
  "lymphocyte" = "T cell",
  "mesenchymal stem cell of adipose tissue" = "Adipose stromal cell",
  "Mesothelial cell" = "Mesothelial cell",
  "myeloid cell" = "Macrophage/Monocyte/DC",
  "natural killer cell" = "NK cell",
  "Prostate cell" = "Adjacent parenchymal cells",
  "T cell" = "T cell",
  "VSMC" = "VSMC"
)

idx_other <- merged_obj$seurat_clusters != "14"
merged_obj$Cellname[idx_other] <-
  unname(celltype_map[merged_obj$ref_celltype[idx_other]])

# Save the final curated reference used for anchor-based label transfer.
saveRDS(merged_obj, file = "adipose_merged_obj_reference.rds")
adipose_reference <- readRDS("adipose_merged_obj_reference.rds")

# The heart reference was passed to create.RCTD()/run.RCTD().
# The adipose reference was passed to FindTransferAnchors()/TransferData().
