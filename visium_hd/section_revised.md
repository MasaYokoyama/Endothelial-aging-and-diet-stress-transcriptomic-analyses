# Part II. Spatial transcriptomics analyses

Heart Visium HD samples use 8-µm bins; adipose samples use segmented cells.
The detailed `visium_hd/` folder provides all eight sample settings, input object
specifications, region-selection tables, and executable neighborhood/LIANA scripts.

## 1. R setup and sample definitions


```r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(Matrix)
  library(hdf5r)
  library(arrow)
  library(harmony)
  library(dplyr)
  library(ggplot2)
})
 
heart_outs <- list(
  Heart_Aged_1  = "/work/Visium/Visium_session1/Heart_Aged_1/outs",
  Heart_Young_1 = "/work/Visium/Visium_session1/Heart_Young_1/outs",
  Heart_Aged_2  = "/work/Visium/Visium_session2/Heart_Aged_2/outs",
  Heart_Young_2 = "/work/Visium/Visium_session2/Heart_Young_2/outs"
)
 
adipose_outs <- list(
  Adipose_Aged_1  = "/work/Visium/Visium_session1/Adipose_Aged_1/outs",
  Adipose_Young_1 = "/work/Visium/Visium_session1/Adipose_Young_1/outs",
  Adipose_Aged_2  = "/work/Visium/Visium_session2/Adipose_Aged_2/outs",
  Adipose_Young_2 = "/work/Visium/Visium_session2/Adipose_Young_2/outs"
)
```


## 2. Reference construction for spatial cell-type annotation


```r
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
```


## 3. Heart Visium HD 8-um bin processing in Seurat


```r
# Input: Space Ranger v4.0.1 binned output, square_008um/filtered_feature_bc_matrix.h5
h5 <- lapply(heart_outs, function(o) {
  file.path(o, "binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
})
stopifnot(all(file.exists(unlist(h5))))
 
odir <- "on_disk_mats_heart"
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
 
mat_list <- lapply(names(h5), function(nm) {
  m <- BPCells::open_matrix_10x_hdf5(path = h5[[nm]], feature_type = "Gene Expression")
  out <- file.path(odir, nm)
  BPCells::write_matrix_dir(mat = m, dir = out, overwrite = FALSE)
  BPCells::open_matrix_dir(out)
})
names(mat_list) <- names(h5)
 
obj_list <- lapply(names(mat_list), function(nm) {
  so <- CreateSeuratObject(counts = mat_list[[nm]], project = nm, assay = "RNA",
                           min.cells = 0, min.features = 0)
  so$sample <- nm
  RenameCells(so, add.cell.id = nm)
})
names(obj_list) <- names(mat_list)
 
heart <- Reduce(function(a, b) merge(a, b), obj_list)
# Initial permissive QC removes very low-count bins before sketching.
heart <- subset(heart, subset = nCount_RNA >= 25)
 
DefaultAssay(heart) <- "RNA"
heart <- NormalizeData(heart, scale.factor = median(heart$nCount_RNA), verbose = FALSE)
heart <- FindVariableFeatures(heart, verbose = FALSE)
heart <- SketchData(heart, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
 
DefaultAssay(heart) <- "sketch"
heart <- NormalizeData(heart, verbose = FALSE)
heart <- FindVariableFeatures(heart, verbose = FALSE)
heart <- ScaleData(heart, verbose = FALSE)
heart <- RunPCA(heart, npcs = 50, verbose = FALSE)
heart <- IntegrateLayers(heart, method = HarmonyIntegration, orig = "pca",
                         new.reduction = "harmony", dims = 1:30)
heart <- RunUMAP(heart, reduction = "harmony", dims = 1:30, reduction.name = "umap")
heart <- FindNeighbors(heart, reduction = "harmony", dims = 1:30)
heart <- FindClusters(heart, resolution = 0.5)
 
heart <- ProjectIntegration(heart, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
heart <- ProjectData(heart, assay = "RNA", full.reduction = "harmony.full",
                     sketched.assay = "sketch", sketched.reduction = "harmony.full",
                     umap.model = "umap", dims = 1:30,
                     refdata = list(cluster_full = "seurat_clusters"))
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony.rds")
```


## 4. Heart gene-symbol conversion and RCTD-based annotation


```r
# For the Seurat object used for visualization and marker-based annotation,
# Ensembl identifiers were renamed to MGI gene symbols using biomaRt.
# Unmapped Ensembl identifiers were retained.
convert_ensembl_to_symbol <- function(seurat_obj, assay = "RNA") {
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  ens_ids <- rownames(seurat_obj[[assay]])
  annot <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                 filters = "ensembl_gene_id", values = ens_ids, mart = mart)
  annot <- annot[annot$mgi_symbol != "", ]
  annot <- annot[!duplicated(annot$ensembl_gene_id), ]
  rownames(seurat_obj[[assay]]) <- make.unique(ifelse(
    ens_ids %in% annot$ensembl_gene_id,
    annot$mgi_symbol[match(ens_ids, annot$ensembl_gene_id)],
    ens_ids
  ))
  seurat_obj
}
 
heart <- readRDS("heart_mergeobj_preProject_full_sketchHarmony.rds")
heart <- convert_ensembl_to_symbol(heart, assay = "RNA")
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept.rds")
 
# Cluster marker-based naming
DefaultAssay(heart) <- "RNA"
Idents(heart) <- heart$cluster_full
heart <- JoinLayers(heart)
markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 
heart$cluster_naming <- NA_character_
cl <- as.character(heart$cluster_full)
heart$cluster_naming[cl == "0"]  <- "Atrial cardiac muscle"
heart$cluster_naming[cl == "16"] <- "Platelet"
heart$cluster_naming[cl == "10"] <- "Red blood cell"
heart$cluster_naming[cl == "23"] <- "Airway epithelium"
heart$cluster_naming[cl == "17"] <- "Airway structural tissue"
heart$cluster_naming[cl == "26"] <- "Cardiac neuron"
npr3 <- FetchData(heart, vars = "Npr3", layer = "data")[colnames(heart), 1]
# Final saved rule: cluster 8 with detectable Npr3 at the selected 0.01 cutoff.
heart$cluster_naming[cl == "8" & npr3 >= 0.01] <- "Endocardium"
heart$cluster_naming[cl == "12"] <- "Epicardium"
heart$cluster_naming[cl == "5"]  <- "Adipose tissue"
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming.rds")
 
# RCTD annotation per sample
library(spacexr)
library(org.Mm.eg.db)
reference_rctd <- readRDS("reference_rctd_heart.rds")
heart_raw <- readRDS("heart_mergeobj_preProject_full_sketchHarmony.rds")
heart_raw <- JoinLayers(heart_raw, assay = "RNA", layers = "counts", new = "counts_joined")
counts <- GetAssayData(heart_raw, assay = "RNA", layer = "counts_joined")
 
run_rctd_one_heart <- function(sample, pos_path, nCount_cutoff = 30) {
  # Pre-mapping total-count QC. RCTD UMI_min independently filters mapped counts.
  cells <- colnames(heart_raw)[heart_raw$sample == sample & heart_raw$nCount_RNA >= nCount_cutoff]
  pos <- read_parquet(pos_path)
  coords_raw <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
                           row.names = pos$barcode)
  cells_base <- sub("^.*_(s_008um_)", "\\1", cells)
  idx <- match(cells_base, rownames(coords_raw))
  stopifnot(all(!is.na(idx)))
  coords_sample <- coords_raw[idx, , drop = FALSE]
  rownames(coords_sample) <- cells
 
  counts_use <- counts[, cells, drop = FALSE]

  # RCTD input was constructed independently from the raw count matrix.
  # Here Ensembl IDs were mapped with org.Mm.eg.db, and counts mapping to the
  # same gene symbol were summed before SpatialRNA object construction.
  # This serves a different purpose from the biomaRt-based renaming above.
  ens_keys <- unique(rownames(counts_use))
  map_df <- AnnotationDbi::select(org.Mm.eg.db, keys = ens_keys, keytype = "ENSEMBL",
                                  columns = "SYMBOL")
  map_df <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", c("ENSEMBL", "SYMBOL")]
  map_df <- map_df[!duplicated(map_df$ENSEMBL), ]
  ens2sym <- setNames(map_df$SYMBOL, map_df$ENSEMBL)
  ens_in <- intersect(rownames(counts_use), names(ens2sym))
  M_ens <- counts_use[ens_in, , drop = FALSE]
  sym_ids <- ens2sym[ens_in]
  symbols <- unique(sym_ids)
  A <- sparseMatrix(i = match(sym_ids, symbols), j = seq_along(sym_ids), x = 1,
                    dims = c(length(symbols), length(sym_ids)))
  M_sym <- A %*% M_ens
  rownames(M_sym) <- symbols
 
  spatialRNA <- SpatialRNA(coords = coords_sample, counts = M_sym, nUMI = Matrix::colSums(M_sym))
  myRCTD <- create.RCTD(spatialRNA = spatialRNA, reference = reference_rctd,
                        max_cores = 8, keep_reference = TRUE, UMI_min = 30)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  saveRDS(myRCTD, file = paste0("RCTD_fullref_cutoff30_", sample, ".rds"))
}
```


## 5. Heart merged cell labels, XY coordinates, and quality filtering


```r
# Merge RCTD predictions, treat rejected spots as Unknown, and attach to Seurat metadata.
paths <- c(
  Heart_Young_1 = "RCTD_fullref_cutoff30_Heart_Young_1.rds",
  Heart_Young_2 = "RCTD_fullref_cutoff30_Heart_Young_2.rds",
  Heart_Aged_1  = "RCTD_fullref_cutoff30_Heart_Aged_1.rds",
  Heart_Aged_2  = "RCTD_fullref_cutoff30_Heart_Aged_2.rds"
)
 
read_one_rctd <- function(sample, fp) {
  r <- readRDS(fp)
  cellid <- colnames(r@spatialRNA@counts)
  W <- r@results$weights
  results_df <- r@results$results_df
  stopifnot(all(cellid %in% rownames(W)), all(cellid %in% rownames(results_df)))
  W <- W[cellid, , drop = FALSE]
  pred_type <- colnames(W)[max.col(W, ties.method = "first")]
  spot_class <- as.character(results_df[cellid, "spot_class"])
  pred_type2 <- pred_type
  pred_type2[spot_class == "reject"] <- "Unknown"
  data.frame(cellid = cellid, sample = sample, spot_class = spot_class,
             rctd_type = pred_type, rctd_type_unknown = pred_type2,
             stringsAsFactors = FALSE)
}
 
rctd_merged <- do.call(rbind, Map(read_one_rctd, names(paths), paths))
stopifnot(!anyDuplicated(rctd_merged$cellid))
saveRDS(rctd_merged, "RCTD_merged_cutoff30_rejectUnknown.rds")
 
obj <- readRDS("heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming.rds")
spotclass_vec <- setNames(rctd_merged$spot_class, rctd_merged$cellid)
type_vec <- setNames(rctd_merged$rctd_type_unknown, rctd_merged$cellid)
obj$rctd_spot_class <- unname(spotclass_vec[colnames(obj)])
obj$rctd_type <- unname(type_vec[colnames(obj)])
 
# Retain units with an RCTD assignment, including the explicit Unknown class.
# Prefer the marker/cluster annotation; otherwise use the RCTD label.
use <- !is.na(obj$rctd_type)
labels <- as.character(obj$cluster_naming)
fallback <- is.na(labels) | labels == ""
labels[fallback] <- as.character(obj$rctd_type)[fallback]
map_new <- c(B_cell_lineage="B cell lineage", Cardiac_muscle_cell="Ventricular cardiac muscle",
 Cardiac_neuron="Cardiac neuron", Endocardial_cell="Endothelial cell",
 Endothelial_cell="Endothelial cell", Erythroid_cell="Red blood cell", Fibroblast="Fibroblast",
 Lymphatic_endothelial="Lymphatic endothelial cell", Macrophage___Monocyte___DC="Macrophage / Monocyte / DC",
 Mast_cell="Mast cell", Pericyte="VSMC / Pericyte", T_cell_lineage="T cell lineage",
 Unknown="Unknown", VSMC="VSMC / Pericyte")
i <- labels %in% names(map_new)
labels[i] <- unname(map_new[labels[i]])
# Each supplied Valve CSV contains manually selected base barcodes in its first column.
# Only fibroblast-labelled units in these selections are reassigned.
valve_files <- list.files("annotations", pattern = "^Heart_.*_Valve\\.csv$", full.names=TRUE)
stopifnot(length(valve_files) == 4)
valve_ids <- unlist(lapply(valve_files, function(f) {
 s <- sub("_Valve\\.csv$", "", basename(f))
 ids <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)[[1]]
 paste0(s, "_", gsub("^'+|'+$", "", ids))
}))
labels[which(colnames(obj) %in% valve_ids & labels == "Fibroblast")] <- "Cardiac valve"
labels[!use] <- NA_character_
obj$Merged_cellname <- labels
obj_use <- subset(obj, subset = !is.na(Merged_cellname))
obj_use$Merged_cellname <- factor(obj_use$Merged_cellname)
saveRDS(obj_use, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming_RCTDmerged_MergedCellname_ASSIGNEDonly.rds")
 
# Per-sample Visium HD 8-um objects were filtered and XY coordinates were added from tissue_positions.parquet.
add_xy_heart <- function(obj, pos_path) {
  pos <- read_parquet(pos_path)
  coords_raw <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
                           row.names = pos$barcode)
  cells <- colnames(obj)
  cells_base <- sub("^.*_(s_008um_)", "\\1", cells)
  idx <- match(cells_base, rownames(coords_raw))
  stopifnot(all(!is.na(idx)))
  obj$x <- coords_raw$x[idx]
  obj$y <- coords_raw$y[idx]
  obj
}
 
# The additional spatial-count filter was applied to Heart_Aged_1 and Heart_Young_1.
# The saved source does not apply this 50 threshold to all four samples.
obj_hd_aged1_50 <- subset(obj_hd_aged1_lab, subset = nCount_Spatial.008um >= 50)
obj_hd_aged1_50 <- add_xy_heart(obj_hd_aged1_50,
  "Heart_Aged_1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet")
saveRDS(obj_hd_aged1_50, "Heart_Aged_1_visiumHD8um_ASSIGNEDonly_withXY.rds")
```


## 6. Adipose Visium HD cell-segmentation processing in Seurat


```r
# Input: Space Ranger v4.0.1 segmented output, segmented_outputs/filtered_feature_cell_matrix.h5
h5_seg <- lapply(adipose_outs, function(o) {
  file.path(o, "segmented_outputs/filtered_feature_cell_matrix.h5")
})
stopifnot(all(file.exists(unlist(h5_seg))))
 
odir <- "on_disk_mats_adipose"
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
 
mat_list <- lapply(names(h5_seg), function(nm) {
  m <- BPCells::open_matrix_10x_hdf5(path = h5_seg[[nm]], feature_type = "Gene Expression")
  out <- file.path(odir, nm)
  BPCells::write_matrix_dir(mat = m, dir = out, overwrite = FALSE)
  BPCells::open_matrix_dir(out)
})
names(mat_list) <- names(h5_seg)
 
obj_list <- lapply(names(mat_list), function(nm) {
  so <- CreateSeuratObject(counts = mat_list[[nm]], project = nm, assay = "RNA",
                           min.cells = 0, min.features = 0)
  so$sample <- nm
  RenameCells(so, add.cell.id = nm)
})
names(obj_list) <- names(mat_list)
 
adipose <- Reduce(function(a, b) merge(a, b), obj_list)
# Initial permissive QC removes very low-count segmented cells before sketching.
adipose <- subset(adipose, subset = nCount_RNA >= 25)
DefaultAssay(adipose) <- "RNA"
adipose <- NormalizeData(adipose, scale.factor = median(adipose$nCount_RNA), verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
adipose <- SketchData(adipose, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
 
DefaultAssay(adipose) <- "sketch"
adipose <- NormalizeData(adipose, verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
adipose <- ScaleData(adipose, verbose = FALSE)
adipose <- RunPCA(adipose, npcs = 50, verbose = FALSE)
adipose <- IntegrateLayers(adipose, method = HarmonyIntegration, orig = "pca",
                           new.reduction = "harmony", dims = 1:30)
adipose <- RunUMAP(adipose, reduction = "harmony", dims = 1:30, reduction.name = "umap")
adipose <- FindNeighbors(adipose, reduction = "harmony", dims = 1:30)
adipose <- FindClusters(adipose, resolution = 0.5)
saveRDS(adipose, "adipose_mergeobj_preProject_sketchHarmony.rds")
 
adipose <- ProjectIntegration(adipose, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
adipose <- ProjectData(adipose, assay = "RNA", full.reduction = "harmony.full",
                       sketched.assay = "sketch", sketched.reduction = "harmony.full",
                       umap.model = "umap", dims = 1:30,
                       refdata = list(cluster_full = "seurat_clusters"))
saveRDS(adipose, "adipose_mergeobj_preProject_full_sketchHarmony.rds")
```


## 7. Adipose label transfer and curated merged cell labels


```r
# Label transfer from adipose scRNA-seq reference
adipose <- readRDS("adipose_mergeobj_preProject_sketchHarmony.fixed_onDiskPath_joinCounts.rds")
ref_obj <- readRDS("adipose_merged_obj_reference.rds")
DefaultAssay(ref_obj) <- "RNA"
DefaultAssay(adipose) <- "RNA"
 
ref_obj <- FindVariableFeatures(ref_obj, verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
features_use <- intersect(VariableFeatures(ref_obj), VariableFeatures(adipose))
features_use <- head(features_use, 500)
 
anchors <- FindTransferAnchors(reference = ref_obj, query = adipose,
                               features = features_use, dims = 1:30)
pred <- TransferData(anchorset = anchors, refdata = ref_obj$Cellname, dims = 1:30)
adipose <- AddMetaData(adipose, metadata = pred)
adipose$Cellname_final <- ifelse(adipose$prediction.score.max < 0.4,
                                 "Unknown", adipose$predicted.id)
saveRDS(adipose, "adipose_mergeobj_preProject_sketchHarmony_labelTransferred.rds")
 
# Generate the SCT assay used for marker review and downstream module scoring.
# Source: Adipose_Combination_5, cell 17 (one-based).
DefaultAssay(adipose) <- "RNA"
adipose <- SCTransform(adipose, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
# Ambiguous labels were reclassified using Harmony SNN reclustering.
DefaultAssay(adipose) <- "SCT"
adipose <- FindNeighbors(adipose, reduction = "harmony", dims = 1:30,
                         k.param = 50, graph.name = "harmony_snn_k50")
adipose <- FindClusters(adipose, graph.name = "harmony_snn_k50", resolution = 0.2,
                        random.seed = 1)
adipose$clusters_k50_r020 <- adipose[["harmony_snn_k50_res.0.2"]][, 1]
 
# Assign broad cluster identities used to resolve mixed/unknown labels.
cl <- as.character(adipose$clusters_k50_r020)
new_name <- rep("Unassigned", length(cl))
new_name[cl == "0"] <- "Adipocyte (mature; lipid-signal mixed)"
new_name[cl == "1"] <- "Endothelial (capillary~general)"
new_name[cl == "2"] <- "Adipose stromal / Fibroblast"
new_name[cl == "3"] <- "Myeloid (Macrophage/DC; MHCII+; incl M2-like)"
new_name[cl == "4"] <- "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)"
new_name[cl == "5"] <- "B cell"
new_name[cl == "8"] <- "Erythroid"
new_name[cl %in% as.character(10:29)] <- "Adipocyte-like / lipid-signal dominated"
adipose$cluster_based_naming <- new_name
 
# Cluster-based curation of mixed/unknown categories.
adipose$Cellname_partial <- as.character(adipose$Cellname_final)
adipose$Cellname_partial[adipose$Cellname_partial == "Adjacent parenchymal cells"] <- "Unknown"
adipose$Cellname_partial[adipose$Cellname_partial == "Endothelial cell"] <- "Endothelial_mixed"
adipose$Cellname_partial[adipose$Cellname_partial == "Macrophage/Monocyte/DC"] <- "Myeloid_mixed"
 
# Final merged categories used downstream.
adipose$Merged_cellname <- as.character(adipose$Cellname_partial)
map_cb_to_merged <- c(
  "Adipocyte (mature; lipid-signal mixed)" = "Adipocyte",
  "Adipocyte-like / lipid-signal dominated" = "Adipocyte",
  "Adipose stromal / Fibroblast" = "Adipose stromal cell",
  "B cell" = "B cell",
  "Endothelial (capillary~general)" = "Endothelial cell",
  "Erythroid" = "Erythroid cell",
  "Unassigned" = "Unknown",
  "Myeloid (Macrophage/DC; MHCII+; incl M2-like)" = "Macrophage / Monocyte / DC",
  "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)" = "VSMC/Pericyte"
)
idx_change <- as.character(adipose$Cellname_partial) %in% c("Endothelial_mixed", "Myeloid_mixed", "Unknown")
cb <- as.character(adipose$cluster_based_naming[idx_change])
new_label <- unname(map_cb_to_merged[cb])
new_label[is.na(new_label) | new_label == ""] <- "Unknown"
adipose$Merged_cellname[idx_change] <- new_label
saveRDS(adipose, "adipose_mergeobj_sketchHarmony_SCT_labelTransferred_k50_r020_Naming_260220.rds")
```


## 8. Adipose segmentation geometry and centroid coordinates


```r
suppressPackageStartupMessages({
  library(sf)
  library(jsonlite)
  library(stringr)
  library(tibble)
  library(dplyr)
})
 
add_segmentation_centroids <- function(seu, geojson_path, id_col = "cell_id") {
  cells_sf <- st_read(geojson_path, quiet = TRUE)
  poly_geoms <- st_geometry(cells_sf)
  invalid_idx <- which(!st_is_valid(poly_geoms))
  if (length(invalid_idx) > 0) {
    cells_sf[invalid_idx, ] <- st_make_valid(cells_sf[invalid_idx, ])
    poly_geoms <- st_geometry(cells_sf)
  }
 
  calc_area_px <- function(g) as.numeric(st_area(g))
  cells_sf$seg_area_px <- vapply(poly_geoms, calc_area_px, numeric(1))
  cents <- st_centroid(poly_geoms)
  coords <- st_coordinates(cents)
  cells_sf$centroid_x <- coords[, 1]
  cells_sf$centroid_y <- coords[, 2]
 
  barcodes <- colnames(seu)
  id_num <- as.integer(str_extract(barcodes, "(?<=_)\\d+(?=-)"))
  stopifnot(sum(is.na(id_num)) == 0)
  map_df <- tibble(barcode = barcodes, !!id_col := id_num)
  cell_meta <- cells_sf |>
    st_drop_geometry() |>
    dplyr::select(!!id_col, seg_area_px, centroid_x, centroid_y)
  cell_meta_bar <- inner_join(map_df, cell_meta, by = id_col)
  cell_meta_named <- column_to_rownames(cell_meta_bar, "barcode")
  cell_meta_named <- as.data.frame(cell_meta_named)[barcodes, , drop = FALSE]
  stopifnot(identical(rownames(cell_meta_named), barcodes))
  AddMetaData(seu, metadata = cell_meta_named)
}
 
# All four annotated samples need their own segmentation geometry.
adipose_annotated <- readRDS("adipose_mergeobj_sketchHarmony_SCT_labelTransferred_k50_r020_Naming_260220.rds")
obj_list <- SplitObject(adipose_annotated, split.by = "sample")
samples <- c("Adipose_Young_1", "Adipose_Young_2", "Adipose_Aged_1", "Adipose_Aged_2")
for (sample in samples) {
  obj_list[[sample]] <- add_segmentation_centroids(
    obj_list[[sample]], file.path("GeoJason", sample, "cell_segmentations.geojson"))
  saveRDS(obj_list[[sample]], paste0(sample, "_visiumHDsementation_ASSIGNEDonly_withXY.rds"))
}
```


## 9. Aging_1 module score and spatial neighborhood definition


```r
gene_list <- c(
  "Ifit3", "Ifi44", "Isg15", "Ifit1", "Oasl2", "Ifitm3", "Bst2",
  "Ifi204", "Psmb8", "Lgals3bp", "Oasl1", "Irgm2", "Herc6", "Usp18",
  "Irgm1", "Trim30a", "Gbp3", "Mndal", "Irf7", "Xaf1", "Lgals9",
  "Phf11d", "Ifit2", "Stat1", "Gbp7", "Iigp1", "Ifi203", "Cmpk2",
  "Rsad2", "Ifit3b", "Rtp4"
)
 
# Aging_1 scores were evaluated within each spatial data representation rather
# than interpreted on a common absolute scale across tissues. Heart data were
# analyzed as fixed 8-um bins, whereas adipose data were analyzed as segmented
# cells. Tissue-specific thresholds were selected to delineate spatially
# restricted Aging_1-high foci clearly from surrounding spatial units.
 
# Heart: score and 100-um neighborhood in LV ROI (8-um bins)
library(RANN)
# Score the full annotated sample first, then subset to the manual LV ROI.
# This preserves the module-score background used in the source workflow.
heart_sample <- readRDS("Heart_Aged_1_visiumHD8um_ASSIGNEDonly_withXY.rds")
DefaultAssay(heart_sample) <- "Spatial.008um"
heart_sample <- NormalizeData(heart_sample, assay = "Spatial.008um")
gene_use <- intersect(gene_list, rownames(heart_sample))
heart_sample <- AddModuleScore(heart_sample, features = list(gene_use), name = "Aging_1_Score")
roi <- read.csv("annotations/Aged_1_LV.csv", stringsAsFactors = FALSE)
roi_col <- intersect(c("barcode", "Barcode", "BARCODE", "cellid", "CellID", "cell_id"), names(roi))[1]
stopifnot(!is.na(roi_col))
roi_barcodes <- unique(na.omit(as.character(roi[[roi_col]])))
base_barcodes <- sub("^.*_(s_008um_)", "\\1", colnames(heart_sample))
keep <- colnames(heart_sample)[base_barcodes %in% roi_barcodes]
stopifnot(length(keep) > 0)
heart_roi <- subset(heart_sample, cells = keep)
 
thr <- 0.3  # threshold used to delineate spatially restricted Aging_1-high bins
radius_um <- 100
bin_um <- 8
score_col <- "Aging_1_Score1"
seed_cells <- colnames(heart_roi)[heart_roi@meta.data[[score_col]] >= thr]
tc <- GetTissueCoordinates(heart_roi)[, c("x", "y"), drop = FALSE]
coords_all <- as.matrix(tc)
coords_seed <- coords_all[seed_cells, , drop = FALSE]
spot_px <- heart_roi@images[[1]]@scale.factors$spot
radius_px <- (radius_um / bin_um) * spot_px
nn <- RANN::nn2(data = coords_seed, query = coords_all, k = 1)
d1 <- nn$nn.dists[, 1]
names(d1) <- rownames(tc)
neighbor_cells <- setdiff(names(d1)[d1 <= radius_px], seed_cells)
heart_roi$Aging1_nbhd100um <- "other"
heart_roi$Aging1_nbhd100um[neighbor_cells] <- "neighbor"
heart_roi$Aging1_nbhd100um[seed_cells] <- "seed"
heart_roi$Aging1_nbhd100um <- factor(heart_roi$Aging1_nbhd100um,
                                     levels = c("other", "neighbor", "seed"))
saveRDS(heart_roi, "rds_out/Heart_Aged_1_LV_ROI_norm_Aging1Score.rds")
 
# Adipose: score and 100-um neighborhood by segmentation centroids (segmented cells)
adipose_a1 <- readRDS("Adipose_Aged_1_visiumHDsementation_ASSIGNEDonly_withXY.rds")
DefaultAssay(adipose_a1) <- "SCT"
adipose_a1 <- AddModuleScore(adipose_a1, features = list(gene_list), name = "Aging_1_Score")
 
thr <- 0.05  # threshold used to delineate spatially restricted Aging_1-high cells
radius_um <- 100
seed_cells <- colnames(adipose_a1)[adipose_a1@meta.data[[score_col]] >= thr]
md <- adipose_a1@meta.data
coords_all <- as.matrix(md[, c("centroid_x", "centroid_y")])
rownames(coords_all) <- rownames(md)
coords_seed <- coords_all[seed_cells, , drop = FALSE]
# Pixel-to-micron conversion estimated from the 6.5-mm capture area.
px_per_um <- mean(c(diff(range(coords_all[,1], na.rm = TRUE)),
                    diff(range(coords_all[,2], na.rm = TRUE)))) / 6500
radius_px <- radius_um * px_per_um
nn <- RANN::nn2(data = coords_seed, query = coords_all, k = 1)
d1 <- nn$nn.dists[, 1]
names(d1) <- rownames(md)
neighbor_cells <- setdiff(names(d1)[d1 <= radius_px], seed_cells)
adipose_a1$Aging1_nbhd100um <- "other"
adipose_a1$Aging1_nbhd100um[neighbor_cells] <- "neighbor"
adipose_a1$Aging1_nbhd100um[seed_cells] <- "seed"
adipose_a1$Aging1_nbhd100um <- factor(adipose_a1$Aging1_nbhd100um,
                                      levels = c("other", "neighbor", "seed"))
saveRDS(adipose_a1, "rds_out/Adipose_Aged_1_norm_Aging1Score.rds")
```


## 10. LIANA export from Seurat objects


```r
# Heart LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_heart <- function(obj, outdir) {
  DefaultAssay(obj) <- "Spatial.008um"
  obj <- JoinLayers(obj, assay = "Spatial.008um", layers = "data")
  keep_ct <- c("Ventricular cardiac muscle", "B cell lineage", "Cardiac neuron", "Endocardium",
               "Endothelial cell", "Epicardium", "Fibroblast", "Lymphatic endothelial cell",
               "Macrophage / Monocyte / DC", "Mast cell", "Platelet", "T cell lineage",
               "VSMC / Pericyte")
  obj_use <- subset(obj, subset = Merged_cellname %in% keep_ct)
  obj_use$Aging1_role <- dplyr::case_when(
    obj_use$Aging1_nbhd100um == "seed" ~ "Core",
    obj_use$Aging1_nbhd100um == "neighbor" ~ "Neighbor",
    obj_use$Aging1_nbhd100um == "other" ~ "Distal",
    TRUE ~ NA_character_
  )
  obj_use$CCI_group <- dplyr::case_when(
    obj_use$Aging1_role == "Core" ~ "Core",
    obj_use$Aging1_role == "Neighbor" ~ paste0("Neighbor_", obj_use$Merged_cellname),
    obj_use$Aging1_role == "Distal" ~ paste0("Distal_", obj_use$Merged_cellname),
    TRUE ~ NA_character_
  )
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = "Spatial.008um", layer = "data")
  stopifnot(identical(colnames(X), colnames(obj_use)))
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- data.frame(cell_id = colnames(obj_use), CCI_group = obj_use$CCI_group,
                    stringsAsFactors = FALSE)
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
 
# Adipose LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_adipose <- function(obj, outdir) {
  assay_use <- "SCT"
  DefaultAssay(obj) <- assay_use
  obj <- tryCatch(JoinLayers(obj, assay = assay_use, layers = "data"), error = function(e) obj)
  exclude_ct <- c("Erythroid cell", "Unknown")
  obj$CCI_group <- NA_character_
  obj$CCI_group[obj$Aging1_nbhd100um == "seed"] <- "Core"
  obj$CCI_group[obj$Aging1_nbhd100um == "neighbor"] <- paste0("Neighbor_", obj$Merged_cellname[obj$Aging1_nbhd100um == "neighbor"])
  obj$CCI_group[obj$Aging1_nbhd100um == "other"] <- paste0("Distal_", obj$Merged_cellname[obj$Aging1_nbhd100um == "other"])
  obj_use <- subset(obj, subset = !is.na(CCI_group) &
                    (CCI_group == "Core" | !(Merged_cellname %in% exclude_ct)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = assay_use, layer = "data")
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- obj_use@meta.data[, c("CCI_group"), drop = FALSE]
  obs$cell_id <- rownames(obs)
  obs <- obs[, c("cell_id", "CCI_group")]
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
```


## 11. Python LIANA rank aggregation and Neighbor-versus-Distal filtering


```python
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import liana as li
 
 
def load_export(outdir: str) -> sc.AnnData:
    X = sio.mmread(f"{outdir}/X.mtx").tocsr().T.tocsr()  # cells x genes
    genes = pd.read_csv(f"{outdir}/genes.txt", header=None)[0].tolist()
    cells = pd.read_csv(f"{outdir}/cells.txt", header=None)[0].tolist()
    obs = pd.read_csv(f"{outdir}/obs.csv").set_index("cell_id").loc[cells]
    adata = sc.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = cells
    adata.obs = obs
    return adata
 
 
def neighbor_distal_keep(res_df: pd.DataFrame) -> pd.DataFrame:
    df = res_df.copy()
    df["lr"] = df["ligand_complex"].astype(str) + "→" + df["receptor_complex"].astype(str)
    nb = df[(df["target"] == "Core") & (df["source"].str.startswith("Neighbor_"))].copy()
    ds = df[(df["target"] == "Core") & (df["source"].str.startswith("Distal_"))].copy()
    nb = nb[~nb["source"].isin({"Neighbor_Erythroid cell", "Neighbor_Unknown"})].copy()
    ds = ds[~ds["source"].isin({"Distal_Erythroid cell", "Distal_Unknown"})].copy()
    nb["ct"] = nb["source"].str.replace("^Neighbor_", "", regex=True)
    ds["ct"] = ds["source"].str.replace("^Distal_", "", regex=True)
    nb_mean = nb.groupby(["ct", "lr"], as_index=False)["lrscore"].mean().rename(columns={"lrscore": "mean_neighbor"})
    ds_mean = ds.groupby(["ct", "lr"], as_index=False)["lrscore"].mean().rename(columns={"lrscore": "mean_distal"})
    m = nb_mean.merge(ds_mean, on=["ct", "lr"], how="left")
    m["mean_distal"] = m["mean_distal"].fillna(0.0)
    m["delta"] = m["mean_neighbor"] - m["mean_distal"]
    return m[m["delta"] > 0].copy()
 
 
def run_liana_one(tag: str, outdir: str) -> pd.DataFrame:
    adata = load_export(outdir)
    li.mt.rank_aggregate(
        adata,
        groupby="CCI_group",
        resource_name="mouseconsensus",
        use_raw=False,
        verbose=True,
    )
    res = adata.uns["liana_res"].copy()
    keep = neighbor_distal_keep(res)
    keep["sample"] = tag
    keep.to_csv(f"{tag}_liana_neighbor_gt_distal.csv", index=False)
    return keep
 
exports = {
    "Heart_Aged1":  "Aged1_LIANA_export",
    "Heart_Aged2":  "Aged2_LIANA_export",
    "Heart_Young1": "Young1_LIANA_export",
    "Heart_Young2": "Young2_LIANA_export",
}
 
heart_liana = pd.concat([run_liana_one(tag, outdir) for tag, outdir in exports.items()], ignore_index=True)
heart_liana.to_csv("Heart_LIANA_neighbor_gt_distal_all_samples.csv", index=False)
 
exports = {
    "Adipose_Aged1":  "Adipose_Aged1_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Aged2":  "Adipose_Aged2_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Young1": "Adipose_Young1_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Young2": "Adipose_Young2_LIANA_export_SCT_CoreNeighborDistal",
}
 
adipose_liana = pd.concat([run_liana_one(tag, outdir) for tag, outdir in exports.items()], ignore_index=True)
adipose_liana.to_csv("Adipose_LIANA_neighbor_gt_distal_all_samples.csv", index=False)
```


## 12. Cross-sample LIANA summary for age-group comparison


```python
def combine_two(k1: pd.DataFrame, k2: pd.DataFrame, prefix: str) -> pd.DataFrame:
    x = k1[["ct", "lr", "mean_neighbor"]].rename(columns={"mean_neighbor": f"{prefix}_1"})
    y = k2[["ct", "lr", "mean_neighbor"]].rename(columns={"mean_neighbor": f"{prefix}_2"})
    out = x.merge(y, on=["ct", "lr"], how="outer").fillna(0.0)
    out[f"{prefix}_present_n"] = (out[[f"{prefix}_1", f"{prefix}_2"]] > 0).sum(axis=1)
    out[f"{prefix}_mean"] = out[[f"{prefix}_1", f"{prefix}_2"]].mean(axis=1)
    return out
 
# Example framework after loading per-sample Neighbor>Distal tables.
# Retain ligand-receptor pairs reproducibly observed in either two Aged samples or two Young samples.
aged = combine_two(k_aged1, k_aged2, "Aged")
young = combine_two(k_young1, k_young2, "Young")
combined = aged.merge(young, on=["ct", "lr"], how="outer").fillna(0.0)
combined["max_mean"] = combined[["Aged_mean", "Young_mean"]].max(axis=1)
combined = combined[(combined["Aged_present_n"] >= 2) | (combined["Young_present_n"] >= 2)].copy()
combined.to_csv("LIANA_age_group_summary_neighbor_gt_distal.csv", index=False)
```
