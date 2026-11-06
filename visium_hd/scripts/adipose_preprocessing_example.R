# Generalized methods example. See README and config/input_objects.tsv for input specifications.
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
