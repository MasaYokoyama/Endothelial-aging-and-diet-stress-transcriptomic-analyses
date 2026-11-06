# Generalized methods example. See README and config/input_objects.tsv for input specifications.
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
