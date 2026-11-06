# Generalized methods example; requires inputs and upstream objects documented in README.
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
