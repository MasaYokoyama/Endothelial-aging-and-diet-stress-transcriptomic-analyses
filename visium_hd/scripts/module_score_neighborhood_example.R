# Generalized methods example; requires inputs and upstream objects documented in README.
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
