# Generalized methods example; requires inputs and upstream objects documented in README.
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
