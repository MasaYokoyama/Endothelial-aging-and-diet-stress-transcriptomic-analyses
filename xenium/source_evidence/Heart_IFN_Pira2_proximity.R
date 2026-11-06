# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 5 (one-based)
young_coord <- GetTissueCoordinates(
  heart_qc[["Young_Heart"]],
  which = "centroids"
)

aged_coord <- GetTissueCoordinates(
  heart_qc[["Aged_Heart"]],
  which = "centroids"
)


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 7 (one-based)
heart_meta <- FetchData(
  heart_qc,
  vars = c(
    "age_group",
    "heart_annotation_refined",
    "Aging_1",
    "Bst2",
    "Pira2",
    "Bst2_EC_display",
    "Pira2_Mac_display"
  )
)

heart_meta$cell <- rownames(heart_meta)


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 8 (one-based)
young_spatial_df <- young_coord %>%
  left_join(
    heart_meta,
    by = "cell"
  )

aged_spatial_df <- aged_coord %>%
  left_join(
    heart_meta,
    by = "cell"
  )


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 11 (one-based)
young_mac_df <- young_spatial_df %>%
  filter(
    heart_annotation_refined == "Monocyte/macrophage lineage"
  )

aged_mac_df <- aged_spatial_df %>%
  filter(
    heart_annotation_refined == "Monocyte/macrophage lineage"
  )

dim(young_mac_df)
dim(aged_mac_df)


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 62 (one-based)
bst2_high_cutoff <- quantile(
  c(
    young_spatial_df$Bst2,
    aged_spatial_df$Bst2
  ),
  probs = 0.90,
  na.rm = TRUE
)

bst2_high_cutoff


# SOURCE Heart_IFN_Pira2_proximity.ipynb CELL 63 (one-based)
young_bst2_high <- young_spatial_df %>%
  filter(
    Bst2 >= bst2_high_cutoff
  )

aged_bst2_high <- aged_spatial_df %>%
  filter(
    Bst2 >= bst2_high_cutoff
  )

dim(young_bst2_high)
dim(aged_bst2_high)
