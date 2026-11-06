# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Liver_IFN_Pira2_proximity.ipynb CELL 5 (one-based)
young_mac_df <- young_spatial_df %>%
  filter(
    liver_annotation_refined ==
      "Monocyte/macrophage lineage"
  )

aged_mac_df <- aged_spatial_df %>%
  filter(
    liver_annotation_refined ==
      "Monocyte/macrophage lineage"
  )


# SOURCE Liver_IFN_Pira2_proximity.ipynb CELL 19 (one-based)
bst2_high_cutoff <- quantile(
  c(
    young_spatial_df$Bst2,
    aged_spatial_df$Bst2
  ),
  probs = 0.90,
  na.rm = TRUE
)

bst2_high_cutoff
