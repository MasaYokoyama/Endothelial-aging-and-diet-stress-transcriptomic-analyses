# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Adipose_IFN_Pira2_proximity.ipynb CELL 5 (one-based)
young_mac_df <- young_spatial_df %>%
  filter(
    adipose_annotation_refined ==
      "Monocyte/macrophage lineage"
  )

aged_mac_df <- aged_spatial_df %>%
  filter(
    adipose_annotation_refined ==
      "Monocyte/macrophage lineage"
  )


# SOURCE Adipose_IFN_Pira2_proximity.ipynb CELL 6 (one-based)
dim(young_mac_df)
dim(aged_mac_df)

c(
  Young_Pira2_positive_pct =
    mean(young_mac_df$Pira2 > 0) * 100,

  Aged_Pira2_positive_pct =
    mean(aged_mac_df$Pira2 > 0) * 100
)


# SOURCE Adipose_IFN_Pira2_proximity.ipynb CELL 17 (one-based)
bst2_high_cutoff <- quantile(
  c(
    young_spatial_df$Bst2,
    aged_spatial_df$Bst2
  ),
  probs = 0.90,
  na.rm = TRUE
)

bst2_high_cutoff

c(
  Young_Bst2_positive_pct =
    mean(young_spatial_df$Bst2 > 0, na.rm = TRUE) * 100,
  Aged_Bst2_positive_pct =
    mean(aged_spatial_df$Bst2 > 0, na.rm = TRUE) * 100
)


# SOURCE Adipose_IFN_Pira2_proximity.ipynb CELL 18 (one-based)
young_bst2_high <- young_spatial_df %>%
  filter(Bst2 >= bst2_high_cutoff)

aged_bst2_high <- aged_spatial_df %>%
  filter(Bst2 >= bst2_high_cutoff)
