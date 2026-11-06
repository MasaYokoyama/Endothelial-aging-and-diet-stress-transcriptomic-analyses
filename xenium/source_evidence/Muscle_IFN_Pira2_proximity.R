# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Muscle_IFN_Pira2_proximity.ipynb CELL 5 (one-based)
young_coord <- GetTissueCoordinates(
  muscle_qc[["Young_Muscle"]],
  which = "centroids"
)

aged_coord <- GetTissueCoordinates(
  muscle_qc[["Aged_Muscle"]],
  which = "centroids"
)

muscle_meta <- FetchData(
  muscle_qc,
  vars = c(
    "age_group",
    "muscle_annotation_refined",
    "Aging_1",
    "Bst2",
    "Pira2"
  )
)

muscle_meta$cell <- rownames(muscle_meta)

young_spatial_df <- young_coord %>%
  left_join(
    muscle_meta,
    by = "cell"
  )

aged_spatial_df <- aged_coord %>%
  left_join(
    muscle_meta,
    by = "cell"
  )


# SOURCE Muscle_IFN_Pira2_proximity.ipynb CELL 6 (one-based)
mac_labels <- c(
  "Monocyte/macrophage lineage"
)

young_mac_df <- young_spatial_df %>%
  filter(
    muscle_annotation_refined %in% mac_labels
  )

aged_mac_df <- aged_spatial_df %>%
  filter(
    muscle_annotation_refined %in% mac_labels
  )

dim(young_mac_df)
dim(aged_mac_df)

c(
  Young_Pira2_positive_pct =
    mean(young_mac_df$Pira2 > 0) * 100,

  Aged_Pira2_positive_pct =
    mean(aged_mac_df$Pira2 > 0) * 100
)


# SOURCE Muscle_IFN_Pira2_proximity.ipynb CELL 18 (one-based)
summary(muscle_meta$Bst2)

quantile(
  muscle_meta$Bst2,
  probs = c(
    0,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)

table(muscle_meta$Bst2 > 0)

mean(muscle_meta$Bst2 > 0) * 100


# SOURCE Muscle_IFN_Pira2_proximity.ipynb CELL 19 (one-based)
bst2_direct <- FetchData(
  muscle_qc,
  vars = "Bst2"
)

summary(bst2_direct$Bst2)

all.equal(
  muscle_meta$Bst2,
  bst2_direct$Bst2[
    match(
      muscle_meta$cell,
      rownames(bst2_direct)
    )
  ]
)


# SOURCE Muscle_IFN_Pira2_proximity.ipynb CELL 20 (one-based)
young_bst2_high <- young_spatial_df %>%
  filter(Bst2 > 0)

aged_bst2_high <- aged_spatial_df %>%
  filter(Bst2 > 0)
