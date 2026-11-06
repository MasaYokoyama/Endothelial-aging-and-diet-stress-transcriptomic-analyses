# Part IV. Xenium spatial analyses

Distances and cumulative proximity curves summarize the spatial organization within
each examined tissue specimen. Each organ–age condition represents one mouse.
The optional independent-cell bootstrap describes within-specimen resampling
variation; biological replication is evaluated at the animal level.

Expression is taken from the normalized Xenium assay. Heart and adipose use the
90th percentile of Bst2 across all cells pooled over Young and Aged within each
organ, followed by selection of endothelial-containing segments. Muscle uses
Bst2-positive endothelial-containing segments. The detailed runner and settings
are in `xenium/`.

## 1. R setup and organ-specific Xenium objects


```r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(RANN)
})

# Final organ-specific Xenium objects after manual annotation refinement.
heart_qc <- readRDS("heart_qc_annotated_refined_updated.rds")
muscle_qc <- readRDS("muscle_qc_annotated_refined.rds")
adipose_qc <- readRDS("adipose_qc_annotated_refined_updated.rds")

# Final annotation columns:
#   heart_qc   : heart_annotation_refined
#   muscle_qc  : muscle_annotation_refined
#   adipose_qc : adipose_annotation_refined
#
# FOV names:
#   Heart   : Young_Heart, Aged_Heart
#   Muscle  : Young_Muscle, Aged_Muscle
#   Adipose : Young_Adipose, Aged_Adipose
```


## 2. Cell-centroid coordinates and Bst2/Pira2 expression


```r
# Xenium centroid coordinates are expressed in micrometers.
# Expression values are taken from the normalized data layer of the Xenium assay.
# The original preprocessing used NormalizeData(), not SCTransform().

make_xenium_spatial_table <- function(obj,
                                      fov,
                                      age,
                                      annotation_col,
                                      assay = "Xenium") {
  coords <- GetTissueCoordinates(obj[[fov]], which = "centroids")

  if (!"cell" %in% colnames(coords)) {
    coords$cell <- rownames(coords)
  }

  stopifnot(all(c("cell", "x", "y") %in% colnames(coords)))

  coords <- coords[, c("cell", "x", "y"), drop = FALSE]
  stopifnot(!anyDuplicated(coords$cell))
  stopifnot(all(is.finite(coords$x)), all(is.finite(coords$y)))

  md <- obj@meta.data
  stopifnot(annotation_col %in% colnames(md))
  stopifnot(all(coords$cell %in% rownames(md)))

  expr <- GetAssayData(
    obj,
    assay = assay,
    layer = "data"
  )

  stopifnot(all(c("Bst2", "Pira2") %in% rownames(expr)))
  stopifnot(all(coords$cell %in% colnames(expr)))

  expr_use <- as.data.frame(
    t(as.matrix(expr[c("Bst2", "Pira2"), coords$cell, drop = FALSE]))
  )
  expr_use$cell <- rownames(expr_use)

  out <- coords %>%
    left_join(expr_use, by = "cell")

  out$Age <- age
  out$annotation <- as.character(md[out$cell, annotation_col])
  stopifnot(all(is.finite(out$Bst2)), all(is.finite(out$Pira2)))
  if ("age_group" %in% colnames(md)) {
    stopifnot(all(as.character(md[out$cell, "age_group"]) == age))
  }

  out
}


heart_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    heart_qc,
    fov = "Young_Heart",
    age = "Young",
    annotation_col = "heart_annotation_refined"
  ),
  make_xenium_spatial_table(
    heart_qc,
    fov = "Aged_Heart",
    age = "Aged",
    annotation_col = "heart_annotation_refined"
  )
)

muscle_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    muscle_qc,
    fov = "Young_Muscle",
    age = "Young",
    annotation_col = "muscle_annotation_refined"
  ),
  make_xenium_spatial_table(
    muscle_qc,
    fov = "Aged_Muscle",
    age = "Aged",
    annotation_col = "muscle_annotation_refined"
  )
)

adipose_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    adipose_qc,
    fov = "Young_Adipose",
    age = "Young",
    annotation_col = "adipose_annotation_refined"
  ),
  make_xenium_spatial_table(
    adipose_qc,
    fov = "Aged_Adipose",
    age = "Aged",
    annotation_col = "adipose_annotation_refined"
  )
)
```


## 3. Organ-specific definition of Bst2-high endothelial-containing segments


```r
# Bst2-high cells are selected only from cells annotated as
# "Endothelial-containing segment".
#
# Heart/Adipose: calculate the 90th percentile from ALL spatial cells pooled
# across the two ages, then restrict reference cells to endothelial segments.
# Muscle: the all-cell 90th percentile was zero (5.62908% Bst2-positive),
# so the original final analysis retained all Bst2-positive endothelial segments.
# This cutoff is an empirical selection rule, not a significance threshold.

define_bst2_high_ec <- function(spatial_df, organ_name,
                                endothelial_label = "Endothelial-containing segment") {
  ec_df <- spatial_df %>%
    filter(annotation == endothelial_label)

  stopifnot(nrow(ec_df) > 0)

  positive_fraction <- mean(spatial_df$Bst2 > 0, na.rm = TRUE)
  q90 <- unname(quantile(spatial_df$Bst2, probs = 0.90, na.rm = TRUE, type = 7))

  if (organ_name == "Muscle") {
    cutoff <- 0
    rule <- "All Bst2-positive endothelial-containing segments"
    high_df <- ec_df %>%
      filter(Bst2 > 0)
  } else {
    stopifnot(q90 > 0) # Do not silently select every zero-expressing segment.
    cutoff <- q90
    rule <- "Bst2 at or above the organ-specific pooled all-cell 90th percentile"
    high_df <- ec_df %>%
      filter(Bst2 > 0, Bst2 >= cutoff)
  }

  summary_df <- data.frame(
    Organ = organ_name,
    Endothelial_n = nrow(ec_df),
    Bst2_positive_n = sum(ec_df$Bst2 > 0, na.rm = TRUE),
    All_spatial_cell_n = nrow(spatial_df),
    All_cell_Bst2_positive_fraction = positive_fraction,
    Cutoff_population = "All Young and Aged spatial cells pooled within organ",
    Bst2_high_cutoff = cutoff,
    Bst2_high_n = nrow(high_df),
    Selection_rule = rule,
    stringsAsFactors = FALSE
  )

  list(
    high_cells = high_df,
    summary = summary_df
  )
}


heart_bst2_definition <- define_bst2_high_ec(
  heart_spatial_df,
  organ_name = "Heart"
)

muscle_bst2_definition <- define_bst2_high_ec(
  muscle_spatial_df,
  organ_name = "Muscle"
)

adipose_bst2_definition <- define_bst2_high_ec(
  adipose_spatial_df,
  organ_name = "Adipose"
)

bst2_high_definition_summary <- bind_rows(
  heart_bst2_definition$summary,
  muscle_bst2_definition$summary,
  adipose_bst2_definition$summary
)

bst2_high_definition_summary

# Confirm the muscle-specific rule used in the analysis.
stopifnot(
  muscle_bst2_definition$summary$All_cell_Bst2_positive_fraction < 0.10
)
```


## 4. Nearest-neighbor distance from monocyte/macrophage-lineage cells to Bst2-high endothelial-containing segments


```r
calculate_nearest_bst2_high_ec_distance <- function(spatial_df,
                                                    bst2_high_ec_df) {
  macrophage_df <- spatial_df %>%
    filter(annotation == "Monocyte/macrophage lineage") %>%
    mutate(
      Pira2_status = ifelse(
        Pira2 > 0,
        "Pira2-positive",
        "Pira2-negative"
      )
    )

  distance_list <- lapply(c("Young", "Aged"), function(age_group) {
    query_df <- macrophage_df %>%
      filter(Age == age_group)

    reference_df <- bst2_high_ec_df %>%
      filter(Age == age_group)

    if (nrow(query_df) == 0 || nrow(reference_df) == 0) {
      stop("Missing macrophages or Bst2-high endothelial segments for ", age_group)
    }

    nn <- RANN::nn2(
      data = as.matrix(reference_df[, c("x", "y")]),
      query = as.matrix(query_df[, c("x", "y")]),
      k = 1
    )

    query_df$nearest_Bst2_high_EC_distance_um <- nn$nn.dists[, 1]
    query_df$nearest_Bst2_high_EC_cell <-
      reference_df$cell[nn$nn.idx[, 1]]

    query_df
  })

  bind_rows(distance_list) %>%
    mutate(
      Age = factor(Age, levels = c("Young", "Aged")),
      Pira2_status = factor(
        Pira2_status,
        levels = c("Pira2-negative", "Pira2-positive")
      )
    )
}
```


## 5. Descriptive proximity curves and optional cell-resampling bands


```r
summarize_cumulative_proximity <- function(distance_df,
                                           thresholds_um = seq(0, 200, 5)) {
  summary_list <- lapply(
    split(distance_df, list(distance_df$Age, distance_df$Pira2_status),
          drop = TRUE),
    function(df) {
      data.frame(
        Age = as.character(df$Age[1]),
        Pira2_status = as.character(df$Pira2_status[1]),
        Distance_um = thresholds_um,
        Cumulative_percent = vapply(
          thresholds_um,
          function(d) {
            100 * mean(df$nearest_Bst2_high_EC_distance_um <= d)
          },
          numeric(1)
        ),
        Cell_n = nrow(df),
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(summary_list) %>%
    mutate(
      Age = factor(Age, levels = c("Young", "Aged")),
      Pira2_status = factor(
        Pira2_status,
        levels = c("Pira2-negative", "Pira2-positive")
      )
    )
}


summarize_proximity_difference <- function(distance_df,
                                           thresholds_um = seq(0, 200, 5)) {
  bind_rows(lapply(c("Young", "Aged"), function(age_group) {
    dd <- distance_df[as.character(distance_df$Age) == age_group, , drop = FALSE]
    dp <- dd$nearest_Bst2_high_EC_distance_um[dd$Pira2_status == "Pira2-positive"]
    dn <- dd$nearest_Bst2_high_EC_distance_um[dd$Pira2_status == "Pira2-negative"]
    if (!length(dp) || !length(dn)) stop("Both Pira2 groups are required: ", age_group)
    data.frame(Age = age_group, Distance_um = thresholds_um,
      Difference_percent = vapply(thresholds_um,
        function(d) 100 * (mean(dp <= d) - mean(dn <= d)), numeric(1)),
      Pira2_positive_n = length(dp), Pira2_negative_n = length(dn))
  }))
}


# Legacy cell-resampling algorithm: pointwise percentile bands conditional on
# observed cells. These are NOT animal-level confidence intervals and do not
# account for spatial dependence. The caller controls whether bands are drawn.
bootstrap_proximity_difference <- function(distance_df,
                                           thresholds_um = seq(0, 200, 5),
                                           n_boot = 1000,
                                           seed = 123) {
  set.seed(seed)

  age_results <- lapply(c("Young", "Aged"), function(age_group) {
    df_age <- distance_df %>%
      filter(Age == age_group)

    d_pos <- df_age$nearest_Bst2_high_EC_distance_um[
      df_age$Pira2_status == "Pira2-positive"
    ]
    d_neg <- df_age$nearest_Bst2_high_EC_distance_um[
      df_age$Pira2_status == "Pira2-negative"
    ]

    if (length(d_pos) == 0 || length(d_neg) == 0) {
      return(NULL)
    }

    observed_difference <- vapply(
      thresholds_um,
      function(d) {
        100 * (
          mean(d_pos <= d) -
            mean(d_neg <= d)
        )
      },
      numeric(1)
    )

    boot_matrix <- replicate(
      n_boot,
      {
        d_pos_boot <- sample(d_pos, size = length(d_pos), replace = TRUE)
        d_neg_boot <- sample(d_neg, size = length(d_neg), replace = TRUE)

        vapply(
          thresholds_um,
          function(d) {
            100 * (
              mean(d_pos_boot <= d) -
                mean(d_neg_boot <= d)
            )
          },
          numeric(1)
        )
      }
    )

    data.frame(
      Age = age_group,
      Distance_um = thresholds_um,
      Difference_percent = observed_difference,
      CI_lower = apply(
        boot_matrix,
        1,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      CI_upper = apply(
        boot_matrix,
        1,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),
      Pira2_positive_n = length(d_pos),
      Pira2_negative_n = length(d_neg),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(age_results) %>%
    mutate(Age = factor(Age, levels = c("Young", "Aged")))
}


plot_proximity_difference <- function(bootstrap_df, organ_name, show_cell_bands = FALSE) {
  age_colors <- c(
    Young = "#F8766D",
    Aged = "#00BFC4"
  )

  ggplot(
    bootstrap_df,
    aes(
      x = Distance_um,
      y = Difference_percent,
      color = Age,
      fill = Age
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "grey40"
    ) +
    {if (show_cell_bands) geom_ribbon(
      aes(ymin = CI_lower, ymax = CI_upper),
      alpha = 0.20, color = NA
    )} +
    geom_line(linewidth = 1) +
    scale_color_manual(values = age_colors) +
    scale_fill_manual(values = age_colors) +
    labs(
      title = organ_name,
      subtitle = if (show_cell_bands)
        "Within-specimen cell-resampling bands; one mouse per age; not animal-level inference"
      else "Descriptive within-specimen comparison; one mouse per age",
      x = "Distance to nearest Bst2-high endothelial-containing segment (µm)",
      y = "Pira2-positive minus Pira2-negative monocyte/macrophage-lineage cells (%)"
    ) +
    theme_classic(base_size = 12)
}
```


## 6. Organ-wise Bst2–Pira2 proximity analysis


```r
run_bst2_pira2_proximity <- function(spatial_df,
                                     bst2_definition,
                                     organ_name,
                                     thresholds_um = seq(0, 200, 5),
                                     n_boot = 1000,
                                     seed = 123,
                                     include_cell_bootstrap = FALSE) {
  distance_df <- calculate_nearest_bst2_high_ec_distance(
    spatial_df = spatial_df,
    bst2_high_ec_df = bst2_definition$high_cells
  )

  cumulative_df <- summarize_cumulative_proximity(
    distance_df = distance_df,
    thresholds_um = thresholds_um
  )

  if (include_cell_bootstrap) {
    bootstrap_df <- bootstrap_proximity_difference(
      distance_df = distance_df, thresholds_um = thresholds_um,
      n_boot = n_boot, seed = seed
    )
  } else {
    # The observed difference is calculated without resampling.
    bootstrap_df <- summarize_proximity_difference(distance_df, thresholds_um)
  }

  proximity_plot <- plot_proximity_difference(
    bootstrap_df = bootstrap_df,
    organ_name = organ_name,
    show_cell_bands = include_cell_bootstrap
  )

  list(
    bst2_definition = bst2_definition$summary,
    bst2_high_ec = bst2_definition$high_cells,
    macrophage_distance = distance_df,
    cumulative_proximity = cumulative_df,
    bootstrap_difference = bootstrap_df,
    proximity_plot = proximity_plot
  )
}


heart_proximity <- run_bst2_pira2_proximity(
  spatial_df = heart_spatial_df,
  bst2_definition = heart_bst2_definition,
  organ_name = "Heart"
)

muscle_proximity <- run_bst2_pira2_proximity(
  spatial_df = muscle_spatial_df,
  bst2_definition = muscle_bst2_definition,
  organ_name = "Muscle"
)

adipose_proximity <- run_bst2_pira2_proximity(
  spatial_df = adipose_spatial_df,
  bst2_definition = adipose_bst2_definition,
  organ_name = "Adipose"
)

heart_proximity$proximity_plot
muscle_proximity$proximity_plot
adipose_proximity$proximity_plot
```


## 7. Save proximity-analysis results


```r
saveRDS(
  heart_proximity,
  file = "Heart_Bst2_Pira2_proximity_analysis.rds"
)

saveRDS(
  muscle_proximity,
  file = "Muscle_Bst2_Pira2_proximity_analysis.rds"
)

saveRDS(
  adipose_proximity,
  file = "Adipose_Bst2_Pira2_proximity_analysis.rds"
)

write.csv(
  bst2_high_definition_summary,
  file = "Bst2_high_definition_summary.csv",
  row.names = FALSE
)
```
