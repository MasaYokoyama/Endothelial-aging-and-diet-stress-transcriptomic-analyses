# Source evidence only. Selected original code cells, not a standalone pipeline.
# These excerpts preserve source dependencies and historical paths.
# Use run_xenium.R for the consolidated proximity workflow.


# SOURCE Muscle_IFN_Pira2_proximity_Figures.ipynb CELL 7 (one-based)
young_bst2_high_ec <- muscle_proximity$young_bst2_high |>
  dplyr::filter(
    muscle_annotation_refined ==
      "Endothelial-containing segment"
  )

aged_bst2_high_ec <- muscle_proximity$aged_bst2_high |>
  dplyr::filter(
    muscle_annotation_refined ==
      "Endothelial-containing segment"
  )

c(
  young_bst2_high_ec = nrow(young_bst2_high_ec),
  aged_bst2_high_ec = nrow(aged_bst2_high_ec)
)


# SOURCE Muscle_IFN_Pira2_proximity_Figures.ipynb CELL 8 (one-based)
mac_existing <- muscle_proximity$muscle_mac_bst2_distance_df

young_mac_bst2_ec_df <- mac_existing |>
  dplyr::filter(age_group == "Young")

aged_mac_bst2_ec_df <- mac_existing |>
  dplyr::filter(age_group == "Aged")

young_nn_bst2_ec <- RANN::nn2(
  data = young_bst2_high_ec[, c("x", "y")],
  query = young_mac_bst2_ec_df[, c("x", "y")],
  k = 1
)

aged_nn_bst2_ec <- RANN::nn2(
  data = aged_bst2_high_ec[, c("x", "y")],
  query = aged_mac_bst2_ec_df[, c("x", "y")],
  k = 1
)

young_mac_bst2_ec_df$nearest_Bst2_high_EC_um <-
  young_nn_bst2_ec$nn.dists[, 1]

aged_mac_bst2_ec_df$nearest_Bst2_high_EC_um <-
  aged_nn_bst2_ec$nn.dists[, 1]


# SOURCE Muscle_IFN_Pira2_proximity_Figures.ipynb CELL 10 (one-based)
distance_grid <- seq(0, 200, by = 5)
n_boot <- 1000
set.seed(123)

bootstrap_bst2_ec_proximity <- function(mac_df, age_group) {

  plus_dist <- mac_df$nearest_Bst2_high_EC_um[
    mac_df$Pira2_status == "Pira2-positive"
  ]

  minus_dist <- mac_df$nearest_Bst2_high_EC_um[
    mac_df$Pira2_status == "Pira2-negative"
  ]

  observed <- sapply(
    distance_grid,
    function(d) {
      100 * (
        mean(plus_dist <= d) -
        mean(minus_dist <= d)
      )
    }
  )

  boot_mat <- replicate(
    n_boot,
    {
      plus_boot <- sample(
        plus_dist,
        length(plus_dist),
        replace = TRUE
      )

      minus_boot <- sample(
        minus_dist,
        length(minus_dist),
        replace = TRUE
      )

      sapply(
        distance_grid,
        function(d) {
          100 * (
            mean(plus_boot <= d) -
            mean(minus_boot <= d)
          )
        }
      )
    }
  )

  data.frame(
    distance = distance_grid,
    observed = observed,
    ci_lower = apply(
      boot_mat,
      1,
      quantile,
      probs = 0.025
    ),
    ci_upper = apply(
      boot_mat,
      1,
      quantile,
      probs = 0.975
    ),
    age_group = age_group
  )
}


# SOURCE Muscle_IFN_Pira2_proximity_Figures.ipynb CELL 11 (one-based)
young_bst2_ec_bootstrap <- bootstrap_bst2_ec_proximity(
  young_mac_bst2_ec_df,
  "Young"
)

aged_bst2_ec_bootstrap <- bootstrap_bst2_ec_proximity(
  aged_mac_bst2_ec_df,
  "Aged"
)

muscle_bst2_ec_proximity_bootstrap_df <- dplyr::bind_rows(
  young_bst2_ec_bootstrap,
  aged_bst2_ec_bootstrap
)
