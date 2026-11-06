# Meaningful synthetic checks only; no study data or external files are loaded.
# Run from repository root: Rscript xenium/tests/check_proximity.R
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(RANN)})
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
source(file.path(dirname(normalizePath(script)), "..", "proximity_functions.R"))
cutoff_df <- data.frame(cell = paste0("c", 1:10), Bst2 = 1:10,
                       annotation = ifelse(1:10 %in% c(2,8), "Endothelial-containing segment", "Other"),
                       Age = rep(c("Young", "Aged"), each = 5))
def <- define_bst2_high_ec(cutoff_df, "Heart")
stopifnot(abs(def$summary$Bst2_high_cutoff - 9.1) < 1e-10, nrow(def$high_cells) == 0)
# A within-EC percentile would incorrectly select c8: this check detects that regression.
df <- data.frame(cell = c("YE", "AE", "YP", "YN", "AP", "AN"),
  Age = c("Young", "Aged", "Young", "Young", "Aged", "Aged"),
  x = c(0,100,3,6,103,106), y = c(0,100,4,8,104,108),
  Bst2 = c(1,1,0,0,0,0), Pira2 = c(0,0,1,0,1,0),
  annotation = c(rep("Endothelial-containing segment",2), rep("Monocyte/macrophage lineage",4)))
muscle <- define_bst2_high_ec(df, "Muscle")
d <- calculate_nearest_bst2_high_ec_distance(df, muscle$high_cells)
stopifnot(identical(as.character(d$cell), c("YP","YN","AP","AN")))
stopifnot(identical(as.numeric(d$nearest_Bst2_high_EC_distance_um), c(5,10,5,10)))
stopifnot(identical(as.character(d$nearest_Bst2_high_EC_cell), c("YE","YE","AE","AE")))
s <- summarize_proximity_difference(d, c(0,5,10))
stopifnot(identical(s$Difference_percent, rep(c(0,100,0),2)))
r <- run_bst2_pira2_proximity(df, muscle, "Muscle", thresholds_um = c(0,5,10))
stopifnot(!"CI_lower" %in% names(r$bootstrap_difference))
# Legacy bootstrap reproducibility with at least two cells in each group.
d2 <- bind_rows(d, d)
b1 <- bootstrap_proximity_difference(d2, thresholds_um = c(0,5,10), n_boot = 20, seed = 123)
b2 <- bootstrap_proximity_difference(d2, thresholds_um = c(0,5,10), n_boot = 20, seed = 123)
stopifnot(identical(b1,b2), all(b1$CI_lower <= b1$CI_upper))
cat("PASS: all-cell cutoff, organ-specific rule, age-restricted neighbors, observed differences, default no CI, bootstrap reproducibility.\n")
