# Small synthetic checks for alignment and correlation-ranking safeguards.
# These are software checks, not analyses of study data.
source("R/regulon_pathway_associations.R")
auc <- matrix(c(1, 2, 3, 4), nrow = 1,
              dimnames = list("Irf7", c("a", "b", "c", "d")))
scores <- data.frame(positive = c(4, 3, 2, 1), negative = c(1, 2, 3, 4),
                      row.names = c("d", "c", "b", "a"))
x <- rank_positive_pearson(auc, scores)
stopifnot(nrow(x) == 1L, x$Pathway == "positive", abs(x$R_squared - 1) < 1e-12)
y <- regulon_pathway_associations(auc, scores, "Irf7", "spearman")
stopifnot(all(abs(y$Coefficient - c(1, -1)) < 1e-12))
bad <- scores; rownames(bad)[1] <- "wrong_cell"
stopifnot(inherits(try(rank_positive_pearson(auc, bad), silent = TRUE), "try-error"))
cat("PASS: shuffled cell alignment, positive-only Pearson ranking, Spearman signs, cell mismatch rejection\n")
