# Synthetic helper checks; parse scripts but do not run their top-level I/O.
# No database queries, reference data, or RCTD fitting are performed.
suppressPackageStartupMessages(library(Matrix))
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(normalizePath(script)), ".."))
files <- list.files(root, pattern = "\\.R$", full.names = TRUE)
for (file in files) {parse(file = file); cat("PARSE_OK", basename(file), "\n")}
extract_functions <- function(path, wanted) {
  expressions <- parse(file = path)
  for (e in expressions) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        is.symbol(e[[2L]]) && as.character(e[[2L]]) %in% wanted &&
        is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
      eval(e, envir = .GlobalEnv)
    }
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)))
}
extract_functions(file.path(root, "02_compare_mappings.R"), c("clean", "first", "sets"))
extract_functions(file.path(root, "03_rerun_paired_rctd.R"), "aggregate_symbols")
mapping_rows <- data.frame(ensembl_gene_id = c("e1","e1","e2","e3","e4","e5","e6"),
                          symbol = c("S1","ALT","S1","S3",NA,"","S6"),
                          stringsAsFactors = FALSE)
m <- first(mapping_rows)
stopifnot(identical(names(m),c("e1","e2","e3","e6")),
          identical(unname(m),c("S1","S1","S3","S6")))
s <- sets(mapping_rows)
stopifnot(identical(s$e1,c("ALT","S1")),length(s$e4)==0L,length(s$e5)==0L)
# Two Ensembl IDs collapse to S1; unannotated e4 counts must be discarded.
x <- Matrix(matrix(c(1,2,3,4,10,20,30,40), nrow = 4,
                   dimnames = list(c("e1","e2","e3","e4"),c("bin1","bin2"))), sparse = TRUE)
y <- aggregate_symbols(x,m)
expected <- matrix(c(3,3,30,30),nrow=2,dimnames=list(c("S1","S3"),c("bin1","bin2")))
stopifnot(identical(dimnames(y),dimnames(expected)),
          isTRUE(all.equal(as.matrix(y),expected)),
          identical(as.numeric(colSums(y)),c(6,60)))
# Row order differences do not change symbol-aligned totals.
y_reordered <- aggregate_symbols(x[c("e3","e2","e1","e4"),,drop=FALSE],m)
stopifnot(isTRUE(all.equal(as.matrix(y_reordered[rownames(y),,drop=FALSE]),as.matrix(y))))
# Mapping order is deliberately preserved, not alphabetically reselected.
m_reversed <- first(mapping_rows[c(2,1,3:7),])
stopifnot(m_reversed[["e1"]]=="ALT")
cat("PASS: first nonempty mapping, one-to-many candidates, unmapped IDs, raw count summation, barcode preservation, row-order invariance.\n")
# Diagnose whether direct comparisons of factor labels need normalization.
a <- factor(c("singlet","doublet"), levels=c("singlet","doublet"))
b <- factor(c("singlet","reject"), levels=c("singlet","reject"))
err <- try(a==b,silent=TRUE)
cat("FACTOR_DIRECT_COMPARE_ERROR",inherits(err,"try-error"),"\n")
stopifnot(identical(as.character(a)==as.character(b),c(TRUE,FALSE)))
