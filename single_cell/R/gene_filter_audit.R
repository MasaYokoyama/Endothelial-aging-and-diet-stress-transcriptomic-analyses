# Historical name-pattern filter, preserved for reproducibility.
# A pattern match is not a verified gene-biotype or noise annotation.
gene_filter_pattern <- paste(
  "mt-", "MT-", "MMU", "SNO", "Gm", "Rn6", "RNA", "7SK",
  "SNORD", "SNORA", "SCARNA", "B3g", "Vmn", "Mir", "Rik",
  "Snora", "Snord", "LOC", "Rn4", "OTTMUSG", "Scarna",
  "RNU", "RMR", "Rpl", "Rps", "Rnu", "Rmr", "RPl", "RPS",
  "AW", "BC", "Malat1", "H2-", "B2m", "Actb",
  "AA", "AB", "AC", "AF", "AI", "AL", "AU", "AV", "AY", "Hist", sep = "|")

audit_gene_filter <- function(gene_names, output_file) {
  terms <- strsplit(gene_filter_pattern, "|", fixed = TRUE)[[1]]
  matches <- vapply(gene_names, function(g) {
    paste(terms[vapply(terms, grepl, logical(1), x = g)], collapse = ";")
  }, character(1))
  out <- data.frame(gene = gene_names,
                     excluded = grepl(gene_filter_pattern, gene_names),
                     matched_terms = matches)
  write.csv(out, output_file, row.names = FALSE)
  invisible(out)
}
# Apply only at the documented marker/correlation steps. Do not silently apply
# this broad filter to SCENIC, IFN signature scoring, or all downstream assays.
# In particular H2-, B2m, and B3g patterns can exclude functional immune genes;
# unanchored terms can match within symbols. Inspect the exported table.
