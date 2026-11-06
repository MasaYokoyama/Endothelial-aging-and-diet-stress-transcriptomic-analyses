# Copy to config.local.R and edit paths on the original analysis computer.
heart_dir <- "/work/Visium/Heart_combined_samples"
output_dir <- "validation_results/mapping_rctd"
samples <- c("Heart_Aged_1", "Heart_Aged_2", "Heart_Young_1", "Heart_Young_2")
# Prefer the original getBM result, with ensembl_gene_id and mgi_symbol columns.
# Preserve its row order: the original pipeline retained the first mapping.
biomart_table <- ""
# Used only when biomart_table is empty. This is a NEW retrieval, not a claim
# to reconstruct the original unpinned Ensembl release. The query is frozen.
biomart_host <- "https://www.ensembl.org"
comparison_seed <- 123L
comparison_cores <- 4L
