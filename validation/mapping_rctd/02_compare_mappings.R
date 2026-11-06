# Same Ensembl universe, recorded one-to-many mappings, original first-row rule.
args <- commandArgs(TRUE)
stopifnot(length(args)==1L, file.exists(args[1]))
source(args[1])
suppressPackageStartupMessages({library(Seurat);library(BPCells);library(Matrix)
 library(AnnotationDbi);library(org.Mm.eg.db);library(biomaRt)})
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
raw <- readRDS(file.path(heart_dir,"heart_mergeobj_preProject_full_sketchHarmony.rds"))
ids <- rownames(raw[["RNA"]])
stopifnot(length(ids)>0, !anyDuplicated(ids))
writeLines(ids,file.path(output_dir,"ensembl_universe.txt"))
# No silent version-suffix removal; this matches the submitted input IDs.
org <- AnnotationDbi::select(org.Mm.eg.db, keys=ids, keytype="ENSEMBL",columns="SYMBOL")
names(org)[match(c("ENSEMBL","SYMBOL"),names(org))] <- c("ensembl_gene_id","symbol")
if(nzchar(biomart_table)) {
 bm <- read.delim(biomart_table,stringsAsFactors=FALSE,check.names=FALSE)
 provenance <- paste("User-supplied BioMart table:",biomart_table)
} else {
 mart <- useMart("ensembl",dataset="mmusculus_gene_ensembl",host=biomart_host)
 bm <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"),
   filters="ensembl_gene_id",values=ids,mart=mart)
 saveRDS(mart,file.path(output_dir,"queried_mart.rds"))
 provenance <- paste("NEW BioMart retrieval; original release not pinned.",
   "Host:",biomart_host,"UTC:",format(Sys.time(),tz="UTC",usetz=TRUE))
}
stopifnot(all(c("ensembl_gene_id","mgi_symbol") %in% names(bm)))
names(bm)[names(bm)=="mgi_symbol"] <- "symbol"
bm <- bm[bm$ensembl_gene_id %in% ids,c("ensembl_gene_id","symbol")]
org <- org[,c("ensembl_gene_id","symbol")]
write_tsv <- function(x,n) write.table(x,file.path(output_dir,n),sep="\t",row.names=FALSE,quote=FALSE,na="")
write_tsv(bm,"biomart_all_rows.tsv"); write_tsv(org,"org_mm_all_rows.tsv")
clean <- function(x) x[!is.na(x$symbol)&nzchar(x$symbol),,drop=FALSE]
first <- function(x) { x<-clean(x); x<-x[!duplicated(x$ensembl_gene_id),];setNames(x$symbol,x$ensembl_gene_id) }
sets <- function(x) { x<-clean(x); lapply(split(x$symbol,x$ensembl_gene_id),function(v)sort(unique(v))) }
b <- first(bm);o <- first(org); bs<-sets(bm);os<-sets(org)
res <- data.frame(ensembl_gene_id=ids,biomart_symbol=unname(b[ids]),org_mm_symbol=unname(o[ids]),stringsAsFactors=FALSE)
res$status <- ifelse(is.na(res$biomart_symbol)&is.na(res$org_mm_symbol),"unmapped_both",
 ifelse(is.na(res$biomart_symbol),"unmapped_biomart",
 ifelse(is.na(res$org_mm_symbol),"unmapped_org_mm",
 ifelse(res$biomart_symbol==res$org_mm_symbol,"concordant","discordant"))))
res$biomart_candidate_n <- vapply(ids,function(id)length(bs[[id]]),integer(1))
res$org_mm_candidate_n <- vapply(ids,function(id)length(os[[id]]),integer(1))
res$candidate_sets_equal <- vapply(ids,function(id)identical(bs[[id]],os[[id]]),logical(1))
write_tsv(res,"mapping_comparison.tsv")
write_tsv(as.data.frame(table(res$status)),"mapping_status_counts.tsv")
write_tsv(res[res$status!="concordant"|res$biomart_candidate_n>1|res$org_mm_candidate_n>1,],"mapping_review.tsv")
saveRDS(list(biomart=b,org_mm=o,ids=ids,provenance=provenance),file.path(output_dir,"frozen_mappings.rds"))
writeLines(c(provenance,"One-to-many rule: first nonempty row, matching original code.",
 "Many Ensembl IDs to one symbol: sum raw counts for BOTH RCTD arms.",
 "Seurat make.unique display names are NOT used as RCTD gene symbols.",
 "Mapping counts alone do not demonstrate unchanged RCTD results."),file.path(output_dir,"mapping_provenance.txt"))
writeLines(capture.output(sessionInfo()),file.path(output_dir,"mapping_sessionInfo.txt"))
