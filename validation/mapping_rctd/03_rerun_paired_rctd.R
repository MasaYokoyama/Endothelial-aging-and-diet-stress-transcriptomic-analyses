# Full paired reanalysis; memory/time can be substantial. Run only after preflight.
# Each arm starts with identical raw bins, reference, config and random seed.
args <- commandArgs(TRUE)
stopifnot(length(args)==1L,file.exists(args[1])); source(args[1])
suppressPackageStartupMessages({library(Seurat);library(BPCells);library(Matrix);library(arrow);library(spacexr)})
mapping <- readRDS(file.path(output_dir,"frozen_mappings.rds"))
heart <- readRDS(file.path(heart_dir,"heart_mergeobj_preProject_full_sketchHarmony.rds"))
heart <- JoinLayers(heart,assay="RNA",layers="counts",new="counts_joined")
counts <- GetAssayData(heart,assay="RNA",layer="counts_joined")
stopifnot(identical(rownames(counts),mapping$ids))
ref <- readRDS(file.path(heart_dir,"reference_rctd_heart.rds"))
aggregate_symbols <- function(x,map) {
 ids<-intersect(rownames(x),names(map)); syms<-unname(map[ids]); levels<-unique(syms)
 A<-sparseMatrix(i=match(syms,levels),j=seq_along(syms),x=1,dims=c(length(levels),length(ids)))
 y<-A %*% as(x[ids,,drop=FALSE],"dgCMatrix");rownames(y)<-levels; y
}
write_tsv <- function(x,n) write.table(x,file.path(output_dir,n),sep="\t",row.names=FALSE,quote=FALSE,na="")
for(s in samples) {
 old<-readRDS(file.path(heart_dir,paste0("RCTD_fullref_cutoff30_",s,".rds")))
 config<-old@config
 original_results<-old@results$results_df
 rm(old);gc()
 # Match original settings by name. The UMI threshold must be known, not guessed.
 stopifnot(!is.null(config$UMI_min))
 create_args<-config[intersect(names(config),names(formals(spacexr::create.RCTD)))]
 create_args[c("spatialRNA","reference","max_cores","keep_reference")]<-NULL
 create_args$reference<-ref;create_args$max_cores<-comparison_cores;create_args$keep_reference<-TRUE
 cells<-colnames(heart)[heart$sample==s & heart$nCount_RNA>=30]
 stopifnot(length(cells)>0)
 pos<-arrow::read_parquet(file.path(heart_dir,"Tissue_Position",s,"tissue_positions.parquet"))
 base<-sub("^.*_(s_008um_)","\\1",cells)
 ix<-match(base,pos$barcode);stopifnot(!anyNA(ix))
 xy<-data.frame(x=pos$pxl_col_in_fullres[ix],y=pos$pxl_row_in_fullres[ix],row.names=cells)
 x<-counts[,cells,drop=FALSE]
 outputs<-list(); retained<-list()
 for(arm in c("org_mm","biomart")) {
  y<-aggregate_symbols(x,mapping[[arm]])
  z<-SpatialRNA(coords=xy,counts=y,nUMI=Matrix::colSums(y))
  set.seed(comparison_seed)
  create_args$spatialRNA<-z
  fit<-do.call(spacexr::create.RCTD,create_args)
  fit<-spacexr::run.RCTD(fit,doublet_mode="doublet")
  outputs[[arm]]<-fit@results$results_df
  retained[[arm]]<-data.frame(arm=arm,input_bins=length(cells),mapped_genes=nrow(y),
    mapped_UMI=sum(y),fitted_bins=nrow(outputs[[arm]]),UMI_min=config$UMI_min)
  saveRDS(fit,file.path(output_dir,paste0(s,"_",arm,"_RCTD.rds")))
  write_tsv(cbind(barcode=rownames(outputs[[arm]]),outputs[[arm]]),paste0(s,"_",arm,"_labels.tsv"))
  rm(y,z,fit);gc()
 }
 write_tsv(do.call(rbind,retained),paste0(s,"_retention.tsv"))
 a<-outputs$org_mm;b<-outputs$biomart; common<-intersect(rownames(a),rownames(b))
 fields<-intersect(c("spot_class","first_type","second_type"),intersect(names(a),names(b)))
 per_bin<-data.frame(barcode=common)
 for(f in fields) {per_bin[[paste0(f,"_org")]]<-as.character(a[common,f]);per_bin[[paste0(f,"_biomart")]]<-as.character(b[common,f])}
 write_tsv(per_bin,paste0(s,"_paired_labels.tsv"))
 stats<-do.call(rbind,lapply(fields,function(f) {
  valid<-!is.na(a[common,f])&!is.na(b[common,f])
  data.frame(field=f,common_bins=length(common),valid_pairs=sum(valid),
    agree_n=sum(as.character(a[common,f][valid])==as.character(b[common,f][valid])),
    agree_fraction=if(sum(valid)) mean(as.character(a[common,f][valid])==as.character(b[common,f][valid])) else NA_real_)
 }))
 write_tsv(stats,paste0(s,"_agreement.tsv"))
 write_tsv(data.frame(barcode=union(setdiff(rownames(a),rownames(b)),setdiff(rownames(b),rownames(a)))),paste0(s,"_arm_exclusive_bins.tsv"))
 # Existing output versus new org.Mm arm separates rerun drift from mapping effects.
 cc<-intersect(rownames(a),rownames(original_results))
 baseline<-do.call(rbind,lapply(fields,function(f) {
  valid<-!is.na(a[cc,f])&!is.na(original_results[cc,f])
  data.frame(field=f,valid_pairs=sum(valid),agree_fraction=if(sum(valid)) mean(as.character(a[cc,f][valid])==as.character(original_results[cc,f][valid])) else NA_real_)
 }))
 write_tsv(baseline,paste0(s,"_original_vs_rerun.tsv"))
 for(f in fields) write_tsv(as.data.frame(table(org=as.character(a[common,f]),biomart=as.character(b[common,f]))),paste0(s,"_",f,"_confusion.tsv"))
 rm(x,outputs,create_args);gc()
}
writeLines(capture.output(sessionInfo()),file.path(output_dir,"RCTD_sessionInfo.txt"))
