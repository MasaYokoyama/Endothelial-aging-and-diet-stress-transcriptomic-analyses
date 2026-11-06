# Rscript scripts/run_neighborhood_export.R VISIUM_DATA_ROOT NEW_OUTPUT_DIRECTORY config/samples.tsv
# Uses saved, already-normalized module-score objects. Does not recalculate scores within a smaller ROI.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Supply VISIUM_DATA_ROOT, a NEW output directory, and samples.tsv")
root <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
if(dir.exists(out)) stop("Output directory already exists; choose a fresh directory")
manifest <- read.delim(args[3], stringsAsFactors=FALSE, check.names=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(RANN); library(dplyr)})
script_arg <- grep("^--file=", commandArgs(), value=TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir,"export_liana.R"))
stopifnot(nrow(manifest)==8, !anyDuplicated(manifest$sample))
input_files <- file.path(root,manifest$scored_rds_relative)
if(any(!file.exists(input_files))) stop(paste("Missing inputs:",paste(input_files[!file.exists(input_files)],collapse="\n")))
dir.create(out,recursive=TRUE)
qc <- list()
for(i in seq_len(nrow(manifest))) {
 p <- manifest[i,]; message(p$sample)
 obj <- readRDS(input_files[i])
 md <- obj[[]]
 stopifnot(all(c("Aging_1_Score1","Merged_cellname") %in% names(md)),
           identical(rownames(md),colnames(obj)), all(is.finite(md$Aging_1_Score1)))
 if(p$tissue=="Heart") {
  tc <- GetTissueCoordinates(obj)
  stopifnot(all(colnames(obj) %in% rownames(tc)))
  coords <- as.matrix(tc[colnames(obj),c("x","y"),drop=FALSE])
  # Source's image-coordinate calibration: 8-um bin width in image pixels.
  spot_px <- obj@images[[1]]@scale.factors$spot
  radius_px <- p$radius_um / 8 * spot_px
  calibration <- "image spot scale / 8 um"
 } else {
  coords <- as.matrix(md[,c("centroid_x","centroid_y"),drop=FALSE])
  # Original estimate, retained transparently. It assumes the observed centroid
  # span represents the 6.5-mm capture area in the original coordinate conversion.
  px_per_um <- mean(apply(coords,2,function(x) diff(range(x))))/6500
  radius_px <- p$radius_um * px_per_um
  calibration <- "mean observed centroid span / 6500 um (original approximation)"
 }
 stopifnot(all(is.finite(coords)),is.finite(radius_px),radius_px>0)
 seeds <- colnames(obj)[md$Aging_1_Score1 >= p$core_threshold]
 if(!length(seeds)) stop(paste(p$sample,"has no core units at the declared threshold"))
 nn <- RANN::nn2(data=coords[seeds,,drop=FALSE],query=coords,k=1)
 labels <- rep("other",ncol(obj)); names(labels)<-colnames(obj)
 labels[nn$nn.dists[,1] <= radius_px] <- "neighbor"
 labels[seeds] <- "seed"
 obj$Aging1_nbhd100um <- factor(labels,levels=c("other","neighbor","seed"))
 expdir <- file.path(out,p$liana_tag)
 if(p$tissue=="Heart") export_liana_heart(obj,expdir) else export_liana_adipose(obj,expdir)
 qc[[i]] <- data.frame(sample=p$sample,tissue=p$tissue,n_total=ncol(obj),
  n_core=sum(labels=="seed"),n_neighbor=sum(labels=="neighbor"),n_distal=sum(labels=="other"),
  core_threshold=p$core_threshold,radius_um=p$radius_um,radius_px=radius_px,calibration=calibration)
 write.csv(qc[[i]],file.path(expdir,"neighborhood_qc.csv"),row.names=FALSE)
 rm(obj);gc()
}
write.csv(do.call(rbind,qc),file.path(out,"neighborhood_qc_all.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
