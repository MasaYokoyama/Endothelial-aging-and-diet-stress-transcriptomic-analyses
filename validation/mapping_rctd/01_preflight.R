# Run on the computer used for the original Visium analysis.
# This does not run RCTD and does not alter any original input.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
source(args[1])
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
pkgs <- c("Seurat", "SeuratObject", "BPCells", "Matrix", "arrow",
          "spacexr", "biomaRt", "AnnotationDbi", "org.Mm.eg.db")
versions <- do.call(rbind, lapply(pkgs, function(p) {
  present <- requireNamespace(p, quietly=TRUE)
  d <- if (present) packageDescription(p) else list()
  data.frame(package=p, version=if(present) as.character(packageVersion(p)) else NA,
    remote_sha=if(is.null(d$RemoteSha)) NA else d$RemoteSha)
}))
write.table(versions, file.path(output_dir,"runtime_versions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir,"sessionInfo.txt"))
files <- c("reference_rctd_heart.rds", "heart_mergeobj_preProject_full_sketchHarmony.rds",
  paste0("RCTD_fullref_cutoff30_",samples,".rds"),
  file.path("Tissue_Position",samples,"tissue_positions.parquet"))
paths <- file.path(heart_dir,files)
write.table(data.frame(file=files,exists=file.exists(paths),bytes=file.info(paths)$size),
 file.path(output_dir,"input_availability.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if (requireNamespace("spacexr",quietly=TRUE)) {
 for (s in samples) {
  f <- file.path(heart_dir,paste0("RCTD_fullref_cutoff30_",s,".rds"))
  if (!file.exists(f)) next
  x <- readRDS(f)
  saveRDS(x@config,file.path(output_dir,paste0(s,"_original_config.rds")))
  writeLines(capture.output(str(x@config)),file.path(output_dir,paste0(s,"_original_config.txt")))
  rm(x); gc()
 }
}
if (anyNA(versions$version)) stop("Missing packages: see runtime_versions.tsv")
message("Environment inventory saved; no RCTD inference performed.")
