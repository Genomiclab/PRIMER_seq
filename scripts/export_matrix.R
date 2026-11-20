#!/usr/bin/env Rscript
# Write the model-input matrix as gzipped TSV with explicit sample and feature keys.
a <- commandArgs(TRUE)
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if (!file.exists(entry)) entry <- gsub("~+~", " ", entry, fixed=TRUE)
source(file.path(dirname(dirname(normalizePath(entry))),"R","io.R"))
if(length(a)!=4L || !a[2L] %in% c("train","validation","blind") || !a[3L] %in% channels ||
   !endsWith(a[4L],".tsv.gz")) stop("Usage: Rscript scripts/export_matrix.R MATRICES train|validation|blind CHANNEL NEW.tsv.gz")
export <- function() {
  if(file.exists(a[4L]) || !dir.exists(dirname(a[4L]))) stop("Provide a new .tsv.gz in an existing directory")
  meta <- read_samples(a[1L])
  x <- scoring_features(a[1L],a[3L],meta$sample_id,NULL,a[2L],meta)
  df <- data.table::as.data.table(x,keep.rownames="sample_id")
  tmp <- tempfile(".matrix-",tmpdir=dirname(a[4L]),fileext=".tsv.gz")
  on.exit(unlink(tmp),add=TRUE)
  data.table::fwrite(df,tmp,sep="\t",compress="gzip",na="NA")
  if(file.exists(a[4L]) || !file.rename(tmp,a[4L])) stop("Could not publish matrix")
}
export()
