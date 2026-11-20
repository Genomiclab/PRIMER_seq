#!/usr/bin/env Rscript
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"scripts","arguments.R"))
source(file.path(root,"R","io.R"))
source(file.path(root,"R","processed.R"))
a <- read_options("out",list(profiles="",tables="",cgi="",em="",bpm="",motif_format="counts",layout="auto"),
    "Rscript scripts/step2_build_matrices.R --profiles SAMPLES.tsv | --tables TABLE_DIRECTORY | --cgi CGI_MANIFEST --em MATRIX --bpm MATRIX --out NEW_MATRICES [--motif_format counts|canonical-frequencies|all-frequencies] [--layout auto|samples-by-features|features-by-samples]")
if(sum(nzchar(c(a$profiles,a$tables,a$cgi)))!=1L) stop("Supply exactly one of --profiles, --tables or --cgi")
if(nzchar(a$cgi)) {
  if(!nzchar(a$em) || !nzchar(a$bpm)) stop("CGI input requires --em and --bpm")
  import_processed(a$cgi,a$em,a$bpm,a$out,a$motif_format,a$layout)
} else {
  if(nzchar(a$em) || nzchar(a$bpm)) stop("Separate motif matrices require --cgi")
  if(nzchar(a$profiles)) read_profiles(a$profiles,a$out) else import_tables(a$tables,a$out)
}
