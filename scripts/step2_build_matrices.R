#!/usr/bin/env Rscript
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"scripts","arguments.R"))
source(file.path(root,"R","io.R"))
a <- read_options("out",list(profiles="",tables=""),
    "Rscript scripts/step2_build_matrices.R --profiles SAMPLES.tsv | --tables TABLE_DIRECTORY --out NEW_MATRICES")
if(sum(nzchar(c(a$profiles,a$tables)))!=1L) stop("Supply exactly one of --profiles or --tables")
if(nzchar(a$profiles)) read_profiles(a$profiles,a$out) else import_tables(a$tables,a$out)
