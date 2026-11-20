#!/usr/bin/env Rscript
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
for(f in c("io","learners","ensemble")) source(file.path(root,"R",paste0(f,".R")))
source(file.path(root,"scripts","arguments.R"))
a <- read_options(c("train","out"),list(markers=file.path(root,"marker_file","markers.tsv"),threads="1",preset="study",select_motifs="false"),
    "Rscript scripts/step3_train_ensemble.R --train MATRICES --out NEW_MODEL [--markers MARKERS.tsv] [--preset demo|study] [--threads N] [--select_motifs true|false]")
if(!a$select_motifs %in% c("true","false")) stop("select_motifs must be true or false")
fit_study(a$train,a$markers,a$out,algorithms,as.integer(a$threads),a$select_motifs=="false",a$preset)
