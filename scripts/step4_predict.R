#!/usr/bin/env Rscript
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
for(f in c("io","learners","ensemble")) source(file.path(root,"R",paste0(f,".R")))
source(file.path(root,"scripts","arguments.R"))
a <- read_options(c("data","model","out"),list(cohort="validation"),
    "Rscript scripts/step4_predict.R --data MATRICES --model TRAINED_MODEL --out NEW_PREDICTIONS [--cohort validation|blind]")
score_study(a$data,a$model,a$out,a$cohort)
