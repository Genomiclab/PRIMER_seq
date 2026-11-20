#!/usr/bin/env Rscript
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
for(f in c("io","evaluation")) source(file.path(root,"R",paste0(f,".R")))
source(file.path(root,"scripts","arguments.R"))
a <- read_options(c("predictions","truth","out"),list(),
    "Rscript scripts/step5_evaluate.R --predictions PREDICTIONS.tsv --truth TRUTH.tsv --out NEW_REPORT")
evaluate(a$predictions,a$truth,a$out)
