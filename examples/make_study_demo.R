#!/usr/bin/env Rscript
# Extend the supplied synthetic methylation tables to the complete marker panel.
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R","io.R"))
source(file.path(root,"scripts","arguments.R"))
a <- read_options("out",list(base=file.path(root,"examples","three_cancers"),
  markers=file.path(root,"marker_file","markers.tsv")),
  "Rscript examples/make_study_demo.R --out NEW_TABLE_DIRECTORY [--base SYNTHETIC_TABLES] [--markers MARKERS.tsv]")
markers <- read_tab(a$markers)
in_new_directory(a$out,function(out) {
  cohorts <- c("train","validation","blind")
  for(cohort in cohorts) {
    dir.create(file.path(out,cohort))
    for(name in c("samples.tsv","truth.tsv")) {
      ok <- file.copy(file.path(a$base,cohort,name),file.path(out,cohort,name))
      if(!ok) stop("Cannot copy synthetic metadata")
    }
  }
  for(ch in channels) {
    train <- read_tab(file.path(a$base,"train",paste0(ch,".tsv.gz")))
    original <- setdiff(names(train),"sample_id")
    selected <- unique(markers$feature[markers$channel==ch])
    if(ch %in% c("EM","BPM")) {
      if(!all(selected %in% original)) stop("Synthetic motif table is missing required hexamers")
      for(cohort in cohorts) file.copy(file.path(a$base,cohort,paste0(ch,".tsv.gz")),
                                      file.path(out,cohort,paste0(ch,".tsv.gz")))
      next
    }
    extra <- setdiff(selected,original)
    set.seed(20261120L+match(ch,channels))
    left <- sample(original,length(extra),replace=TRUE)
    right <- sample(original,length(extra),replace=TRUE)
    weight <- runif(length(extra),.2,.8)
    for(cohort in cohorts) {
      d <- read_tab(file.path(a$base,cohort,paste0(ch,".tsv.gz")))
      x <- as.matrix(d[,original,drop=FALSE])
      set.seed(20261120L+100L*match(ch,channels)+match(cohort,cohorts))
      values <- sweep(x[,left,drop=FALSE],2,weight,"*") +
                sweep(x[,right,drop=FALSE],2,1-weight,"*") +
                matrix(rnorm(nrow(x)*length(extra),sd=.02),nrow(x),length(extra))
      values[values<0] <- 0
      values[values>1] <- 1
      colnames(values) <- extra
      x <- cbind(x,values)[,selected,drop=FALSE]
      con <- gzfile(file.path(out,cohort,paste0(ch,".tsv.gz")),"wt")
      tryCatch(write.table(data.frame(sample_id=d$sample_id,x,check.names=FALSE),con,
                            sep="\t",quote=FALSE,row.names=FALSE),finally=close(con))
    }
  }
  writeLines("Synthetic values generated from the supplied demonstration tables; not participant measurements.",
             file.path(out,"README.txt"))
})
