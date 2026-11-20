#!/usr/bin/env Rscript
mode <- commandArgs(TRUE)
if(length(mode)!=1L || !mode %in% c("bam","train","all"))
  stop("Usage: Rscript scripts/check_environment.R bam|train|all")
bam <- c("data.table","Rsamtools","Biostrings","GenomeInfoDb","GenomicAlignments","BSgenome.Hsapiens.UCSC.hg19")
train <- c("data.table","caret","glmnet","gbm","e1071","xgboost","catboost","doParallel","foreach")
pkgs <- switch(mode,bam=bam,train=train,all=unique(c(bam,train)))
available <- vapply(pkgs,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))
print(data.frame(package=pkgs,version=vapply(pkgs,function(p)
  if(available[[p]]) as.character(packageVersion(p)) else "MISSING","")),row.names=FALSE)
if(any(!available)) stop("Install the missing R packages shown above")
if("xgboost" %in% pkgs && (packageVersion("xgboost") < "1.7.0" || packageVersion("xgboost") >= "2.0.0"))
  stop("This implementation uses the XGBoost 1.7 R API; see docs/ENVIRONMENT.md")
if(mode %in% c("bam","all")) {
  bins <- Sys.which(c("mhaptools","mhaptk","samtools","tabix","python3"))
  print(bins)
  if(any(!nzchar(bins))) stop("Put the missing executables on PATH")
}
