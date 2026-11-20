#!/usr/bin/env Rscript
# Regression for the study's two normalization conventions and absent-motif representation.
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if (!file.exists(entry)) entry <- gsub("~+~", " ", entry, fixed=TRUE)
source(file.path(dirname(dirname(normalizePath(entry))),"R","io.R"))
run_test <- function() {
  tmp <- tempfile("primer-motif-test-")
  dir.create(tmp); on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  profile <- file.path(tmp,"profile"); dir.create(profile)
  write_tab(data.frame(chr="chr1",start=10L,end=20L,nReads=0L,MM=0,PDR=0,
                       CHALM=0,MHL=0,Entropy=0),file.path(profile,"cgi.tsv"))
  k <- do.call(paste0,expand.grid(rep(list(c("A","C","G","T","N")),6L),stringsAsFactors=FALSE))
  m <- data.frame(kmer=k,EM=integer(length(k)),BPM=integer(length(k)))
  i <- match(c("AAAAAA","CCCCCC","NNNNNN"),k)
  m$EM[i] <- c(3L,1L,4L); m$BPM[i] <- c(1L,2L,1L)
  write_tab(m,file.path(profile,"motifs.tsv"))
  manifest <- file.path(tmp,"samples.tsv"); dat <- file.path(tmp,"matrices")
  write_tab(data.frame(sample_id="001",profile="profile",sex="F"),manifest)
  read_profiles(manifest,dat)
  meta <- read_samples(dat)
  stopifnot(identical(meta$sample_id,"001"), meta$EM_valid_fraction==0.5,meta$BPM_valid_fraction==0.75)
  sel <- c("AAAAAA","CCCCCC","TTTTTT")
  for (s in c("train","validation","blind")) {
    em <- scoring_features(dat,"EM","001",sel,s,meta)
    bpm <- scoring_features(dat,"BPM","001",sel,s,meta)
    stopifnot(isTRUE(all.equal(c(em),c(.75,.25,0)*if(s=="blind") .5 else 1)),
              isTRUE(all.equal(c(bpm),c(1/3,2/3,0)*if(s=="blind") .75 else 1)))
  }
  stopifnot(is.na(load_features(dat,"EM","001","TTTTTT",impute=FALSE)[1]),
            is.na(load_features(dat,"MM","001",impute=FALSE)[1]),
            load_features(dat,"MM","001")[1]==0)
}
run_test()
uniform <- matrix(1/12,4L,12L,dimnames=list(NULL,tissues))
adjusted <- sex_adjust(uniform,c("F","M","U",NA_character_))
stopifnot(adjusted[1L,"Prostate"]==0,adjusted[2L,"Ovary"]==0,adjusted[2L,"Breast"]==0,
          isTRUE(all.equal(adjusted[3:4,,drop=FALSE],uniform[3:4,,drop=FALSE])))
cat("Motif denominators, sample keys and known/unknown sex handling passed.\n")
