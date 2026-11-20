#!/usr/bin/env Rscript
# Raw counts and both frequency conventions must produce the same model inputs.
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R","io.R"))
source(file.path(root,"R","processed.R"))
expect_error <- function(expr,pattern) {
  message <- tryCatch({force(expr);NA_character_},error=function(e)conditionMessage(e))
  if(is.na(message) || !grepl(pattern,message,fixed=TRUE)) stop("Expected error: ",pattern)
}
run_test <- function() {
  tmp <- tempfile("primer-processed-test-"); dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  k <- do.call(paste0,expand.grid(rep(list(c("A","C","G","T")),6L),stringsAsFactors=FALSE))
  x <- matrix(0,2L,4097L,dimnames=list(c("001","002"),c(k,"NNNNNN")))
  x[,c("AAAAAA","CCCCCC","NNNNNN")] <- rbind(c(3,1,4),c(1,2,1))
  x[,"TTTTTT"] <- NA_real_
  expected <- normalize_motif_matrix(x,"counts")
  stopifnot(identical(unname(expected$fraction),c(.5,.75)),is.na(expected$matrix[1,"TTTTTT"]))
  canonical <- expected$matrix; canonical[is.na(canonical)] <- 0
  all_freq <- x/rowSums(x,na.rm=TRUE)
  for(v in list(normalize_motif_matrix(canonical,"canonical-frequencies",expected$fraction),
                normalize_motif_matrix(all_freq,"all-frequencies"),
                normalize_motif_matrix(all_freq[,k],"all-frequencies"),
                normalize_motif_matrix(x[,k],"counts",expected$fraction))) {
    stopifnot(isTRUE(all.equal(v$matrix,expected$matrix)),isTRUE(all.equal(v$fraction,expected$fraction)))
  }
  expect_error(normalize_motif_matrix(canonical,"canonical-frequencies"),"requires the original")
  expect_error(normalize_motif_matrix(x,"counts",c(.9,.9)),"disagrees")
  expect_error(normalize_motif_matrix(x/2,"counts"),"integers")
  expect_error(normalize_motif_matrix(canonical/2,"canonical-frequencies",c(1,1)),"sum to one")

  # Preserve leading-zero identifiers, matrix orientation and manifest row order.
  ids <- c("002","001")
  for(ext in c("tsv","csv","rds","RData")) {
    file <- file.path(tmp,paste0("matrix.",ext))
    if(ext=="rds") saveRDS(x,file)
    else if(ext=="RData") save(x,file=file)
    else write.table(data.frame(sample_id=rownames(x),x,check.names=FALSE),file,
                     sep=if(ext=="csv") "," else "\t",quote=ext=="csv",row.names=FALSE)
    stopifnot(isTRUE(all.equal(read_motif_matrix(file,ids),x[ids,])))
  }
  transposed <- file.path(tmp,"transposed.tsv")
  write.table(data.frame(motif=colnames(x),t(x),check.names=FALSE),transposed,
              sep="\t",quote=FALSE,row.names=FALSE)
  stopifnot(isTRUE(all.equal(read_motif_matrix(transposed,ids),x[ids,])))
  bad <- file.path(tmp,"bad.rds"); saveRDS(x[,-1],bad)
  expect_error(read_motif_matrix(bad,ids),"4096")
  save(x,canonical,file=file.path(tmp,"multiple.RData"))
  expect_error(read_motif_matrix(file.path(tmp,"multiple.RData"),ids),"exactly one")

  # A processed CGI+matrix import must agree with the BAM-profile assembly.
  profile <- file.path(tmp,"profile"); dir.create(profile)
  write_tab(data.frame(chr=c("chr1","chr1"),start=c(10,30),end=c(20,40),nReads=c(4,0),
                       MM=c(.5,0),PDR=c(.25,0),CHALM=c(.5,0),MHL=c(.1,0),Entropy=c(.2,0)),file.path(profile,"cgi.tsv"))
  all_k <- do.call(paste0,expand.grid(rep(list(c("A","C","G","T","N")),6L),stringsAsFactors=FALSE))
  motifs <- data.frame(kmer=all_k,EM=0,BPM=0)
  i <- match(c("AAAAAA","CCCCCC","NNNNNN"),all_k)
  motifs$EM[i] <- c(3,1,4); motifs$BPM[i] <- c(1,2,1)
  write_tab(motifs,file.path(profile,"motifs.tsv"))
  write_tab(data.frame(sample_id="001",sex="U",profile="profile"),file.path(tmp,"profile.tsv"))
  read_profiles(file.path(tmp,"profile.tsv"),file.path(tmp,"from_profile"))
  write_tab(data.frame(sample_id="001",sex="U",cgi="profile/cgi.tsv"),file.path(tmp,"cgi.tsv"))
  for(ch in c("EM","BPM")) {
    mat <- matrix(motifs[[ch]],1L,dimnames=list("001",all_k))
    saveRDS(mat,file.path(tmp,paste0(ch,".rds")))
  }
  import_processed(file.path(tmp,"cgi.tsv"),file.path(tmp,"EM.rds"),file.path(tmp,"BPM.rds"),file.path(tmp,"from_processed"))
  for(ch in channels) stopifnot(identical(readRDS(file.path(tmp,"from_profile",paste0(ch,".rds"))),
                                        readRDS(file.path(tmp,"from_processed",paste0(ch,".rds")))))
  stopifnot(identical(read_samples(file.path(tmp,"from_profile")),read_samples(file.path(tmp,"from_processed"))))
  meta <- read_samples(file.path(tmp,"from_processed"))
  blind <- scoring_features(file.path(tmp,"from_processed"),"EM","001",c("AAAAAA","CCCCCC"),"blind",meta)
  stopifnot(isTRUE(all.equal(as.numeric(blind),c(3/8,1/8))))
}
run_test()
cat("Processed matrices, frequency conventions, identifiers and profile equivalence passed.\n")
