channels <- c("MM","PDR","CHALM","MHL","Entropy","EM","BPM")
tissues <- c("Lung","Colorectum","Stomach","Liver","Thyroid","Breast","Lymphoma","Esophagus","Pancreas","Bladder","Prostate","Ovary")
algorithms <- c("LR","GBM","SVM_linear","XGBoost","Catboost")
archive_channel <- function(channel) switch(channel, EM="mer6_raw_4bp", BPM="BPM_raw_4bp", channel)
read_tab <- function(path) {
  hdr <- names(read.delim(path,nrows=0L,check.names=FALSE))
  cl <- rep(NA_character_,length(hdr)); cl[hdr %in% c("sample_id","participant_id","profile","sex","group","stage","feature","kmer")] <- "character"
  read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="",colClasses=cl)
}
write_tab <- function(x,path) write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
study_sex <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x==""] <- "U"
  if(any(!x %in% c("F","M","U"))) stop("Use F, M or U (unknown) for sex")
  x
}
in_new_directory <- function(dest, fn) {
  if (file.exists(dest)) stop("Use a new output directory: ",dest)
  parent <- dirname(dest)
  if (!dir.exists(parent) && !dir.create(parent,recursive=TRUE)) stop("Cannot create ",parent)
  tmp <- tempfile(".primer-",tmpdir=parent)
  if (!dir.create(tmp)) stop("Cannot create working directory")
  on.exit(unlink(tmp,recursive=TRUE), add=TRUE)
  fn(tmp)
  if (file.exists(dest) || !file.rename(tmp,dest)) stop("Cannot publish completed output to ",dest)
}
read_profiles <- function(manifest, outfile) {
  meta <- read_tab(manifest)
  needed <- c("sample_id","profile","sex")
  if (!all(needed %in% names(meta)) || !nrow(meta)) stop("Manifest requires sample_id, profile, sex")
  ids <- as.character(meta$sample_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("sample_id must be unique and nonempty")
  meta$sex <- study_sex(meta$sex)
  if ("group" %in% names(meta)) {
    meta$group[which(meta$group=="")] <- NA_character_
    if(any(!is.na(meta$group) & !meta$group %in% c("Normal",tissues))) stop("Unknown diagnosis in group")
  }
  paths <- ifelse(startsWith(meta$profile, "/"), meta$profile,file.path(dirname(manifest),meta$profile))
  paths <- vapply(paths,normalizePath,"",mustWork=TRUE,USE.NAMES=FALSE)
  if (anyDuplicated(paths)) stop("Repeated sample profile")
  matrices <- list()
  for (ch in c("EM","BPM")) meta[[paste0(ch,"_valid_fraction")]] <- NA_real_
  for (i in seq_along(ids)) {
    x <- data.table::fread(file.path(paths[i],"cgi.tsv"),data.table=FALSE)
    needed <- c("chr","start","end","nReads",channels[1:5])
    if (!all(needed %in% names(x))) stop("Invalid CGI table for ",ids[i])
    x <- x[x$chr %in% paste0("chr",1:22),,drop=FALSE]
    loci <- paste(x$chr,x$start,x$end,sep=".")
    if (!length(loci) || anyDuplicated(loci)) stop("Invalid or duplicate autosomal CGI coordinates")
    m <- read_tab(file.path(paths[i],"motifs.tsv"))
    if (!all(c("kmer","EM","BPM") %in% names(m)) || nrow(m)!=15625L || anyDuplicated(m$kmer) ||
        any(!grepl("^[ACGTN]{6}$",m$kmer))) stop("Expected all 15,625 A/C/G/T/N 6-mers, including unobserved motifs")
    for (key in channels) {
      motif <- key %in% c("EM","BPM")
      nm <- if (motif) m$kmer else loci
      val <- if (motif) as.numeric(m[[key]]) else as.numeric(x[[key]])
      if (motif) {
        if (any(!is.finite(val) | val<0 | val!=floor(val)) || sum(val)==0) stop("Invalid or empty ",key," counts for ",ids[i])
        canonical <- !grepl("N",nm,fixed=TRUE)
        total <- sum(val)
        val <- val[canonical]; nm <- nm[canonical]
        if (!sum(val)) stop("No unambiguous ",key," motifs for ",ids[i])
        meta[[paste0(key,"_valid_fraction")]][i] <- sum(val)/total
        val <- val/sum(val)
        val[val==0] <- NA_real_  # absent 6-mers were NA in the study count matrix
      } else {
        val[x$nReads==0 | is.nan(val)] <- NA_real_
        if (any(!is.na(val) & !is.finite(val))) stop("Invalid ",key," values")
      }
      if (i==1L) matrices[[key]] <- matrix(NA_real_,length(ids),length(nm),dimnames=list(ids,nm))
      ix <- match(colnames(matrices[[key]]),nm)
      if (length(ix)!=length(nm) || anyNA(ix)) stop("Inconsistent feature coordinates in ",ids[i])
      matrices[[key]][i,] <- val[ix]
    }
  }
  in_new_directory(outfile,function(dir) {
    for (key in channels) saveRDS(matrices[[key]],file.path(dir,paste0(key,".rds")),compress=FALSE)
    meta$profile <- NULL
    meta$sample_id <- ids
    write_tab(meta,file.path(dir,"samples.tsv"))
  })
}
load_features <- function(datadir,channel,ids,features=NULL,impute=TRUE) {
  x <- readRDS(file.path(datadir,paste0(channel,".rds")))
  if (!is.matrix(x) || !is.numeric(x) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      is.null(rownames(x)) || is.null(colnames(x)) || !all(ids %in% rownames(x))) stop("Invalid ",channel," matrix")
  if (is.null(features)) features <- colnames(x)
  absent <- setdiff(features,colnames(x))
  if (length(absent)) stop(channel," is missing ",length(absent)," trained features, e.g. ",absent[1L],". Check hg19 coordinates and BED.")
  x <- x[ids,features,drop=FALSE]
  if (impute) x[is.na(x)] <- 0
  x
}
scoring_features <- function(datadir,channel,ids,features,cohort,meta) {
  x <- load_features(datadir,channel,ids,features)
  if (cohort=="blind" && channel %in% c("EM","BPM")) {
    f <- meta[ids,paste0(channel,"_valid_fraction")]
    if(length(f)!=length(ids) || any(!is.finite(f) | f<=0 | f>1))
      stop("Missing valid-motif denominator for ",channel)
    x <- x*f
  }
  x
}
sex_adjust <- function(p,sex) {
  if(length(sex)!=nrow(p)) stop("Need one sex code per probability row")
  sex <- study_sex(sex)
  p[sex=="M",intersect(c("Ovary","Breast"),colnames(p))] <- 0
  p[sex=="F",intersect("Prostate",colnames(p))] <- 0
  s <- rowSums(p)
  if (any(!is.finite(s) | s<=0)) stop("Invalid TOO probability sum")
  p/s
}
read_samples <- function(dir,labelled=FALSE) {
  x <- read_tab(file.path(dir,"samples.tsv"))
  if(!all(c("sample_id","sex") %in% names(x)) || !nrow(x) || anyNA(x$sample_id) ||
     any(!nzchar(x$sample_id)) || anyDuplicated(x$sample_id)) stop("Invalid samples.tsv")
  x$sample_id <- as.character(x$sample_id); x$sex <- study_sex(x$sex)
  if("participant_id" %in% names(x) && (anyNA(x$participant_id) || any(!nzchar(x$participant_id)) || anyDuplicated(x$participant_id)))
    stop("Supply one sample per participant and nonempty participant_id values")
  if(labelled && (!"group" %in% names(x) || any(!x$group %in% c("Normal",tissues)))) stop("Training requires a known group for every sample")
  rownames(x) <- x$sample_id
  x
}
import_tables <- function(source,dest) {
  meta <- read_samples(source)
  in_new_directory(dest,function(out) {
    for(ch in channels) {
      d <- read_tab(file.path(source,paste0(ch,".tsv.gz")))
      if(!"sample_id" %in% names(d) || anyDuplicated(names(d)) ||
         anyNA(d$sample_id) || anyDuplicated(d$sample_id) || !setequal(d$sample_id,meta$sample_id))
        stop("Invalid sample identifiers or feature columns in ",ch)
      d <- d[match(meta$sample_id,d$sample_id),,drop=FALSE]
      x <- as.matrix(d[,setdiff(names(d),"sample_id"),drop=FALSE])
      rownames(x) <- d$sample_id
      if(!is.numeric(x) || ncol(x)<2L || any(!is.na(x) & !is.finite(x))) stop("Invalid numeric matrix in ",ch)
      if(ch %in% c("EM","BPM")) {
        if(any(x<0,na.rm=TRUE) || any(abs(rowSums(x,na.rm=TRUE)-1)>1e-6))
          stop("EM and BPM tables must contain all 4096 canonical frequencies, summing to one")
        if(ncol(x)!=4096L || any(!grepl("^[ACGT]{6}$",colnames(x)))) stop("Expected 4096 distinct A/C/G/T hexamers")
        f <- meta[[paste0(ch,"_valid_fraction")]]
        if(length(f)!=nrow(meta) || any(!is.finite(f) | f<=0 | f>1))
          stop("Supply the valid-motif count fraction for ",ch)
        x[x==0] <- NA_real_
      }
      saveRDS(x,file.path(out,paste0(ch,".rds")),compress=FALSE)
    }
    write_tab(meta,file.path(out,"samples.tsv"))
  })
}
