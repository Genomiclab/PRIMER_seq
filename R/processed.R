# Matrix import for existing mHap-derived statistics and EM/BPM tables.
read_motif_matrix <- function(path,ids,layout="auto") {
  if(!layout %in% c("auto","samples-by-features","features-by-samples")) stop("Invalid matrix layout")
  if(grepl("[.]rds$",path,ignore.case=TRUE)) {
    x <- readRDS(path)
  } else if(grepl("[.](rdata|rda)$",path,ignore.case=TRUE)) {
    env <- new.env(parent=emptyenv()); keys <- load(path,env)
    if(length(keys)!=1L) stop("RData input must contain exactly one matrix or data frame")
    x <- env[[keys]]
  } else {
    comma <- grepl("[.]csv([.]gz)?$",path,ignore.case=TRUE)
    sep <- if(comma) "," else "\t"; quote <- if(comma) "\"" else ""
    hdr <- read.table(path,header=TRUE,nrows=0L,sep=sep,quote=quote,comment.char="",row.names=1L,check.names=FALSE)
    x <- read.table(path,header=TRUE,sep=sep,check.names=FALSE,stringsAsFactors=FALSE,
                    quote=quote,comment.char="",row.names=1L,colClasses=c("character",rep("numeric",ncol(hdr))))
  }
  if(is.data.frame(x) && "sample_id" %in% names(x)) {
    rownames(x) <- as.character(x$sample_id); x$sample_id <- NULL
  }
  if(is.data.frame(x)) x <- as.matrix(x)
  if(!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
     anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop("Expected a named numeric motif matrix: ",path)
  row_match <- all(ids %in% rownames(x)); col_match <- all(ids %in% colnames(x))
  if(layout=="auto") {
    if(row_match==col_match) stop("Cannot infer matrix orientation from sample keys: ",path)
    layout <- if(row_match) "samples-by-features" else "features-by-samples"
  }
  if(layout=="features-by-samples") x <- t(x)
  if(!all(ids %in% rownames(x))) stop("Motif matrix is missing requested samples")
  x <- x[ids,,drop=FALSE]
  if(any(!grepl("^[ACGTN]{6}$",colnames(x)))) stop("Motif columns must be six-base A/C/G/T/N sequences")
  canonical <- !grepl("N",colnames(x),fixed=TRUE)
  if(sum(canonical)!=4096L) stop("Supply all 4096 distinct canonical hexamers")
  if(any(!is.na(x) & (!is.finite(x) | x<0))) stop("Invalid motif values")
  x
}

normalize_motif_matrix <- function(x,format,fraction=NULL) {
  if(!format %in% c("counts","canonical-frequencies","all-frequencies")) stop("Invalid motif format")
  canonical <- !grepl("N",colnames(x),fixed=TRUE)
  n <- nrow(x)
  if(!is.null(fraction) && (length(fraction)!=n || any(!is.finite(fraction) | fraction<=0 | fraction>1)))
    stop("Valid fractions must lie in (0,1]")
  x[is.na(x)] <- 0
  total <- rowSums(x); valid <- rowSums(x[,canonical,drop=FALSE])
  if(any(total<=0 | valid<=0)) stop("Empty motif sample")
  if(format=="counts") {
    if(any(x!=floor(x))) stop("Raw motif counts must be integers")
  } else if(format=="canonical-frequencies" || any(!canonical)) {
    if(any(abs(total-1)>1e-6)) stop("Frequency rows must sum to one")
  } else if(any(total>1+1e-6)) stop("All-count frequencies cannot sum above one")
  if(format=="canonical-frequencies" && any(!canonical)) stop("Canonical-frequency input must have exactly 4096 columns")
  if(all(canonical)) {
    if(format=="all-frequencies") {
      f <- pmin(total,1)
      if(!is.null(fraction) && any(abs(f-fraction)>1e-6)) stop("Supplied valid fraction disagrees with motif frequencies")
    } else {
      if(is.null(fraction)) stop("Canonical-only input requires the original EM/BPM_valid_fraction in sample metadata; use 1 only if no ambiguous ends were counted")
      f <- fraction
    }
  } else {
    f <- valid/total
    if(!is.null(fraction) && any(abs(f-fraction)>1e-6)) stop("Supplied valid fraction disagrees with motif counts")
  }
  x <- x[,canonical,drop=FALSE]/valid
  x[x==0] <- NA_real_
  list(matrix=x,fraction=f)
}

read_cgi_files <- function(meta,paths) {
  matrices <- list(); ids <- meta$sample_id
  for(i in seq_along(ids)) {
    d <- data.table::fread(paths[i],data.table=FALSE)
    needed <- c("chr","start","end","nReads",channels[1:5])
    if(!all(needed %in% names(d))) stop("Invalid CGI table for ",ids[i])
    d <- d[d$chr %in% paste0("chr",1:22),,drop=FALSE]
    loci <- paste(d$chr,d$start,d$end,sep=".")
    if(!length(loci) || anyDuplicated(loci) || any(!is.finite(d$nReads) | d$nReads<0)) stop("Invalid CGI coordinates or coverage")
    for(ch in channels[1:5]) {
      value <- as.numeric(d[[ch]])
      value[d$nReads==0 | is.nan(value)] <- NA_real_
      if(any(!is.na(value) & !is.finite(value))) stop("Invalid methylation value")
      if(i==1L) matrices[[ch]] <- matrix(NA_real_,length(ids),length(loci),dimnames=list(ids,loci))
      ix <- match(colnames(matrices[[ch]]),loci)
      if(length(ix)!=length(loci) || anyNA(ix)) stop("CGI region sets differ across samples")
      matrices[[ch]][i,] <- value[ix]
    }
  }
  matrices
}

import_processed <- function(manifest,em,bpm,dest,format="counts",layout="auto") {
  meta <- read_tab(manifest)
  if(!all(c("sample_id","sex","cgi") %in% names(meta)) || !nrow(meta) || anyNA(meta$sample_id) ||
     any(!nzchar(meta$sample_id)) || anyDuplicated(meta$sample_id)) stop("CGI manifest requires unique sample_id, sex and cgi columns")
  meta$sex <- study_sex(meta$sex)
  if("group" %in% names(meta) && any(!is.na(meta$group) & !meta$group %in% c("", "Normal",tissues)))
    stop("Unknown diagnosis in group")
  if("participant_id" %in% names(meta) && (anyNA(meta$participant_id) || any(!nzchar(meta$participant_id)) || anyDuplicated(meta$participant_id)))
    stop("Use one specimen per participant")
  paths <- ifelse(startsWith(meta$cgi,"/"),meta$cgi,file.path(dirname(manifest),meta$cgi))
  paths <- vapply(paths,normalizePath,"",mustWork=TRUE,USE.NAMES=FALSE)
  if(anyDuplicated(paths)) stop("Repeated CGI input file")
  in_new_directory(dest,function(out) {
    methylation <- read_cgi_files(meta,paths)
    for(ch in names(methylation)) saveRDS(methylation[[ch]],file.path(out,paste0(ch,".rds")),compress=FALSE)
    rm(methylation)
    for(ch in c("EM","BPM")) {
      x <- read_motif_matrix(if(ch=="EM") em else bpm,meta$sample_id,layout)
      field <- paste0(ch,"_valid_fraction")
      values <- normalize_motif_matrix(x,format,meta[[field]])
      meta[[field]] <- values$fraction
      saveRDS(values$matrix,file.path(out,paste0(ch,".rds")),compress=FALSE)
    }
    meta$cgi <- NULL
    write_tab(meta,file.path(out,"samples.tsv"))
    write_tab(data.frame(motif_format=format,matrix_layout=layout),file.path(out,"input_format.tsv"))
  })
}
