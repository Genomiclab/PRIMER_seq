#!/usr/bin/env Rscript
# Generate hg19 paired alignments with simulated bisulfite conversion and CpG methylation.
args <- commandArgs(TRUE)
if(length(args)!=1L || file.exists(args[1L]))
  stop("Usage: Rscript examples/make_demo_bam.R NEW_DIRECTORY")
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R","io.R"))
ref <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
markers <- read_tab(file.path(root,"examples","three_cancers","markers.tsv"))
regions <- unique(markers$feature[!markers$channel %in% c("EM","BPM")])
parts <- strsplit(regions,".",fixed=TRUE)
bed <- data.frame(chr=vapply(parts,"[","",1L),start=as.integer(vapply(parts,"[","",2L)),
                  end=as.integer(vapply(parts,"[","",3L)))
targets <- paste0("chr",1:22)
bed <- bed[order(match(bed$chr,targets),bed$start),]
settings <- read_tab(file.path(root,"examples","simulation.tsv"))
setting <- setNames(settings$value,settings$parameter)
set.seed(as.integer(setting[["seed"]]))
in_new_directory(args[1L],function(out) {
  sam <- file.path(out,"synthetic.sam")
  stream <- file(sam,"wt"); on.exit(close(stream),add=TRUE)
  lengths <- GenomeInfoDb::seqlengths(ref)[targets]
  writeLines(c("@HD\tVN:1.6\tSO:unsorted",sprintf("@SQ\tSN:%s\tLN:%d",targets,as.integer(lengths))),stream)
  cpgs <- list(); pair_id <- 0L
  for(i in seq_len(nrow(bed))) {
    region <- bed[i,]
    lo <- max(1L,region$start-1000L); hi <- min(lengths[region$chr],region$end+1000L)
    seq <- as.character(Biostrings::getSeq(ref,region$chr,start=as.integer(lo),end=as.integer(hi)))
    sites <- gregexpr("CG",seq,fixed=TRUE)[[1L]]
    sites <- sites[sites>0L]+lo-1L
    cpgs[[i]] <- data.frame(chr=region$chr,start=sites,end=sites+1L)
    for(j in seq_len(as.integer(setting[["pairs_per_region"]]))) {
      frag <- max(as.integer(setting[["read_length"]])+10L,
                  as.integer(round(rnorm(1L,setting[["fragment_length"]],setting[["fragment_sd"]]))))
      positions <- seq.int(region$start+1L,max(region$start+1L,region$end-frag))
      left <- as.integer(positions[sample.int(length(positions),1L)])
      n <- as.integer(setting[["read_length"]]); right <- left+frag-n
      minus <- j%%2L==0L
      retained <- sites[runif(length(sites))<runif(1L,.15,.85)]
      pair_id <- pair_id+1L
      for(mate in 1:2) {
        pos <- if((mate==1L)!=minus) left else right
        next_pos <- if(pos==left) right else left
        bases <- strsplit(substr(seq,pos-lo+1L,pos-lo+n),"",fixed=TRUE)[[1L]]
        coords <- pos+seq_along(bases)-1L
        protected <- if(minus) retained+1L else retained
        target <- if(minus) "G" else "C"; converted <- if(minus) "A" else "T"
        bases[bases==target & !coords %in% protected] <- converted
        flag <- if(minus) c(83L,163L)[mate] else c(99L,147L)[mate]
        row <- c(sprintf("SIMFRAG%06d",pair_id),flag,region$chr,pos,60L,paste0(n,"M"),"=",next_pos,
                 if(pos==left) frag else -frag,paste0(bases,collapse=""),strrep("I",n))
        writeLines(paste(row,collapse="\t"),stream)
      }
    }
  }
  flush(stream)
  Rsamtools::asBam(sam,destination=file.path(out,"synthetic"),overwrite=FALSE,indexDestination=TRUE)
  annotation <- unique(do.call(rbind,cpgs))
  annotation <- annotation[order(annotation$chr,annotation$start),]
  txt <- file.path(out,"CpG.tsv")
  write.table(annotation,txt,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  Rsamtools::bgzip(txt,file.path(out,"hg19_demo_CpG.gz"),overwrite=FALSE)
  Rsamtools::indexTabix(file.path(out,"hg19_demo_CpG.gz"),seq=1L,start=2L,end=3L)
  write.table(bed,file.path(out,"regions.bed"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write_tab(data.frame(pairs=pair_id,regions=nrow(bed),read_length=setting[["read_length"]]),
            file.path(out,"summary.tsv"))
})
