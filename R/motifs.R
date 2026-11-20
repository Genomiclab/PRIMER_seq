# The study's endogenous-end calculation. Coordinates below are one-based.
count_motifs <- function(bam, reference = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chunk = 250000L) {
  grid <- expand.grid(rep(list(c("A","T","C","G","N")), 6L), stringsAsFactors = FALSE)
  kmers <- do.call(paste0, rev(grid))
  em <- bpm <- numeric(length(kmers))
  lengths <- GenomeInfoDb::seqlengths(reference)
  targets <- paste0("chr", 1:22)
  declared <- Rsamtools::scanBamHeader(bam)[[1L]]$targets
  common <- intersect(targets,names(declared))
  if (!length(common) || any(declared[common] != lengths[common]))
    stop("BAM autosomal names and lengths must match UCSC hg19")
  flags <- Rsamtools::scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
      isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE, isSecondMateRead=TRUE)
  param <- Rsamtools::ScanBamParam(flag=flags, what=c("rname","flag","pos","cigar"))
  f <- Rsamtools::BamFile(bam, yieldSize=as.integer(chunk))
  open(f)
  on.exit(close(f))
  repeat {
    x <- Rsamtools::scanBam(f, param=param)[[1L]]
    if (!length(x$pos)) break
    chr <- as.character(x$rname)
    rev_read <- bitwAnd(x$flag, 16L) != 0L
    start <- x$pos
    end <- start + GenomicAlignments::cigarWidthAlongReferenceSpace(x$cigar) - 1L
    # R2 on the minus strand uses the reference to the right of its last base.
    # R2 on the plus strand uses the reference to the left of its first base.
    em_start <- ifelse(rev_read, end + 1L, start - 6L)
    bp_start <- ifelse(rev_read, end + 4L, start - 9L)
    site_start <- ifelse(rev_read, end + 4L, start - 7L)
    keep <- which(chr %in% targets & pmin(em_start, bp_start, site_start) >= 1L &
                  pmax(em_start+5L, bp_start+5L, site_start+3L) <= lengths[chr])
    if (!length(keep)) next
    seq_at <- function(index, pos, width, orient=FALSE) {
      s <- Biostrings::getSeq(reference, chr[index], start=as.integer(pos[index]),
                              width=as.integer(width))
      if (orient) s[rev_read[index]] <- Biostrings::reverseComplement(s[rev_read[index]])
      as.character(s)
    }
    keep <- keep[seq_at(keep, site_start, 4L) != "CCGG"]
    if (!length(keep)) next
    em <- em + tabulate(match(seq_at(keep, em_start, 6L, TRUE), kmers), length(kmers))
    bpm <- bpm + tabulate(match(seq_at(keep, bp_start, 6L, TRUE), kmers), length(kmers))
  }
  data.frame(kmer=kmers, EM=em, BPM=bpm)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args)!=2L) stop("Usage: Rscript motifs.R input.bam motifs.tsv")
  if (file.exists(args[2L])) stop("Output already exists")
  write.table(count_motifs(args[1L]), args[2L], sep="\t", quote=FALSE, row.names=FALSE)
}
