# Argument parsing shared by the numbered command-line steps.
read_options <- function(required,defaults=list(),usage) {
  values <- commandArgs(TRUE)
  if(length(values)==0L || any(values %in% c("-h","--help"))) {
    cat(usage,"\n")
    quit(status=0L)
  }
  if(length(values)%%2L) stop(usage)
  keys <- sub("^--","",values[seq(1L,length(values),2L)])
  if(anyDuplicated(keys) || any(!keys %in% c(required,names(defaults))))
    stop("Unknown or repeated option. ",usage)
  ans <- defaults
  for(i in seq_along(keys)) ans[[keys[i]]] <- values[i*2L]
  if(any(!required %in% names(ans))) stop(usage)
  ans
}
