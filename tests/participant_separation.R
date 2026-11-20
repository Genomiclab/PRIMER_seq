#!/usr/bin/env Rscript
# A changed sample key must not bypass participant separation at inference.
entry <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
if(!file.exists(entry)) entry <- gsub("~+~"," ",entry,fixed=TRUE)
root <- dirname(dirname(normalizePath(entry)))
source(file.path(root,"R","io.R"))
source(file.path(root,"R","ensemble.R"))
expect_error <- function(expr,pattern) {
  message <- tryCatch({force(expr);NA_character_},error=function(e)conditionMessage(e))
  if(is.na(message) || !grepl(pattern,message,fixed=TRUE)) stop("Expected error: ",pattern)
}
run_test <- function() {
  tmp <- tempfile("primer-participant-test-")
  dir.create(tmp); on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  model <- file.path(tmp,"model"); data <- file.path(tmp,"data")
  dir.create(model); dir.create(file.path(model,"training")); dir.create(data)
  saveRDS(list(tissues=c("Lung","Colorectum","Liver"),threshold=.9),file.path(model,"spec.rds"))
  write_tab(data.frame(sample_id="001"),file.path(model,"training","detection_oof.tsv"))
  write_tab(data.frame(sample_id="001",participant_id="0001",sex="U",group="Normal"),
            file.path(model,"training","samples.tsv"))
  write_tab(data.frame(sample_id="002",participant_id="0001",sex="U"),file.path(data,"samples.tsv"))
  expect_error(score_study(data,model,file.path(tmp,"pred"),"validation"),"participants overlap training")
  stopifnot(read_samples(data)$participant_id=="0001")
  write_tab(data.frame(sample_id="002",sex="U"),file.path(data,"samples.tsv"))
  expect_error(score_study(data,model,file.path(tmp,"pred"),"blind"),"Supply participant_id")
  write_tab(data.frame(sample_id=c("002","003"),participant_id=c("0002","0002"),sex=c("F","M")),
            file.path(data,"samples.tsv"))
  expect_error(read_samples(data),"one sample per participant")
  p <- matrix(c(.2,.3,.5),3L,3L,byrow=TRUE,dimnames=list(NULL,c("Lung","Colorectum","Liver")))
  stopifnot(isTRUE(all.equal(p,sex_adjust(p,c("F","M","U")))))
}
run_test()
cat("Participant separation, text identifiers and three-class sex adjustment passed.\n")
