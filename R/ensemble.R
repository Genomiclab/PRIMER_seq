# Ordered feature selection, out-of-fold probabilities and the two independent stacks.
select_motifs <- function(x,group,task,classes) {
  if (task=="origin") {
    keep <- group!="Normal"
    x <- x[keep,,drop=FALSE]; group <- group[keep]
  }
  x <- x[,colSums(is.na(x))<50L,drop=FALSE]
  if (task=="detection") x[is.na(x)] <- 0
  contrasts <- if (task=="detection") lapply(classes,function(t)c("Normal",t)) else combn(classes,2,simplify=FALSE)
  selected <- character()
  for (pair in contrasts) {
    a <- x[group==pair[1L],,drop=FALSE]; b <- x[group==pair[2L],,drop=FALSE]
    p <- vapply(seq_len(ncol(x)),function(j) {
      u <- a[,j]; v <- b[,j]
      if (!any(is.finite(u)) || !any(is.finite(v))) return(NA_real_)
      suppressWarnings(stats::wilcox.test(u,v)$p.value)
    },numeric(1))
    selected <- union(selected,colnames(x)[which(p.adjust(p,"BH")<0.01)])
  }
  if (length(selected)<2L) stop("Too few motifs passed the training-data filter for ",task)
  selected
}
model_file <- function(dir,task,channel,algo,fold) {
  tag <- if(task=="detection") "cancer_detection" else "TOO"
  cl <- if(task=="detection") "Binaryclass" else "Multiclass"
  if (channel=="stack") name <- sprintf("Level2_multimodel_%s_model_%d_%s_%s.RData",tag,fold,cl,algo)
  else name <- sprintf("Level1_multimodel_%s_%s_model_%s_%d_%s_%s.RData",
          if(channel %in% c("EM","BPM")) "motif" else "CGI",tag,archive_channel(channel),fold,cl,algo)
  file.path(dir,"weights",name)
}
read_model <- function(path) {
  e <- new.env(parent=emptyenv()); keys <- load(path,e)
  if (!"model" %in% keys) stop("No fitted model in ",basename(path))
  e$model
}
fit_folds <- function(x,group,task,paths,algo,threads,classes,preset) {
  cl <- if(task=="detection") c("Normal","Cancer") else classes
  strata <- factor(group,levels=if(task=="detection") c("Normal",classes) else classes)
  label <- if(task=="detection") ifelse(group=="Normal","Normal","Cancer") else group
  y <- factor(label,levels=cl)
  folds <- caret::createFolds(strata,k=5L)
  out <- matrix(NA_real_,nrow(x),if(task=="detection") 1L else length(cl),
                dimnames=list(rownames(x),if(task=="detection") "Cancer" else cl))
  fold_id <- integer(nrow(x))
  for(i in seq_along(folds)) {
    hold <- folds[[i]]
    model <- fit_learner(algo,x[-hold,,drop=FALSE],y[-hold],threads,preset)
    out[hold,] <- predict_learner(algo,model,x[hold,,drop=FALSE],cl)
    save(model,file=paths[i])
    fold_id[hold] <- i
  }
  if (anyNA(out) || any(fold_id==0L)) stop("Incomplete out-of-fold predictions")
  list(p=out,fold=fold_id)
}
average_folds <- function(x,task,modeldir,channel,algo,classes) {
  cl <- if (task=="origin") classes else c("Normal","Cancer")
  total <- NULL
  for(i in seq_len(5L)) {
    model <- read_model(model_file(modeldir,task,channel,algo,i))
    p <- predict_learner(algo,model,x,cl)
    total <- if(is.null(total)) p else total+p
  }
  total/5
}
stack_names <- function(p,channel,algo,task) {
  nm <- paste(archive_channel(channel),algo,sep="_")
  colnames(p) <- if(task=="detection") paste0(nm,"_score") else paste(nm,colnames(p),"score",sep="_")
  p
}
stack_input <- function(p,algos) {
  # The manuscript scripts merge methylation probabilities first, then the motifs.
  order <- unlist(lapply(list(channels[1:5],channels[6:7]),function(block)
    unlist(lapply(algos,function(a)paste(a,block,sep="/")),use.names=FALSE)),use.names=FALSE)
  do.call(cbind,unname(p[order]))
}
fit_study <- function(datadir,marker_path,dest,algos=algorithms,threads=1L,fixed_motifs=TRUE,preset="study") {
  samples <- read_samples(datadir,TRUE)
  classes <- tissues[tissues %in% samples$group]
  if(length(classes)<2L || !"Normal" %in% samples$group) stop("Training requires controls and at least two cancer types")
  if(any(table(factor(samples$group,levels=c("Normal",classes)))<5L)) stop("Five-fold training requires at least five samples per group")
  if(!preset %in% c("study","demo")) stop("Choose study or demo preset")
  if(!length(algos) || anyDuplicated(algos) || any(!algos %in% algorithms) || is.na(threads) || threads<1L) stop("Invalid learner or thread count")
  markers <- read_tab(marker_path)
  if(!all(c("feature","task","channel") %in% names(markers))) stop("Invalid marker table")
  if(any(!markers$task %in% c("detection","origin")) || any(!markers$channel %in% channels) ||
     anyNA(markers$feature) || any(!nzchar(markers$feature))) stop("Invalid marker task, channel or feature")
  in_new_directory(dest,function(out) {
    dir.create(file.path(out,"weights")); dir.create(file.path(out,"training"))
    spec <- list(algorithms=algos,features=list(),tissues=classes,version=2L,preset=preset,
                 seeds=c(detection_level1=35L,origin_level1=21L,detection_stack=2L,origin_stack=100L))
    for(task in c("detection","origin")) {
      ids <- sort(samples$sample_id[task=="detection" | samples$group!="Normal"])
      features <- setNames(vector("list",length(channels)),channels)
      for(ch in channels) {
        if (ch %in% c("EM","BPM") && !fixed_motifs) {
          x <- load_features(datadir,ch,ids,impute=FALSE)
          features[[ch]] <- select_motifs(x,samples[ids,"group"],task,classes)
        } else features[[ch]] <- as.character(markers$feature[markers$task==task & markers$channel==ch])
        if(length(features[[ch]])<2L || anyDuplicated(features[[ch]])) stop("Invalid feature list for ",task,"/",ch)
      }
      spec$features[[task]] <- features
      probs <- list(); foldmap <- list()
      for(a in algos) {
        set.seed(if(task=="detection") 35L else 21L)
        for(ch in channels) {
          message(task, " / ",ch," / ",a)
          ordered <- ids
          if (!ch %in% c("EM","BPM") && "methylation_block" %in% names(samples))
            ordered <- ids[order(samples[ids,"methylation_block"],ids)]
          x <- load_features(datadir,ch,ordered,features[[ch]])
          g <- samples[ordered,"group"]
          paths <- vapply(1:5,function(i)model_file(out,task,ch,a,i),"")
          res <- fit_folds(x,g,task,paths,a,threads,classes,preset)
          if(task=="origin") res$p <- sex_adjust(res$p,samples[ordered,"sex"])
          key <- paste(a,ch,sep="/")
          probs[[key]] <- stack_names(res$p[ids,,drop=FALSE],ch,a,task)
          foldmap[[key]] <- setNames(res$fold,ordered)[ids]
        }
      }
      x <- stack_input(probs,algos)
      set.seed(if(task=="detection") 2L else 100L)
      res <- fit_folds(x,samples[ids,"group"],task,
             vapply(1:5,function(i)model_file(out,task,"stack","SVM_linear",i),""),"SVM_linear",threads,classes,preset)
      foldmap$stack <- res$fold
      write_tab(data.frame(sample_id=ids,as.data.frame(foldmap),check.names=FALSE),
                file.path(out,"training",paste0(task,"_folds.tsv")))
      if(task=="detection") {
        spec$threshold <- unname(quantile(res$p[samples[ids,"group"]=="Normal",1L],0.97,type=1))
        oof <- data.frame(sample_id=ids,score=res$p[,1L],called=res$p[,1L]>spec$threshold)
      } else oof <- data.frame(sample_id=ids,sex_adjust(res$p,samples[ids,"sex"]),check.names=FALSE)
      write_tab(oof,file.path(out,"training",paste0(task,"_oof.tsv")))
    }
    write_tab(samples,file.path(out,"training","samples.tsv"))
    write_tab(markers,file.path(out,"markers_input.tsv"))
    saveRDS(spec,file.path(out,"spec.rds"))
    write_tab(data.frame(threshold=spec$threshold),file.path(out,"threshold.tsv"))
    writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
  })
}
score_study <- function(datadir,modeldir,dest,cohort) {
  if(!cohort %in% c("validation","blind")) stop("Choose validation or blind")
  meta <- read_samples(datadir); ids <- meta$sample_id
  spec <- readRDS(file.path(modeldir,"spec.rds"))
  classes <- spec$tissues
  if(length(classes)<2L || any(!classes %in% tissues) || anyDuplicated(classes) || length(spec$threshold)!=1L || !is.finite(spec$threshold)) stop("Invalid model specification")
  train_ids <- file.path(modeldir,"training","detection_oof.tsv")
  if(file.exists(train_ids) && any(ids %in% read_tab(train_ids)$sample_id)) stop("Scoring cohort overlaps training")
  training <- read_tab(file.path(modeldir,"training","samples.tsv"))
  if("participant_id" %in% names(training)) {
    if(!"participant_id" %in% names(meta)) stop("Supply participant_id to check separation from the training participants")
    if(any(meta$participant_id %in% training$participant_id)) stop("Scoring participants overlap training")
  }
  score_task <- function(task,rows) {
    probs <- list()
    for(ch in channels) {
      x <- scoring_features(datadir,ch,rows,spec$features[[task]][[ch]],cohort,meta)
      for(a in spec$algorithms) {
        p <- average_folds(x,task,modeldir,ch,a,classes)
        if(task=="origin") p <- sex_adjust(p,meta[rows,"sex"])
        probs[[paste(a,ch,sep="/")]] <- stack_names(p,ch,a,task)
      }
    }
    p <- average_folds(stack_input(probs,spec$algorithms),task,modeldir,"stack","SVM_linear",classes)
    if(task=="origin") sex_adjust(p,meta[rows,"sex"]) else p
  }
  p <- score_task("detection",ids)[,1L]
  positive <- if(cohort=="blind") p>=spec$threshold else p>spec$threshold
  out <- data.frame(sample_id=ids,score=unname(p),called=positive,threshold=spec$threshold)
  origin <- matrix(NA_real_,length(ids),length(classes),dimnames=list(ids,classes))
  if(any(positive)) origin[positive,] <- score_task("origin",ids[positive])
  out <- cbind(out,origin)
  in_new_directory(dest,function(d) write_tab(out,file.path(d,"predictions.tsv")))
}
