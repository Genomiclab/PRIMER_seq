auroc <- function(score,case) {
  n1 <- sum(case); n0 <- sum(!case)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score,ties.method="average")[case])-n1*(n1+1)/2)/(n1*n0)
}
evaluate <- function(predictions,truth,dest) {
  p <- read_tab(predictions); y <- read_tab(truth)
  if (!all(c("sample_id","score","called") %in% names(p)) || !all(c("sample_id","group") %in% names(y))) stop("Check prediction and truth column names")
  if (!nrow(p) || anyNA(p$sample_id) || anyNA(y$sample_id) || anyDuplicated(p$sample_id) || anyDuplicated(y$sample_id) || !setequal(p$sample_id,y$sample_id)) stop("Predictions and truth must contain the same unique IDs")
  y <- y[match(p$sample_id,y$sample_id),,drop=FALSE]
  if (any(!y$group %in% c("Normal",tissues)) || any(!is.finite(p$score) | p$score<0 | p$score>1)) stop("Invalid group or score")
  classes <- tissues[tissues %in% names(p)]
  pos <- as.logical(p$called); case <- y$group!="Normal"
  if (anyNA(pos) || !any(case) || !any(!case)) stop("Need calls, cancer cases and controls")
  tp <- sum(pos & case); tn <- sum(!pos & !case)
  summary <- data.frame(n=nrow(p),TP=tp,FN=sum(case)-tp,TN=tn,FP=sum(!case)-tn,
     sensitivity=tp/sum(case),specificity=tn/sum(!case),AUC=auroc(p$score,case))
  get_sens <- function(name,ix) data.frame(stratum=name,n=sum(ix),detected=sum(ix & pos),sensitivity=if(any(ix)) mean(pos[ix]) else NA_real_)
  strata <- do.call(rbind,lapply(unique(y$group[y$group!="Normal"]),function(t)get_sens(t,y$group==t)))
  if ("stage" %in% names(y)) for(s in c("I","II","III","IV","I-II")) {
    ix <- y$stage %in% (if(s=="I-II") c("I","II") else s)
    strata <- rbind(strata,get_sens(paste("Stage",s),case & ix))
  }
  cuts <- order(p$score,decreasing=TRUE)
  ends <- !duplicated(p$score[cuts],fromLast=TRUE)
  roc <- data.frame(fpr=c(0,cumsum(!case[cuts])[ends]/sum(!case)),tpr=c(0,cumsum(case[cuts])[ends]/sum(case)))
  in_new_directory(dest,function(dir) {
    write_tab(summary,file.path(dir,"detection.tsv"))
    write_tab(strata,file.path(dir,"sensitivity.tsv"))
    write_tab(roc,file.path(dir,"roc.tsv"))
    grDevices::pdf(file.path(dir,"performance.pdf"),width=5.6,height=5.3)
    on.exit(grDevices::dev.off(),add=TRUE)
    plot(roc$fpr,roc$tpr,type="s",xlim=c(0,1),ylim=c(0,1),asp=1,
         xlab="False-positive rate",ylab="Sensitivity",lwd=2,col="#276B78",
         main=sprintf("Cancer detection   |   AUC %.3f",summary$AUC))
    abline(0,1,lty=3,col="grey65")
    points(1-summary$specificity,summary$sensitivity,pch=19,col="#A04A32")
    if(length(classes)>=2L && tp>0) {
      q <- as.matrix(p[case & pos,classes,drop=FALSE]); g <- y$group[case & pos]
      if (any(!is.finite(q)) || any(q<0) || any(abs(rowSums(q)-1)>1e-5)) stop("Invalid TOO probabilities for detected cancer cases")
      ranks <- t(vapply(seq_len(nrow(q)),function(i) classes[order(-q[i,],seq_along(classes))],character(length(classes))))
      k <- seq_len(min(5L,length(classes))); correct <- vapply(k,function(j) sum(rowSums(ranks[,seq_len(j),drop=FALSE]==g)>0),integer(1))
      write_tab(data.frame(k=k,correct=correct,total=tp,accuracy=correct/tp),file.path(dir,"origin_top_k.tsv"))
      tab <- table(truth=factor(g,levels=classes),prediction=factor(ranks[,1L],levels=classes))
      write_tab(as.data.frame(tab),file.path(dir,"origin_confusion.tsv"))
      barplot(correct/tp,names.arg=paste("Top",k),ylim=c(0,1),ylab="Accuracy among detected cancer cases",col="#276B78",border=NA,main="Tissue of origin")
    }
  })
  print(summary,row.names=FALSE)
}
