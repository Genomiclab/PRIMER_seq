#!/usr/bin/env Rscript
# Check the matrix interface against the formula used in the original scripts.
set.seed(204L)
x <- matrix(rnorm(60L*12L),60L,12L,dimnames=list(paste0('SIM_',1:60),paste0('feature',1:12)))
for(k in c(2L,3L)) {
  y <- factor(rep(paste0('class',1:k),length.out=nrow(x)))
  d <- as.data.frame(x);d$Group <- y
  set.seed(73L)
  formula <- e1071::svm(Group~.,data=d,type='C-classification',kernel='linear',scale=TRUE,probability=TRUE)
  set.seed(73L)
  direct <- e1071::svm(x=x,y=y,type='C-classification',kernel='linear',scale=TRUE,probability=TRUE)
  a <- attr(predict(formula,as.data.frame(x),probability=TRUE),'probabilities')
  b <- attr(predict(direct,as.data.frame(x),probability=TRUE),'probabilities')
  stopifnot(identical(formula$index,direct$index),identical(formula$coefs,direct$coefs),
            isTRUE(all.equal(a,b,tolerance=1e-12)))
}
cat('Binary and multiclass SVM input interfaces agree.\n')
