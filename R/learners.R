# Model specifications follow the five active branches of the study's shared R functions.
fit_learner <- function(a,x,y,threads=1L,preset="study") {
  if (!is.factor(y) || anyNA(y)) stop("Invalid class labels")
  multi <- nlevels(y)>2L
  if (a=="LR") {
    response <- as.integer(y)-1L
    return(glmnet::cv.glmnet(x,response,family=if(multi) "multinomial" else "binomial",alpha=0.5,nfolds=if(preset=="demo") 5L else 10L))
  }
  if (a=="GBM") {
    if (threads>1L) {
      cl <- parallel::makeCluster(threads)
      doParallel::registerDoParallel(cl)
      on.exit({parallel::stopCluster(cl); foreach::registerDoSEQ()})
    }
    control <- caret::trainControl(allowParallel=threads>1L,verboseIter=FALSE)
    options <- list(x=as.data.frame(x),y=y,method="gbm",trControl=control,verbose=FALSE)
    if(preset=="demo") {
      options$trControl <- caret::trainControl(method="cv",number=3L,allowParallel=threads>1L)
      options$tuneGrid <- data.frame(n.trees=75L,interaction.depth=2L,shrinkage=0.1,n.minobsinnode=5L)
    }
    return(do.call(caret::train,options))
  }
  if (a=="SVM_linear") {
    df <- as.data.frame(x); df$Group <- y
    return(e1071::svm(Group~.,data=df,type="C-classification",kernel="linear",scale=TRUE,probability=TRUE))
  }
  if (a=="XGBoost") {
    p <- list(booster="gbtree",objective=if(multi) "multi:softprob" else "binary:logistic",
              eval_metric=if(multi) "mlogloss" else "logloss",nthread=threads)
    if (multi) p$num_class <- nlevels(y)
    return(xgboost::xgb.train(p,xgboost::xgb.DMatrix(x,label=as.integer(y)-1L),nrounds=if(preset=="demo") 100L else 1000L,verbose=0))
  }
  if (a=="Catboost") {
    # Historical labels were the one-based factor codes. Use the same ordering.
    pool <- catboost::catboost.load_pool(x,label=as.integer(y))
    return(catboost::catboost.train(pool,NULL,params=list(loss_function=if(multi) "MultiClass" else "Logloss",
        iterations=if(preset=="demo") 150L else 2000L,thread_count=threads,verbose=0L,allow_writing_files=FALSE)))
  }
  stop("Unknown learner: ",a)
}
predict_learner <- function(a,model,x,classes) {
  pkg <- c(LR="glmnet", GBM="caret", SVM_linear="e1071", XGBoost="xgboost", Catboost="catboost")[[a]]
  if (!requireNamespace(pkg,quietly=TRUE)) stop("Install ",pkg)
  multi <- length(classes)>2L
  if (a=="SVM_linear") {
    p <- attr(stats::predict(model,as.data.frame(x),probability=TRUE),"probabilities")
  } else if (a=="LR") {
    p <- stats::predict(model,x,type="response",s=model$lambda.min)
    if (multi) p <- matrix(p,nrow=nrow(x),dimnames=list(NULL,classes))
  } else if (a=="GBM") {
    p <- stats::predict(model,as.data.frame(x),type="prob")
  } else if (a=="XGBoost") {
    p <- stats::predict(model,xgboost::xgb.DMatrix(x),reshape=TRUE)
    if (multi) p <- matrix(p,nrow=nrow(x),dimnames=list(NULL,classes))
  } else if (a=="Catboost") {
    p <- catboost::catboost.predict(model,catboost::catboost.load_pool(x),prediction_type="Probability")
    if (multi) p <- matrix(p,nrow=nrow(x),dimnames=list(NULL,classes))
  } else stop("Unknown learner: ",a)
  p <- as.matrix(p)
  if (!multi && ncol(p)==1L) {
    p <- cbind(1-p[,1L],p[,1L])
    colnames(p) <- classes
  }
  want <- if(identical(classes,c("Normal","Cancer"))) "Cancer" else classes
  if (nrow(p)!=nrow(x) || !all(want %in% colnames(p))) stop("Unexpected ",a," probability shape")
  p <- p[,want,drop=FALSE]
  if (any(!is.finite(p) | p<0 | p>1)) stop("Invalid probabilities from ",a)
  rownames(p) <- rownames(x)
  p
}
