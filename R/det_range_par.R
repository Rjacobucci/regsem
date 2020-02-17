#'
#'
#' Determine the initial range for stability selection, parallel version
#'
#' This function perform regsem on bootstrap samples to determine the initial range for stability selection. Interquartile range of the bootstrap optimal regularization amounts are uesd as the final range. Parallelization is used to achieve faster performance.
#' @param data data frame
#' @param model lavaan output object.
#' @param times number of bootstrap samples used.
#' @param ... Any additional arguments to pass to regsem() or cv_regsem().
#' @return result the lambda values and the upper bound and lower bound of the interquartile range.
#' @export


det_range_par<-function(data,
                        model,
                        times=50,
                        ...){
  #library(parallel)
  nsize<-dim(data)[1]
  num_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(num_cores)
  parallel::clusterEvalQ(cl, {
    #library(regsem)
  })
  lam<-parallel::parSapply(cl,1:times,function(i,model,data,nsize,...){
    ids = sample(1:nsize,nsize,replace=T)
    datasub.boot <- data[ids,]
    est_model_boot <- sem(model, data = datasub.boot)
    lambda<-rep(NA,length(coef(est_model_boot)))
    try(cv.out.boot <- cv_regsem(est_model_boot,...))
    try(lambda<-cv.out.boot$fit[which.min(cv.out.boot$fits[cv.out.boot$fit[,2]==0,4])])
    return(lambda)},model,data,nsize,...)
  stopCluster(cl)
  
  if (sum(lam!=0,na.rm=T)==0){
    warning("0 penalty is selected by all runs.")
    lb=0;ub=0;
  }else{
    lb<-quantile(lam[lam!=0],0.25,na.rm=T)
    ub<-quantile(lam[lam!=0],0.75,na.rm=T)
  }
  Lam<-c(lb,ub)
  result<-list()
  #result$times<-times
  #result$n.lambda<-n.lambda
  #result$jump<-jump
  result$lambdas<-lam
  result$lb<-lb
  result$ub<-ub
  result$zero_removed<-min(lam,na.rm=T)==0
  return(result)
}