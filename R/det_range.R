#'
#'
#' Determine the initial range for stability selection
#'
#' This function perform regsem on bootstrap samples to determine the initial range for stability selection. Interquartile range of the bootstrap optimal regularization amounts are uesd as the final range.
#' @param data data frame
#' @param model lavaan output object.
#' @param times number of bootstrap samples used.
#' @param ... Any additional arguments to pass to regsem() or cv_regsem().
#' @return result the lambda values and the upper bound and lower bound of the interquartile range.
#' @export


det_range<-function(data,
                    model,
                    times=50,
                    ...){
  nsize<-dim(data)[1]
  lam=rep(NA,times)
  for(i in 1:times){
    ids = sample(1:nsize,nsize,replace=T)
    datasub.boot <- data[ids,]
    est_model_boot <- sem(model, data = datasub.boot)
    try(cv.out.boot <- cv_regsem(est_model_boot,...))
    try(lam[i]<-cv.out.boot$fit[which.min(cv.out.boot$fits[cv.out.boot$fit[,2]==0,4])])
  }
  if (sum(lam!=0,na.rm=T)==0){
    warning("0 penalty is selected by all runs.")
    lb=0;ub=0;
  }else{
    lb<-min(lam[lam!=0],na.rm=T)
    ub<-max(lam[lam!=0],na.rm=T)
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
