#'
#'
#' The main function that ties together and runs the models.
#' @param model lavaan output object.
#' @param n.lambda number of penalization values to test.
#' @param mult.start Logical. Whether to use multi_optim() (TRUE) or
#'         regsem() (FALSE).
#' @param multi.iter maximum number of random starts for multi_optim
#' @param jump Amount to increase penalization each iteration.
#' @param type penalty type.
#' @param fit.ret Fit indices to return.
#' @param fit.ret2 Return fits using only dataset "train" or bootstrap "boot"? Have to
#'        do 2 sample CV manually.
#' @param data Optional dataframe. Only required for missing="fiml".
#' @param optMethod solver to use.
#' @param gradFun gradient function to use.
#' @param hessFun hessian function to use.
#' @param test.cov Covariance matrix from test dataset. Necessary for CV=T
#' @param parallel Logical. whether to parallelize the processes running models for all
#'        values of lambda.
#' @param ncore Number of cores to use when parallel=TRUE
#' @param Start type of starting values to use.
#' @param subOpt type of optimization to use in the optimx package.
#' @param longMod longitudinal model?
#' @param optNL type of optimization to use in the NLopt package.
#' @param fac.type using cfa or efa type of model.
#' @param matrices function to use for extracting RAM matrices.
#' @param pars_pen parameter indicators to penalize.
#' @param diff_par parameter values to deviate from.
#' @param LB lower bound vector.
#' @param UB upper bound vector
#' @param calc type of calc function to use with means or not.
#' @param nlminb.control list of control values to pass to nlminb
#' @param warm.start Whether start values are based on previous iteration.
#'        This is not recommended.
#' @param missing How to handle missing data. Current options are "listwise"
#'        and "fiml".
#' @param ... Any additional arguments to pass to regsem() or multi_optim().
#' @keywords optim calc
#' @export
#' @examples
#' \dontrun{
#' library(lavaan)
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod,HS)
#'
#' cv.out = cv_regsem(outt,type="ridge",gradFun="none",n.lambda=100)
#'}



cv_regsem = function(model,
                     n.lambda=100,
                     mult.start=TRUE,
                     multi.iter=100,
                     jump=0.002,
                     type="none",
                     fit.ret=c("rmsea","BIC"),
                     fit.ret2 = "train",
                     data=NULL,
                     optMethod="nlminb",
                    gradFun="ram",
                    hessFun="none",
                    test.cov=NULL,
                    parallel=FALSE,
                    ncore=2,
                    Start="default",
                    subOpt="nlminb",
                    longMod=F,
                    optNL="NLOPT_LN_NEWUOA_BOUND",
                    fac.type="cfa",
                    matrices="extractMatrices",
                    pars_pen=NULL,
                    diff_par=NULL,
                    LB=-Inf,
                    UB=Inf,
                    calc="normal",
                    nlminb.control=list(),
                    warm.start=FALSE,
                    missing="listwise",
                    ...){


#if(fit.ret2 == "test"){
#  ids <-  sample(nrow(dat),nrow(dat)/2)
#  dat.train <- dat[ids,]
#  dat.test <- dat[-ids,]
#}



if(parallel==FALSE){
par.matrix <- matrix(0,n.lambda,model@Fit@npar)
fits <- matrix(NA,n.lambda,length(fit.ret)+2)
SHRINK = 0
count = 0
counts=n.lambda
#res2 <- data.frame(matrix(NA,counts,3))
#coefs = rep(1,14)

while(count < counts){

  count = count + 1
  print(count)
  SHRINK <- jump*(count-1) # 0.01 works well & 0.007 as well with 150 iterations

if(mult.start==FALSE){
  out <- regsem(model=model,lambda=SHRINK,type=type,data=data,
                   optMethod=optMethod,
                   gradFun=gradFun,hessFun=hessFun,
                   parallel=parallel,Start=Start,
                   subOpt=subOpt,
                   longMod=longMod,
                   optNL=optNL,
                   fac.type=fac.type,
                   matrices=matrices,
                   pars_pen=pars_pen,
                   diff_par=diff_par,
                   LB=LB,
                   UB=UB,
                   calc=calc,
                   nlminb.control=nlminb.control,
                   missing=missing)


  }else if(mult.start==TRUE){
   out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                      LB=LB,UB=UB,type=type,optMethod=optMethod,
                      gradFun=gradFun,hessFun=hessFun,nlminb.control=nlminb.control,
                      pars_pen=pars_pen,diff_par=NULL,warm.start=warm.start)
  }


  #if(any(fit.ret2 == "test")==TRUE){
  #  fits[[count]]$test = NA #fit_indices(out,CV=TRUE)[fit.ret]
  #}else
  if(fit.ret2 == "train"){
    fitt = try(fit_indices(out,CV=FALSE)$fits[fit.ret],silent=T)
    if(inherits(fitt, "try-error")) {
      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{
      fits[count,3:ncol(fits)] = fitt
    }

  }else if(fit.ret2 == "test"){
   # stop("fit.ret2=test is currently not implemented")
    fitt = try(fit_indices(out,CovMat=test.cov,CV=TRUE)$fits[fit.ret],silent=T)
    if(inherits(fitt, "try-error")) {
      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{
      fits[count,3:ncol(fits)] = fitt
    }
  }else if(fit.ret2 == "boot"){
    fitt = try(fit_indices(out,CV="boot")$fits[fit.ret],silent=T)
    if(inherits(fitt, "try-error")) {
      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{
      fits[count,3:ncol(fits)] = fitt
    }
  }
  fits[count,1] <- SHRINK
  fits[count,2] <- out$out$convergence

  if(is.null(out$coefficients)==TRUE){
    break
  }
  par.matrix[count,] = as.matrix(out$coefficients)

  colnames(par.matrix) = names(out$coefficients)
  colnames(fits) <- c("lambda","conv",fit.ret)
  out <- list(par.matrix,fits)
 # ret

}
}else if(parallel==TRUE){



  par.matrix <- matrix(0,n.lambda,model@Fit@npar)
  fits <- matrix(NA,n.lambda,length(fit.ret)+2)
  SHRINK = 0
  count = 0
  counts=n.lambda
  #res2 <- data.frame(matrix(NA,counts,3))
  #coefs = rep(1,14)

  #library(snowfall)

  cv_parallel <- function(SHRINK){

    if(mult.start==FALSE){
      out <- regsem(model=model,lambda=SHRINK,type=type,data=data,
                    optMethod=optMethod,
                    gradFun=gradFun,hessFun=hessFun,
                    parallel=parallel,Start=Start,
                    subOpt=subOpt,
                    longMod=longMod,
                    optNL=optNL,
                    fac.type=fac.type,
                    matrices=matrices,
                    pars_pen=pars_pen,
                    diff_par=diff_par,
                    LB=LB,
                    UB=UB,
                    calc=calc,
                    nlminb.control=nlminb.control,
                    missing=missing)


    }else if(mult.start==TRUE){
      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                         LB=LB,UB=UB,type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,nlminb.control=nlminb.control,
                         pars_pen=pars_pen,diff_par=NULL,warm.start=warm.start)
    }


    #if(any(fit.ret2 == "test")==TRUE){
    #  fits[[count]]$test = NA #fit_indices(out,CV=TRUE)[fit.ret]
    #}else
    if(fit.ret2 == "train"){
      fitt = try(fit_indices(out,CV=FALSE)$fits[fit.ret],silent=T)
      if(inherits(fitt, "try-error")) {
        fitss = rep(NA,ncol(fits)-2)
      }else{
        fitss = fitt
      }

    }else if(fit.ret2 == "test"){
      # stop("fit.ret2=test is currently not implemented")
      fitt = try(fit_indices(out,CovMat=test.cov,CV=TRUE)$fits[fit.ret],silent=T)
      if(inherits(fitt, "try-error")) {
        fitss = rep(NA,ncol(fits)-2)
      }else{
        fitss = fitt
      }
    }else if(fit.ret2 == "boot"){
      fitt = try(fit_indices(out,CV="boot")$fits[fit.ret],silent=T)
      if(inherits(fitt, "try-error")) {
        fitss = rep(NA,ncol(fits)-2)
      }else{
        fitss = fitt
      }
    }
    fitss <- matrix(fitss,1,length(fit.ret))
    return(data.frame(SHRINK,conv=out$out$convergence,fitss,out$coefficients))
  }



  snowfall::sfLibrary(regsem)
  snowfall::sfInit(parallel=TRUE, cpus=ncore)
  snowfall::sfExport("model","type","data",
                     "optMethod",
                     "gradFun","hessFun",
                     "parallel","Start",
                     "subOpt",
                     "longMod",
                     "optNL",
                     "fac.type",
                     "matrices",
                     "pars_pen",
                     "diff_par",
                     "LB",
                     "UB",
                     "calc",
                     "nlminb.control",
                     "warm.start",
                     "missing")




  lambdas <- seq(0,by=jump,length.out=n.lambda)
  ret = sfLapply(lambdas,cv_parallel)
  snowfall::sfStop()

  #out

  out <- unlist(ret)
  out <- matrix(out,nrow=n.lambda,ncol=length(ret[[1]]),byrow=T)
  nam <- names(extractMatrices(model)$parameters)
  colnames(out) <- c("lambda","conv",fit.ret,nam)
  out



}
#fits = fit_indices(out,CV=FALSE)

out

}
