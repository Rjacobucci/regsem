#'
#'
#' Stability selection
#' @param data data frame
#' @param model lavaan syntax model.
#' @param det.range Whether to determine the range of penalization values for stability selection through bootstrapping. Default is FALSE, from and to arguments are then needed. If set to TRUE, then jump, times and detr.nlambda arguments will be needed.
#' @param from Minimum value of penalization values for stability selection.
#' @param to Maximum value of penalization values for stability selection.
#' @param times Number of bootstrapping sample used to determine the range. Default is 50.
#' @param jump Amount to increase penalization each iteration. Default is 0.01
#' @param detr.nlambda Number of penalization values to test for determining range.
#' @param n.lambda Number of penalization values to test for stability selection.
#' @param n.boot Number of bootstrap samples needed for stability selection.
#' @param det.thr Whether to determine the probability threshold value. Default is FALSE, p is then needed. If set to TRUE, p.from, p.to, p.method arguments will be needed.
#' @param p Probability threshold: above which selection probability is the path kept in the modle. Default value is 0.8.
#' @param p.from Lower bound of probability threshold to test. Default is 0.5.
#' @param p.to Upper bound of probability threshold to test. Default is 1.
#' @param p.jump Amount to increase threshold each iteration. Default is 0.05.
#' @param p.method Which fit index to use to choose a final model?
#' @param type Penalty type
#' @param pars_pen Parameter indicators to penalize.
#' @param ... Any additional arguments to pass to regsem() or cv_regsem().
#' @examples
#' \dontrun{
#' library(regsem)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ 1*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' x1 ~~ r1*x2;x1 ~~ r2*x3;x1 ~~ r3*x4;x1 ~~ r4*x5
#' '
#' outt = cfa(mod, HS)
#'
#'stabsel.out = stabsel(data=HS,model=mod,det.range=T,detr.nlambda=20,n.lambda=5,
#'                     n.boot=10,p=0.9,type="alasso", p.method="aic",
#'                     pars_pen=c("r1","r2","r3","r4"))
#' stabsel.out$selection_results
#'
#' }
#' @export


stabsel<-function(data,
                  model,
                  det.range=FALSE,
                  from,
                  to,
                  times=50,
                  jump=0.01,
                  detr.nlambda=20,
                  n.lambda=40,
                  n.boot=100,
                  det.thr=FALSE,
                  p=0.8,
                  p.from=0.5,
                  p.to=1,
                  p.jump=0.05,
                  p.method="aic",
                  type="lasso",
                  pars_pen="regressions",
                  ...){
  rtn<-list()
  #determine range for lambda (by prestage or user specify)
  if (det.range==TRUE){#use a prestep to determine lambda range
    lam<-det_range(data,model,times,n.lambda=detr.nlambda,jump=jump,type=type,pars_pen=pars_pen,...)
    lb<-lam$lb;ub<-lam$ub
    #rtn$det.range<-lam
  }else{
    lb<-from;ub<-to
  }
  #determine lambda values:
  Lam<-seq(from=lb,to=ub,length.out=n.lambda)
  if (ub==0){
    return("determine range fails")
  }else if (lb==0){
    warning("removed lambda = 0")
    Lam<-Lam[Lam!=0]
  }


  nlambda<-length(Lam)
  #
  est_model <- sem(model, data)
  parT<-parTable(est_model)
  regsem.out.ss <- regsem(est_model,lambda =1e-10,type=type,pars_pen=pars_pen,...)
  pars.pen<-regsem.out.ss$pars_pen
  nm<-names(regsem.out.ss$coefficients)
  n.pen<-length(pars.pen)
  #stability selection:
  nsize<-dim(data)[1]
  p.select<-data.frame(matrix(NA,1,n.pen))#probabilities of being selected
  for (k in 1:nlambda){
    for(i in 1:n.boot){
      ids = sample(1:nsize,nsize,replace=T)
      datasub.ss <- data[ids,]
      data.hold.ss<-data[-unique(ids),]

      try(est_model_ss <- sem(model, data = datasub.ss))
      try(regsem.out.ss <- regsem(est_model_ss,lambda = Lam[k],type=type,pars_pen=pars_pen,...))
      coeff<-rep(NA,n.pen)
      try(coeff<-regsem.out.ss$coefficients[pars.pen])
      if (i==1){
        pars = coeff
      }else{
        pars = rbind(pars,coeff)
      }
    }
    pars.ss.1<-pars!=0
    p.select[k,]<-colSums(pars.ss.1,na.rm=T)/dim(na.omit(pars.ss.1))[1]
  }
  names(p.select)<-nm[pars.pen]

  #Selection results:

  if (det.thr==FALSE){
    #library(matrixStats)
    max.prob<-matrixStats::colMaxs(as.matrix(p.select))
    sel.res<-max.prob>=p#true false
    names(sel.res)<-colnames(pars.ss.1)
    sel.nm<-colnames(p.select)[sel.res]#names of selected paths
    sel.pars<-pars.pen[sel.res]#label of selected paths in number
    rm.path.nm<-colnames(p.select)[sel.res==0]#names of paths to be removed
    rm.path<-pars.pen[sel.res==0]

    #coefficients from relaxed lasso:
    new.mod<-pen_mod(est_model,nm,rm.path)
    new_est_mod<-sem(new.mod,data)
    coefficient<-coef(new_est_mod)

    #returns:
    rtn$data<-data
    rtn$model<-model
    rtn$sem_model<-est_model
    rtn$nm<-nm
    rtn$pars_pen<-pars.pen
    rtn$Lambda<-Lam
    rtn$probabilities<-p.select
    rtn$threshold<-p
    rtn$selection_results<-sel.res
    rtn$remove_path<-rm.path.nm
    rtn$new_model<-new.mod
    rtn$new_model_est<-new_est_mod
    rtn$coefficients<-coefficient
    rtn
  }else {
    thr<-stabsel_thr(data=data,model=model,est_model=est_model,prob=p.select,nm=nm,pars.pen=pars.pen,from=p.from,to=p.to,jump=p.jump,method=p.method)

    #returns:
    rtn$data<-data
    rtn$model<-model
    rtn$sem_model<-est_model
    rtn$nm<-names(regsem.out.ss$coefficients)
    rtn$pars_pen<-pars.pen
    rtn$Lambda<-Lam
    rtn$probabilities<-p.select
    rtn$test_thresholds<-thr$test_thresholds
    rtn$p.method<-thr$method
    rtn$fit<-thr$fit
    rtn$opt_threshold<-thr$opt_threshold
    rtn$opt_sel_results<-thr$opt_sel_results
    rtn$remove_path<-thr$remove_path
    rtn$new_model<-thr$opt_model
    rtn$new_model_est<-thr$opt_model_est
    rtn$coefficients<-thr$coefficients
    rtn
  }
}
