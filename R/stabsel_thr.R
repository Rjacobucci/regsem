#'
#'
#' Tuning the probability threshold.
#'
#' This function tune the probability threshold parameter.
#' @param stabsel output object from stabsel function. If specified, data, model, est_model, prob, nm, and pars.pen parameters are not needed.
#' @param data data frame
#' @param model lavaan syntax model.
#' @param est_model lavaan output object.
#' @param prob matrix of selection probabilities.
#' @param nm names(regsemOutput$coefficients).
#' @param pars.pen a vector of numbers corresponding to paths to be removed (same sequence as regsemOutput$coefficients).
#' @param from starting value of the threshold parameter.
#' @param to end value of the threshold parameter.
#' @param jump increment of the threshold parameter.
#' @param method fit indices uesd to tune the parameter.
#' @return rtn results using the optimal threshold.
#' @export


stabsel_thr<-function(stabsel=NULL,
                      data=NULL,
                      model=NULL,
                      est_model=NULL,
                      prob=NULL,
                      nm=NULL,
                      pars.pen=NULL,
                      from=0.5,
                      to=1,
                      jump=0.01,
                      method="aic"){
  if (!is.null(stabsel)){
    data<-stabsel$data
    model<-stabsel$model
    est_model<-stabsel$sem_model
    prob<-stabsel$probabilities
    nm<-stabsel$nm
    pars.pen<-stabsel$pars_pen
  }
  if (is.null(est_model)){
    est_model<-sem(data=data, model=model)
  }
  thr<-seq(from=from,to=to,by=jump)
  #library(matrixStats)
  max.prob<-matrixStats::colMaxs(as.matrix(prob))
  fit<-list()
  for (i in 1:length(thr)){
    p<-thr[i]
    sel.res<-max.prob>=p#true false
    rm.path<-pars.pen[sel.res==0]
    #coefficients from relaxed lasso:
    new.mod<-pen_mod(est_model,nm=nm,rm.path)
    new_est_mod<-sem(new.mod,data)
    fit[[i]]<-fitMeasures(new_est_mod)
  }
  m<-unlist(lapply(fit, '[[', which(names(fitMeasures(new_est_mod))==method)))
  if (method %in% c("aic","bic","bic2","rmsea","rmr","srmr")){
    id<-max(which(m==min(m)))
  }else if (method %in% c("cfi","nfi","nnfi","gfi","tli","rni","ifi","mfi","pgfi","agfi")){
    id<-max(which(m==max(m)))
  }else{
    warning(paste0(method," is not in {aic,bic,bic2,rmsea,rmr,srmr,cfi,nfi,nnfi,gfi,tli,rni,ifi,mfi,pgfi,agfi}. ",'Best result is chosen based on smallest value of ',method, ", which may be incorrect."))
    id<-max(which(m==min(m)))
  }

  opt.p<-thr[id]
  opt.sel.res<-max.prob>=opt.p#true false
  opt.rm.path<-pars.pen[opt.sel.res==0]
  opt.rm.path.nm<-colnames(prob)[opt.sel.res==0]#names of paths to be removed

  #coefficients from relaxed lasso:
  opt.new.mod<-pen_mod(est_model,nm,opt.rm.path)
  opt.new_est_mod<-sem(opt.new.mod,data)
  coefficient<-coef(opt.new_est_mod)

  #return:
  rtn<-list()
  rtn$data<-data
  rtn$model<-model
  rtn$sem_model<-est_model
  rtn$test_thresholds<-thr
  rtn$method<-method
  rtn$fit<-m
  rtn$opt_threshold<-opt.p
  rtn$opt_sel_results<-opt.sel.res
  rtn$remove_path<-opt.rm.path.nm
  rtn$opt_model<-opt.new.mod
  rtn$opt_model_est<-opt.new_est_mod
  rtn$coefficients<-coefficient
  rtn
}
