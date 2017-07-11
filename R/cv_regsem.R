#'
#'
#' The main function that ties together and runs the models.
#' @param model lavaan output object.
#' @param n.lambda number of penalization values to test.
#' @param pars_pen parameter indicators to penalize.
#' @param mult.start Logical. Whether to use multi_optim() (TRUE) or
#'         regsem() (FALSE).
#' @param multi.iter maximum number of random starts for multi_optim
#' @param jump Amount to increase penalization each iteration.
#' @param lambda.start What value to start the penalty at
#' @param alpha Mixing for elastic net
#' @param type Penalty type. Options include "none", "lasso", "ridge",
#'        "enet" for the elastic net,
#'        "alasso" for the adaptive lasso
#'        and "diff_lasso". diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par. Two methods for
#'        sparser results than lasso are the smooth clipped absolute deviation,
#'        "scad", and the minimum concave penalty, "mcp".
#' @param fit.ret Fit indices to return.
#' @param fit.ret2 Return fits using only dataset "train" or bootstrap "boot"? Have to
#'        do 2 sample CV manually.
#' @param n.boot Number of bootstrap samples if fit.ret2="boot"
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
#' @param diff_par parameter values to deviate from.
#' @param LB lower bound vector.
#' @param UB upper bound vector
#' @param par.lim Vector of minimum and maximum parameter estimates. Used to
#'        stop optimization and move to new starting values if violated.
#' @param block Whether to use block coordinate descent
#' @param full Whether to do full gradient descent or block
#' @param calc Type of calc function to use with means or not. Not recommended
#'        for use.
#' @param nlminb.control list of control values to pass to nlminb
#' @param max.iter Number of iterations for coordinate descent
#' @param tol Tolerance for coordinate descent
#' @param solver Whether to use solver for coord_desc
#' @param quasi Whether to use quasi-Newton
#' @param solver.maxit Max iterations for solver in coord_desc
#' @param alpha.inc Whether alpha should increase for coord_desc
#' @param step Step size
#' @param momentum Momentum for step sizes
#' @param step.ratio Ratio of step size between A and S. Logical
#' @param line.search Use line search for optimization. Default is no, use fixed step size
#' @param warm.start Whether start values are based on previous iteration.
#'        This is not recommended.
#' @param missing How to handle missing data. Current options are "listwise"
#'        and "fiml".
#' @param ... Any additional arguments to pass to regsem() or multi_optim().
#' @keywords optim calc
#' @export
#' @examples
#' \dontrun{
#' library(regsem)
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod,HS)
#'
#' cv.out = cv_regsem(outt,type="ridge",pars_pen=c(1:2,6:8),n.lambda=100)
#' # check parameter numbers
#' extractMatrices(outt)["A"]
#' # equivalent to
#' mod <- '
#' f =~ 1*x1 + l1*x2 + l2*x3 + l3*x4 + l4*x5 + l5*x6 + l6*x7 + l7*x8 + l8*x9
#' '
#' outt = cfa(mod,HS)
#'
#' cv.out = cv_regsem(outt,type="ridge",pars_pen=c("l1","l2","l6","l7","l8"),
#'          n.lambda=100)
#' plot(cv.out,show.minimum="BIC")
#' }



cv_regsem = function(model,
                     n.lambda=100,
                     pars_pen,
                     mult.start=FALSE,
                     multi.iter=10,
                     jump=0.002,
                     lambda.start=0,
                     alpha=.5,
                     type="lasso",
                     fit.ret=c("rmsea","BIC"),
                     fit.ret2 = "train",
                     n.boot=20,
                     data=NULL,
                     optMethod="coord_desc",
                    gradFun="ram",
                    hessFun="none",
                    test.cov=NULL,
                    parallel=FALSE,
                    ncore=2,
                    Start="lavaan",
                    subOpt="nlminb",
                    diff_par=NULL,
                    LB=-Inf,
                    UB=Inf,
                    par.lim=c(-Inf,Inf),
                    block=TRUE,
                    full=TRUE,
                    calc="normal",
                    max.iter=2000,
                    tol=1e-5,
                    solver=FALSE,
                    quasi=FALSE,
                    solver.maxit=5,
                    alpha.inc=FALSE,
                    step=.1,
                    momentum=FALSE,
                    step.ratio=FALSE,
                    line.search=FALSE,
                    nlminb.control=list(),
                    warm.start=TRUE,
                    missing="listwise",
                    ...){


#if(fit.ret2 == "test"){
#  ids <-  sample(nrow(dat),nrow(dat)/2)
#  dat.train <- dat[ids,]
#  dat.test <- dat[-ids,]
#}
fits.var=NA

if(is.null(pars_pen) & type!="none"){
  stop("for cv_regsem(), pars_pen needs to be specified")
}

if(is.null(pars_pen)==FALSE & is.numeric(pars_pen)==FALSE){
  pars_pen <- parse_parameters(pars_pen,model)
}


if(parallel == TRUE){
  stop("parallel is not currently supported")
}

if(parallel==FALSE){
par.matrix <- matrix(0,n.lambda,length(extractMatrices(model)$parameters))
fits <- matrix(NA,n.lambda,length(fit.ret)+2)
fit.reg <- rep(NA,n.lambda)
fitt.var <- matrix(NA,n.lambda,length(fit.ret))
SHRINK2 = lambda.start
count = 0
counts=n.lambda
#res2 <- data.frame(matrix(NA,counts,3))
#coefs = rep(1,14)

while(count < counts){
  count = count + 1
  print(count)
  SHRINK <- SHRINK2 + jump*(count-1) # 0.01 works well & 0.007 as well with 150 iterations

  if(count > 1 & all(abs(par.matrix[count-1,pars_pen])<.001)){
    break
  }

if(mult.start==FALSE){

  if(warm.start==FALSE | count == 1){
    itt = 0
    Start="lavaan"
  }else if(fits[count-1,2] == 0){
    itt = 0
    Start = par.matrix[count-1,]
    Start[pars_pen] = Start[pars_pen]-step*jump
  }else if(fits[count-1,2] == 99){
    Start="lavaan"
  }else{
    itt = itt + 1
    Start = par.matrix[count-itt-1,]
    Start[pars_pen] = Start[pars_pen]-itt*jump
  }


  if(fit.ret2 == "train"){
    out <- regsem(model=model,lambda=SHRINK,type=type,data=data,
                  optMethod=optMethod,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  solver=solver,
                  quasi=quasi,
                  solver.maxit=solver.maxit,
                  alpha.inc=alpha.inc,
                  step=step,
                  max.iter=max.iter,
                  line.search=line.search,
                  momentum=momentum,
                  step.ratio=step.ratio,
                  nlminb.control=nlminb.control,
                  missing=missing)
  }else if(fit.ret2=="boot"){

    fitt <- matrix(NA,n.boot,length(fit.ret))

    out <- regsem(model=model,lambda=SHRINK,type=type,data=NULL,
                  optMethod=optMethod,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  solver=solver,
                  quasi=quasi,
                  solver.maxit=solver.maxit,
                  alpha.inc=alpha.inc,
                  step=step,
                  max.iter=max.iter,
                  line.search=line.search,
                  momentum=momentum,
                  step.ratio=step.ratio,
                  nlminb.control=nlminb.control,
                  missing=missing)

    for(i in 1:n.boot){
      set.seed(i)
      data <- as.data.frame(model@Data@X)

      ids1 <- sample(1:nrow(data),nrow(data),replace=TRUE)

      train <- data[ids1,]
      test <- data[-ids1,]

      colnames(train) <- model@pta$vnames$ov[[1]]
      colnames(test) <- model@pta$vnames$ov[[1]]


      mod1 <- lavaan(parTable(model),train)

    out2 <- regsem(model=mod1,lambda=SHRINK,type=type,data=NULL,
                  optMethod=optMethod,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  solver=solver,
                  quasi=quasi,
                  solver.maxit=solver.maxit,
                  alpha.inc=alpha.inc,
                  step=step,
                  max.iter=max.iter,
                  line.search=line.search,
                  momentum=momentum,
                  step.ratio=step.ratio,
                  nlminb.control=nlminb.control,
                  missing=missing)

    if(out$convergence==0){
      fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=cov(test))$fits[fit.ret]

    }else{
      fitt[i,] = NA
    }
    }
    fits[count,3:ncol(fits)] <- apply(fitt, 2, function(x) mean(x, trim = .2,na.rm=TRUE))
    fitt.var[count,1:length(fit.ret)] <- apply(fitt, 2, function(x) var(x,na.rm=TRUE))
  }else if(fit.ret2=="cv"){


    fitt <- matrix(NA,5,length(fit.ret))

    out <- regsem(model=model,lambda=SHRINK,type=type,data=NULL,
                  optMethod=optMethod,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  solver=solver,
                  quasi=quasi,
                  solver.maxit=solver.maxit,
                  alpha.inc=alpha.inc,
                  step=step,
                  max.iter=max.iter,
                  line.search=line.search,
                  momentum=momentum,
                  step.ratio=step.ratio,
                  nlminb.control=nlminb.control,
                  missing=missing)


    for(i in 1:5){
      set.seed(i)
      data <- as.data.frame(model@Data@X)

      ids1 <- sample(1:nrow(data),round(nrow(data)*.80,1),replace=FALSE)

      train <- data[ids1,]
      test <- data[-ids1,]

      colnames(train) <- model@pta$vnames$ov[[1]]
      colnames(test) <- model@pta$vnames$ov[[1]]


      mod1 <- lavaan(parTable(model),train)

      out2 <- regsem(model=mod1,lambda=SHRINK,type=type,data=NULL,
                     optMethod=optMethod,
                     gradFun=gradFun,hessFun=hessFun,
                     parallel=parallel,Start=Start,
                     subOpt=subOpt,
                     alpha=alpha,
                     pars_pen=pars_pen,
                     diff_par=diff_par,
                     LB=LB,
                     UB=UB,
                     par.lim=par.lim,
                     block=block,
                     full=full,
                     calc=calc,
                     tol=tol,
                     solver=solver,
                     quasi=quasi,
                     solver.maxit=solver.maxit,
                     alpha.inc=alpha.inc,
                     step=step,
                     line.search=line.search,
                     max.iter=max.iter,
                     momentum=momentum,
                     step.ratio=step.ratio,
                     nlminb.control=nlminb.control,
                     missing=missing)

      if(out$convergence==0){
        fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=cov(test))$fits[fit.ret]
      }else{
        fitt[i,] = NA
      }
    }
    fits[count,3:ncol(fits)] <- apply(fitt, 2, function(x) mean(x, trim = .2,na.rm=TRUE))
    fitt.var[count,1:length(fit.ret)] <- apply(fitt, 2, function(x) var(x,na.rm=TRUE))
  }




  }else if(mult.start==TRUE){

    if(warm.start==FALSE | count == 1 | count == 99){
      itt = 0
      Start2=NULL
    }else if(fits[count-1,2] == 0){
      itt = 0
      Start2 = par.matrix[count-1,]
      Start2[pars_pen] = Start2[pars_pen]-step*jump
    }else if(fits[count-1,2] == 99){
      Start="lavaan"
    }else{
      itt = itt + 1
      Start2 = par.matrix[count-itt-1,]
      Start2[pars_pen] = Start2[pars_pen]-itt*jump
    }




    if(fit.ret2 == "train"){
      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                      LB=LB,UB=UB,par.lim=par.lim,
                      type=type,optMethod=optMethod,
                      gradFun=gradFun,hessFun=hessFun,
                      tol=tol,
                      alpha=alpha,
                      solver=solver,
                      quasi=quasi,
                      solver.maxit=solver.maxit,
                      max.iter=max.iter,
                      full=full,
                      block=block,
                      alpha.inc=alpha.inc,
                      line.search=line.search,
                      step=step,
                      momentum=momentum,
                      Start2=Start2,
                      step.ratio=step.ratio,nlminb.control=nlminb.control,
                      pars_pen=pars_pen,diff_par=NULL)

    }else if(fit.ret2=="boot"){
      fitt <- matrix(NA,n.boot,length(fit.ret))

      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                         LB=LB,UB=UB,par.lim=par.lim,
                         type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,
                         tol=tol,
                         alpha=alpha,
                         solver=solver,
                         quasi=quasi,
                         solver.maxit=solver.maxit,
                         max.iter=max.iter,
                         full=full,
                         block=block,
                         alpha.inc=alpha.inc,
                         line.search=line.search,
                         step=step,
                         momentum=momentum,
                         Start2=Start2,
                         step.ratio=step.ratio,nlminb.control=nlminb.control,
                         pars_pen=pars_pen,diff_par=NULL)


      for(i in 1:n.boot){
        set.seed(i)
        data <- as.data.frame(model@Data@X)

        ids1 <- sample(1:nrow(data),nrow(data),replace=TRUE)

        train <- data[ids1,]
        test <- data[-ids1,]

        colnames(train) <- model@pta$vnames$ov[[1]]
        colnames(test) <- model@pta$vnames$ov[[1]]

        mod1 <- lavaan(parTable(model),train)

        out2 <- multi_optim(model=mod1,max.try=multi.iter,lambda=SHRINK,
                           LB=LB,UB=UB,par.lim=par.lim,
                           type=type,optMethod=optMethod,
                           gradFun=gradFun,hessFun=hessFun,
                           tol=tol,
                           alpha=alpha,
                           solver=solver,
                           quasi=quasi,
                           solver.maxit=solver.maxit,
                           max.iter=max.iter,
                           full=full,
                           block=block,
                           alpha.inc=alpha.inc,
                           line.search=line.search,
                           step=step,
                           momentum=momentum,
                           Start2=Start2,
                           step.ratio=step.ratio,nlminb.control=nlminb.control,
                           pars_pen=pars_pen,diff_par=NULL)


        if(out$convergence==0){
          fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=cov(test))$fits[fit.ret]
        }else{
          fitt[i,] = NA
        }
      }
      fits[count,3:ncol(fits)] <- apply(fitt, 2, function(x) mean(x, trim = .2,na.rm=TRUE))
      fitt.var[count,1:length(fit.ret)] <- apply(fitt, 2, function(x) var(x,na.rm=TRUE))
    }else if(fit.ret2=="cv"){


      fitt <- matrix(NA,5,length(fit.ret))

      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                         LB=LB,UB=UB,par.lim=par.lim,
                         type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,
                         tol=tol,
                         alpha=alpha,
                         solver=solver,
                         quasi=quasi,
                         solver.maxit=solver.maxit,
                         max.iter=max.iter,
                         full=full,
                         block=block,
                         alpha.inc=alpha.inc,
                         step=step,
                         line.search=line.search,
                         momentum=momentum,
                         Start2=Start2,
                         step.ratio=step.ratio,nlminb.control=nlminb.control,
                         pars_pen=pars_pen,diff_par=NULL)


      for(i in 1:5){
        set.seed(i)
        data <- as.data.frame(model@Data@X)

        ids1 <- sample(1:nrow(data),round(nrow(data)*.80,1),replace=FALSE)

        train <- data[ids1,]
        test <- data[-ids1,]

        colnames(train) <- model@pta$vnames$ov[[1]]
        colnames(test) <- model@pta$vnames$ov[[1]]


        mod1 <- lavaan(parTable(model),train)

        out2 <- multi_optim(model=mod1,max.try=multi.iter,lambda=SHRINK,
                            LB=LB,UB=UB,par.lim=par.lim,
                            type=type,optMethod=optMethod,
                            gradFun=gradFun,hessFun=hessFun,
                            tol=tol,
                            alpha=alpha,
                            solver=solver,
                            quasi=quasi,
                            solver.maxit=solver.maxit,
                            max.iter=max.iter,
                            full=full,
                            block=block,
                            alpha.inc=alpha.inc,
                            step=step,
                            line.search=line.search,
                            momentum=momentum,
                            Start2=Start2,
                            step.ratio=step.ratio,nlminb.control=nlminb.control,
                            pars_pen=pars_pen,diff_par=NULL)

        if(out$convergence==0){
          fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=cov(test))$fits[fit.ret]
        }else{
          fitt[i,] = NA
        }
      }
      fits[count,3:ncol(fits)] <- apply(fitt, 2, function(x) mean(x, trim = .2,na.rm=TRUE))
      fitt.var[count,1:length(fit.ret)] <- apply(fitt, 2, function(x) var(x,na.rm=TRUE))
    }

  }
  #print(pars_pen)
 # pars_pen <- out$pars_pen
  #if(any(fit.ret2 == "test")==TRUE){
  #  fits[[count]]$test = NA #fit_indices(out,CV=TRUE)[fit.ret]
  #}else
  if(fit.ret2 == "train"){
    fitt = try(fit_indices(out,CV=FALSE)$fits[fit.ret],silent=T)
    fit.reg[count] <- out$optim_fit
    if(inherits(fitt, "try-error")) {
      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{
      fits[count,3:ncol(fits)] = fitt
    }

  }else if(fit.ret2 == "test"){
   # stop("fit.ret2=test is currently not implemented")
    #print(summary(out))

    fitt = try(fit_indices(out,CovMat=test.cov,CV=TRUE)$fits[fit.ret],silent=T)
    if(inherits(fitt, "try-error")) {

      fits[count,3:ncol(fits)] = rep(NA,ncol(fits)-2)
    }else{

      fits[count,3:ncol(fits)] = fitt
    }
  }
  fits[count,1] <- SHRINK

#  if(class(out$convergence)=="numeric"){
  #print(class(out$convergence));print(1)
  #print(out$convergence);print(class(out$convergence))
    fits[count,2] <- out$convergence

 # }else{
 #   fits[count,2] <- 99
    #out$convergence <- 99
 # }


  if(is.null(out$coefficients)==TRUE){
    break
  }
  par.matrix[count,] = as.matrix(out$coefficients)

  colnames(par.matrix) = names(out$coefficients)
  colnames(fits) <- c("lambda","conv",fit.ret)
  out2 <- list(par.matrix,fits,pars_pen,fitt.var,fit.reg)
 # ret

}
}else if(parallel==TRUE){

  stop("Parallel is not currently recommended")

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
                    pars_pen=pars_pen,
                    diff_par=diff_par,
                    LB=LB,
                    alpha=alpha,
                    UB=UB,
                    calc=calc,
                    nlminb.control=nlminb.control,
                    tol=tol,
                    full=full,
                    block=block,
                    solver=solver,
                    quasi=quasi,
                    solver.maxit=solver.maxit,
                    alpha.inc=alpha.inc,
                    line.search=line.search,
                    step=step,
                    momentum=momentum,
                    step.ratio=step.ratio,
                    missing=missing)


    }else if(mult.start==TRUE){
      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                         LB=LB,UB=UB,type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,nlminb.control=nlminb.control,
                         tol=tol,
                         full=full,
                         alpha=alpha,
                         block=block,
                         solver=solver,
                         quasi=quasi,
                         solver.maxit=solver.maxit,
                         alpha.inc=alpha.inc,
                         step=step,
                         line.search=line.search,
                         momentum=momentum,
                         step.ratio=step.ratio,
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
    data.frame(SHRINK,conv=out$convergence,fitss,out$coefficients)
  }



  snowfall::sfLibrary(regsem)
  snowfall::sfInit(parallel=TRUE, cpus=ncore)
  snowfall::sfExport("model","type","data",
                     "optMethod",
                     "gradFun","hessFun",
                     "parallel","Start",
                     "subOpt",
                     "pars_pen",
                     "diff_par",
                     "LB",
                     "block",
                     "solver",
                     "quasi",
                     "full",
                     "line.search",
                     "UB",
                     "calc",
                     "nlminb.control",
                     "warm.start",
                     "missing")




  lambdas <- seq(0,by=jump,length.out=n.lambda)
  ret = snowfall::sfLapply(lambdas,cv_parallel)
  snowfall::sfStop()

  #out
  pars_pen <- out$pars_pen

  out2 <- unlist(ret)
  out2 <- matrix(out,nrow=n.lambda,ncol=length(ret[[1]]),byrow=T)
  nam <- names(extractMatrices(model)$parameters)
  colnames(out2) <- c("lambda","conv",fit.ret,nam)
  out2



}
#fits = fit_indices(out,CV=FALSE)
out2$pars_pen <- pars_pen
out2$call <- match.call()
class(out2) <- "cvregsem"
out2

}
