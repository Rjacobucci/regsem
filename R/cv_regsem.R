#'
#'
#' The main function that runs multiple penalty values.
#'
#' @param model Lavaan output object. This is a model that was previously
#'        run with any of the lavaan main functions: cfa(), lavaan(), sem(),
#'        or growth(). It also can be from the efaUnrotate() function from
#'        the semTools package. Currently, the parts of the model which cannot
#'        be handled in regsem is the use of multiple group models, missing
#'        other than listwise, thresholds from categorical variable models,
#'        the use of additional estimators other than
#'        ML, most notably WLSMV for categorical variables. Note: the model
#'        does not have to actually run (use do.fit=FALSE), converge etc...
#'        regsem() uses the lavaan object as more of a parser and to get
#'        sample covariance matrix.
#' @param n.lambda number of penalization values to test.
#' @param pars_pen Parameter indicators to penalize. There are multiple ways to specify.
#'        The default is to penalize all regression parameters ("regressions"). Additionally,
#'        one can specify all loadings ("loadings"), or both c("regressions","loadings").
#'        Next, parameter labels can be assigned in the lavaan syntax and passed to pars_pen.
#'        See the example.Finally, one can take the parameter numbers from the A or S matrices and pass these
#'        directly. See extractMatrices(lav.object)$A.
#' @param metric Which fit index to use to choose a final model?
#'        Note that it chooses the best fit that also achieves convergence
#'        (conv=0).
#' @param mult.start Logical. Whether to use multi_optim() (TRUE) or
#'         regsem() (FALSE).
#' @param multi.iter maximum number of random starts for multi_optim
#' @param jump Amount to increase penalization each iteration.
#' @param lambda.start What value to start the penalty at
#' @param alpha Mixture for elastic net. 1 = ridge, 0 = lasso
#' @param gamma Additional penalty for MCP and SCAD
#' @param type Penalty type. Options include "none", "lasso", "ridge",
#'        "enet" for the elastic net,
#'        "alasso" for the adaptive lasso
#'        and "diff_lasso". diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par. Two methods for
#'        sparser results than lasso are the smooth clipped absolute deviation,
#'        "scad", and the minimum concave penalty, "mcp". Last option is "rlasso"
#'        which is the randomised lasso to be used for stability selection.
#' @param random.alpha Alpha parameter for randomised lasso. Has to be between
#'        0 and 1, with a default of 0.5. Note this is only used for
#'        "rlasso", which pairs with stability selection.
#' @param fit.ret Fit indices to return.
#' @param fit.ret2 Return fits using only dataset "train" or bootstrap "boot"? Have to
#'        do 2 sample CV manually.
#' @param n.boot Number of bootstrap samples if fit.ret2="boot"
#' @param data Optional dataframe. Only required for missing="fiml".
#' @param optMethod Solver to use. Two main options for use: rsoolnp and coord_desc.
#'        Although slightly slower, rsolnp works much better for complex models.
#'        coord_desc uses gradient descent with soft thresholding for the type of
#'        of penalty. Rsolnp is a nonlinear solver that doesn't rely on gradient
#'        information. There is a similar type of solver also available for use,
#'        slsqp from the nloptr package. coord_desc can also be used with hessian
#'        information, either through the use of quasi=TRUE, or specifying a hess_fun.
#'        However, this option is not recommended at this time.
#' @param gradFun Gradient function to use. Recommended to use "ram",
#'        which refers to the method specified in von Oertzen & Brick (2014).
#'        Only for use with optMethod="coord_desc".
#' @param hessFun hessian function to use. Currently not recommended.
#' @param test.cov Covariance matrix from test dataset. Necessary for CV=T
#' @param test.n.obs Number of observations in test set. Used when CV=T
#' @param prerun Logical. Use rsolnp to first optimize before passing to
#'        gradient descent? Only for use with coord_desc
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
#' @param round Number of digits to round results to
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
#' @param verbose Print progress bar?
#' @param ... Any additional arguments to pass to regsem() or multi_optim().
#' @keywords optim calc
#' @export
#' @examples
#' \dontrun{
#' library(regsem)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod, HS)
#' # increase to > 25
#' cv.out = cv_regsem(outt,type="lasso", pars_pen=c(1:2,6:8),
#'           n.lambda=5,jump=0.01)
#' # check parameter numbers
#' extractMatrices(outt)["A"]
#' # equivalent to
#' mod <- '
#' f =~ 1*x1 + l1*x2 + l2*x3 + l3*x4 + l4*x5 + l5*x6 + l6*x7 + l7*x8 + l8*x9
#' '
#' outt = cfa(mod,HS)
#' # increase to > 25
#' cv.out = cv_regsem(outt, type="lasso", pars_pen=c("l1","l2","l6","l7","l8"),
#'          n.lambda=5,jump=0.01)
#' summary(cv.out)
#' plot(cv.out, show.minimum="BIC")
#'
#' mod <- '
#'f =~ x1 + x2 + x3 + x4 + x5 + x6
#''
#'outt = cfa(mod, HS)
#'# can penalize all loadings
#'cv.out = cv_regsem(outt,type="lasso", pars_pen="loadings",
#'                   n.lambda=5,jump=0.01)
#'
#'mod2 <- '
#'f =~ x4+x5+x3
#'#x1 ~ x7 + x8 + x9 + x2
#'x1 ~ f
#'x2 ~ f
#''
#'outt2 = cfa(mod2, HS)
#'extractMatrices(outt2)$A
#' # if no pars_pen specification, defaults to all
#' # regressions
#'cv.out = cv_regsem(outt2,type="lasso",
#'                   n.lambda=15,jump=0.03)
#'# check
#'cv.out$pars_pen
#' }



cv_regsem = function(model,
                     n.lambda=40,
                     pars_pen="regressions",
                     metric=ifelse(fit.ret2=="train","BIC","chisq"),
                     mult.start=FALSE,
                     multi.iter=10,
                     jump=0.01,
                     lambda.start=0,
                     alpha=.5,
                     gamma=3.7,
                     type="lasso",
                     random.alpha=0.5,
                     fit.ret=c("rmsea","BIC","chisq"),
                     fit.ret2 = "train",
                     n.boot=20,
                     data=NULL,
                     optMethod="rsolnp",
                    gradFun="ram",
                    hessFun="none",
                    test.cov=NULL,
                    test.n.obs = NULL,
                    prerun=FALSE,
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
                    round=3,
                    solver=FALSE,
                    quasi=FALSE,
                    solver.maxit=5,
                    alpha.inc=FALSE,
                    step=.1,
                    momentum=FALSE,
                    step.ratio=FALSE,
                    line.search=FALSE,
                    nlminb.control=list(),
                    warm.start=FALSE,
                    missing="listwise",
                    verbose=TRUE,
                    ...){


#if(fit.ret2 == "test"){
#  ids <-  sample(nrow(dat),nrow(dat)/2)
#  dat.train <- dat[ids,]
#  dat.test <- dat[-ids,]
#}
fits.var=NA
mats <- extractMatrices(model)

if(metric %in% fit.ret == FALSE){
  stop("Need to specify metric= to match one index in fit.ret")
}

pars_pen2 = NULL

if(any(pars_pen=="regressions") & is.null(mats$regressions)){
  stop("No regression parameters to regularize")
}

if(any(pars_pen == "loadings")){
  pars_pen2 = mats$loadings
}else if(any(pars_pen == "regressions") | is.null(pars_pen)){
  pars_pen2 = c(pars_pen2,mats$regressions)
 # if(is.na(mats$name.factors)==TRUE){
 #   if(any(colnames(mats$A) == "1")){
  #    IntCol = which(colnames(mats$A) == "1")
  #    A_minusInt = mats$A[,-IntCol]
  #    A_pen = A_minusInt != 0
  #    pars_pen2 = c(A_minusInt[A_pen],pars_pen2)
  #  }else{
  #    A_pen = mats$A != 0
  #    pars_pen2 = c(mats$A[A_pen],pars_pen2)
  #  }
 # }else{
    # remove factor loadings
 #   if(any(colnames(mats$A) == "1")){
 #     IntCol = which(colnames(A) == "1" | colnames(mats$A) != mats$name.factors)
 #     A_minusInt = mats$A[,-IntCol]
  #    A_pen = A_minusInt != 0
  #    pars_pen2 = c(A_minusInt[A_pen],pars_pen2)
  #  }else{
  #    inds2 = mats$A[,colnames(mats$A) != mats$name.factors]
#
   #   pars_pen2 = c(inds2[inds2 != 0],pars_pen2)
  #  }
 # }
}else if(is.null(pars_pen)==FALSE & is.numeric(pars_pen)==FALSE){
  #pars_pen2 <- parse_parameters(pars_pen,model)
  ids = which(mats$pars.align[,2] %in% pars_pen)
  pars_pen2 = as.numeric(mats$pars.align[ids,1])
}else if(is.numeric(pars_pen)){
  pars_pen2 = pars_pen
}#else if(is.null(pars_pen)==TRUE){
#  if(any(colnames(mats$A) == "1")){
#    IntCol = which(colnames(mats$A) == "1")
#    A_minusInt = mats$A[,-IntCol]
#    A_pen = A_minusInt != 0
#    pars_pen2 = A_minusInt[A_pen]
#  }else{
#    A_pen = mats$A != 0
#    pars_pen2 = mats$A[A_pen]
#  }
#}



pars_pen = as.numeric(pars_pen2)



if(is.null(pars_pen) & type!="none"){
  stop("for cv_regsem(), pars_pen needs to be specified")
}

if(fit.ret2 == "test" && is.null(test.n.obs)){
  stop("Please provide a sample size for the test sample")
}


#if(is.null(pars_pen)==FALSE & is.numeric(pars_pen)==FALSE){
#  pars_pen <- parse_parameters(pars_pen,model)
#}
if(quasi==TRUE){
  warnings("The quasi-Newton method is currently not recommended")
}

if(parallel == TRUE){
  stop("parallel is not currently supported")
}

if(type == "scad" | type == "mcp" & jump < 0.1){
  warnings("For both scad and mcp it is recommended to increase jump > 0.1")
}

if(parallel==FALSE){
par.matrix <- matrix(0,n.lambda,length(extractMatrices(model)$parameters))
fits <- matrix(NA,n.lambda,length(fit.ret)+2)
fit.reg <- rep(NA,n.lambda)
fitt.var <- matrix(NA,n.lambda,length(fit.ret))
SHRINK2 = lambda.start
dfs = rep(NA,n.lambda)
count = 0
counts=n.lambda
#res2 <- data.frame(matrix(NA,counts,3))
#coefs = rep(1,14)
if(verbose==TRUE){
  pb <- txtProgressBar(min = 0, max = counts, style = 3)
}


while(count < counts){
  count = count + 1

  # create progress bar
if(verbose==TRUE){
  setTxtProgressBar(pb, count)
}



  SHRINK <- SHRINK2 + jump*(count-1) # 0.01 works well & 0.007 as well with 150 iterations

  if(count > 1 & all(abs(par.matrix[count-1,pars_pen])<.001)){
    break
  }

if(mult.start==FALSE){

  if(warm.start==FALSE | count == 1){
    itt = 0
    Start=Start
  }else if(fits[count-1,2] == 0){
    itt = 0
    Start = par.matrix[count-1,]
    Start[pars_pen] = Start[pars_pen]-step*jump
  }else if(fits[count-1,2] == 99){
    Start=Start
  }else{
    itt = itt + 1
    Start = par.matrix[count-itt-1,]
    Start[pars_pen] = Start[pars_pen]-itt*jump
  }


  if(fit.ret2 == "train" || fit.ret2 == "test"){
    out <- regsem(model=model,lambda=SHRINK,type=type,data=data,
                  random.alpha=random.alpha,
                  optMethod=optMethod,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  gamma=gamma,
                  prerun=prerun,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  round=round,
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
                  random.alpha=random.alpha,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,
                  subOpt=subOpt,
                  alpha=alpha,
                  gamma=gamma,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,prerun=prerun,
                  Start=Start,
                  UB=UB,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  round=round,
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
                  random.alpha=random.alpha,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,
                  subOpt=subOpt,
                  gamma=gamma,
                  alpha=alpha,prerun=prerun,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,
                  Start=Start,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  round=round,
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


      if(mats$mean == TRUE){
        mm = mats$A[,"1"]

        SampMean <- colMeans(test)
        ss = match(names(mm[mm > 0]),model@Data@ov$name)
        SampMean[-c(ss)] = 0

        SampCov=cov(test)
        SampCov2 <- SampCov + SampMean%*%t(SampMean)

        # try changing size of SampCov
        SampCov3 = cbind(SampCov2,SampMean)
        SampCov = rbind(SampCov3,append(SampMean,1))

      }else if(mats$mean == FALSE){
        # SampCov <- model@SampleStats@cov[][[1]]
        SampMean = NULL
        SampCov=cov(test)
      }

      fitt.out = try(fit_indices(out2,CV=TRUE,CovMat=SampCov,n.obs=nrow(test))$fits[fit.ret],silent=T)

        if(inherits(fitt.out, "try-error")) {
          fitt[i,] = NA
        }else{
          fitt[i,] = fitt.out
        }


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
                  random.alpha=random.alpha,
                  gradFun=gradFun,hessFun=hessFun,
                  parallel=parallel,Start=Start,
                  subOpt=subOpt,
                  alpha=alpha,
                  gamma=gamma,
                  pars_pen=pars_pen,
                  diff_par=diff_par,
                  LB=LB,
                  UB=UB,prerun=prerun,
                  par.lim=par.lim,
                  block=block,
                  full=full,
                  calc=calc,
                  tol=tol,
                  round=round,
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
                     random.alpha=random.alpha,
                     gradFun=gradFun,hessFun=hessFun,
                     parallel=parallel,Start=Start,
                     subOpt=subOpt,
                     alpha=alpha,
                     gamma=gamma,
                     pars_pen=pars_pen,
                     diff_par=diff_par,
                     LB=LB,
                     UB=UB,
                     par.lim=par.lim,
                     block=block,prerun=prerun,
                     full=full,
                     calc=calc,
                     tol=tol,
                     round=round,
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


        if(mats$mean == TRUE){
          mm = mats$A[,"1"]

          SampMean <- colMeans(test)
          ss = match(names(mm[mm > 0]),model@Data@ov$name)
          SampMean[-c(ss)] = 0

          SampCov=cov(test)
          SampCov2 <- SampCov + SampMean%*%t(SampMean)

          # try changing size of SampCov
          SampCov3 = cbind(SampCov2,SampMean)
          SampCov = rbind(SampCov3,append(SampMean,1))

        }else if(mats$mean == FALSE){
          # SampCov <- model@SampleStats@cov[][[1]]
          SampMean = NULL
          SampCov=cov(test)
        }

        fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=SampCov,n.obs=nrow(test))$fits[fit.ret]

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
      Start2="lavaan"
    }else if(fits[count-1,2] == 0){
      itt = 0
      Start2 = par.matrix[count-1,]
      Start2[pars_pen] = Start2[pars_pen]-step*jump
    }else if(fits[count-1,2] == 99){
      Start2="lavaan"
    }else{
      itt = itt + 1
      Start2 = par.matrix[count-itt-1,]
      Start2[pars_pen] = Start2[pars_pen]-itt*jump
    }




    if(fit.ret2 == "train" || fit.ret2 == "test"){
      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                      LB=LB,UB=UB,par.lim=par.lim,
                      random.alpha=random.alpha,
                      type=type,optMethod=optMethod,
                      gradFun=gradFun,hessFun=hessFun,
                      tol=tol,
                      round=round,
                      alpha=alpha,
                      gamma=gamma,
                      solver=solver,
                      quasi=quasi,
                      solver.maxit=solver.maxit,
                      max.iter=max.iter,
                      full=full,prerun=prerun,
                      block=block,
                      alpha.inc=alpha.inc,
                      line.search=line.search,
                      step=step,
                      momentum=momentum,
                      Start2=Start2,
                      step.ratio=step.ratio,nlminb.control=nlminb.control,
                      pars_pen=pars_pen,diff_par=diff_par)

    }else if(fit.ret2=="boot"){
      fitt <- matrix(NA,n.boot,length(fit.ret))

      out <- multi_optim(model=model,max.try=multi.iter,lambda=SHRINK,
                         LB=LB,UB=UB,par.lim=par.lim,
                         random.alpha=random.alpha,
                         type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,
                         tol=tol,
                         round=round,
                         alpha=alpha,
                         gamma=gamma,
                         solver=solver,prerun=prerun,
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
                         pars_pen=pars_pen,diff_par=diff_par)


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
                           random.alpha=random.alpha,
                           type=type,optMethod=optMethod,
                           gradFun=gradFun,hessFun=hessFun,
                           tol=tol,
                           round=round,
                           alpha=alpha,prerun=prerun,
                           gamma=gamma,
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
                           pars_pen=pars_pen,diff_par=diff_par)


        if(out$convergence==0){


          if(mats$mean == TRUE){
            mm = mats$A[,"1"]

            SampMean <- colMeans(test)
            ss = match(names(mm[mm > 0]),model@Data@ov$name)
            SampMean[-c(ss)] = 0

            SampCov=cov(test)
            SampCov2 <- SampCov + SampMean%*%t(SampMean)

            # try changing size of SampCov
            SampCov3 = cbind(SampCov2,SampMean)
            SampCov = rbind(SampCov3,append(SampMean,1))

          }else if(mats$mean == FALSE){
            # SampCov <- model@SampleStats@cov[][[1]]
            SampMean = NULL
            SampCov=cov(test)
          }

          fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=SampCov,n.obs=nrow(test))$fits[fit.ret]

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
                         random.alpha=random.alpha,
                         type=type,optMethod=optMethod,
                         gradFun=gradFun,hessFun=hessFun,
                         tol=tol,
                         round=round,
                         alpha=alpha,
                         gamma=gamma,
                         solver=solver,
                         quasi=quasi,prerun=prerun,
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
                            random.alpha=random.alpha,
                            type=type,optMethod=optMethod,
                            gradFun=gradFun,hessFun=hessFun,
                            tol=tol,
                            round=round,
                            alpha=alpha,prerun=prerun,
                            gamma=gamma,
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


          if(mats$mean == TRUE){
            mm = mats$A[,"1"]

            SampMean <- colMeans(test)
            ss = match(names(mm[mm > 0]),model@Data@ov$name)
            SampMean[-c(ss)] = 0

            SampCov=cov(test)
            SampCov2 <- SampCov + SampMean%*%t(SampMean)

            # try changing size of SampCov
            SampCov3 = cbind(SampCov2,SampMean)
            SampCov = rbind(SampCov3,append(SampMean,1))

          }else if(mats$mean == FALSE){
            # SampCov <- model@SampleStats@cov[][[1]]
            SampMean = NULL
            SampCov=cov(test)
          }

          fitt[i,] = fit_indices(out2,CV=TRUE,CovMat=SampCov,n.obs=nrow(test))$fits[fit.ret]

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

    fitt = try(fit_indices(out,CovMat=test.cov,CV=TRUE, n.obs = test.n.obs)$fits[fit.ret],silent=T)
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
  dfs[count] = out$df

  colnames(par.matrix) = names(out$coefficients)
  colnames(fits) <- c("lambda","conv",fit.ret)
  fit.index = fits[,metric]
  conv = fits[,"conv"]
  if(metric=="TLI" | metric=="CFI"){
    loc = which(abs(fit.index)==max(abs(fit.index[conv==0 & is.nan(fit.index) == FALSE & is.na(conv)==FALSE])))[1]
  }else{
    loc = which(abs(fit.index)==min(abs(fit.index[conv==0 & is.nan(fit.index) == FALSE & is.na(conv)==FALSE])))[1]
  }

  final_pars = par.matrix[loc,]

  out2 <- list(par.matrix,fits,final_pars,pars_pen,dfs,metric) #fitt_var
 # ret

}
}else if(parallel==TRUE){

  stop("Parallel is not currently supported")

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
                    random.alpha=random.alpha,
                    gradFun=gradFun,hessFun=hessFun,
                    parallel=parallel,Start=Start,
                    subOpt=subOpt,
                    pars_pen=pars_pen,
                    gamma=gamma,
                    diff_par=diff_par,
                    LB=LB,
                    alpha=alpha,
                    UB=UB,prerun=prerun,
                    calc=calc,
                    nlminb.control=nlminb.control,
                    tol=tol,
                    round=round,
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
                         random.alpha=random.alpha,
                         gradFun=gradFun,hessFun=hessFun,nlminb.control=nlminb.control,
                         tol=tol,
                         round=round,
                         full=full,
                         alpha=alpha,
                         gamma=gamma,
                         block=block,prerun=prerun,
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
      fitt = try(fit_indices(out,CovMat=test.cov,CV=TRUE, n.obs = test.n.obs)$fits[fit.ret],silent=T)
      if(inherits(fitt, "try-error")) {
        fitss = rep(NA,ncol(fits)-2)
      }else{
        fitss = fitt
      }
    }else if(fit.ret2 == "boot"){
      fitt = try(fit_indices(out,CV="boot",n.obs=model@SampleStats@nobs[[1]][1])$fits[fit.ret],silent=T)
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
                     "optMethod","random.alpha",
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
#out2$pars_pen <- pars_pen
out2$call <- match.call()
class(out2) <- "cvregsem"
names(out2) <- c("parameters","fits","final_pars","pars_pen","df","metric","call")#"fit_variance"
out2

#close(pb)

}
