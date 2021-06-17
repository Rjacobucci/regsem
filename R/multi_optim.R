#'
#'
#' Multiple starts for Regularized Structural Equation Modeling
#'
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
#' @param max.try number of starts to try before convergence.
#' @param lambda Penalty value. Note: higher values will result in additional
#'        convergence issues.
#' @param alpha Mixture for elastic net.
#' @param gamma Additional penalty for MCP and SCAD
#' @param random.alpha Alpha parameter for randomised lasso. Has to be between
#'        0 and 1, with a default of 0.5. Note this is only used for
#'        "rlasso", which pairs with stability selection.
#' @param LB lower bound vector. Note: This is very important to specify
#'        when using regularization. It greatly increases the chances of
#'        converging.
#' @param UB Upper bound vector
#' @param par.lim Vector of minimum and maximum parameter estimates. Used to
#'        stop optimization and move to new starting values if violated.
#' @param block Whether to use block coordinate descent
#' @param full Whether to do full gradient descent or block
#' @param type Penalty type. Options include "none", "lasso",
#'        "enet" for the elastic net,
#'        "alasso" for the adaptive lasso
#'        and "diff_lasso". If ridge penalties are desired, use type="enet" and
#'        alpha=1. diff_lasso penalizes the discrepency between
#'        parameter estimates and some pre-specified values. The values
#'        to take the deviation from are specified in diff_par. Two methods for
#'        sparser results than lasso are the smooth clipped absolute deviation,
#'        "scad", and the minimum concave penalty, "mcp". Last option is "rlasso"
#'        which is the randomised lasso to be used for stability selection.
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
#' @param pars_pen Parameter indicators to penalize. There are multiple ways to specify.
#'        The default is to penalize all regression parameters ("regressions"). Additionally,
#'        one can specify all loadings ("loadings"), or both c("regressions","loadings").
#'        Next, parameter labels can be assigned in the lavaan syntax and passed to pars_pen.
#'        See the example.Finally, one can take the parameter numbers from the A or S matrices and pass these
#'        directly. See extractMatrices(lav.object)$A.
#' @param diff_par Parameter values to deviate from. Only used when
#'        type="diff_lasso".
#' @param hessFun Hessian function to use. Currently not recommended.
#' @param prerun Logical. Use rsolnp to first optimize before passing to
#'        gradient descent? Only for use with coord_desc.
#' @param tol Tolerance for coordinate descent
#' @param round Number of digits to round results to
#' @param solver Whether to use solver for coord_desc
#' @param quasi Whether to use quasi-Newton. Currently not recommended.
#' @param solver.maxit Max iterations for solver in coord_desc
#' @param alpha.inc Whether alpha should increase for coord_desc
#' @param line.search Use line search for optimization. Default is no, use fixed step size
#' @param step Step size
#' @param momentum Momentum for step sizes
#' @param step.ratio Ratio of step size between A and S. Logical
#' @param verbose Whether to print iteration number.
#' @param warm.start Whether start values are based on previous iteration.
#'        This is not recommended.
#' @param Start2 Provided starting values. Not required
#' @param nlminb.control list of control values to pass to nlminb
#' @param max.iter Number of iterations for coordinate descent
#' @keywords multiple optim
#' @export
#' @examples
#' \dontrun{
#' # Note that this is not currently recommended. Use cv_regsem() instead
#' library(regsem)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[ ,7:15]))
#' mod <- '
#' f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#' outt = cfa(mod, HS, meanstructure=TRUE)
#'
#' fit1 <- multi_optim(outt, max.try=40,
#'                    lambda=0.1, type="lasso")
#'
#'
#'# growth model
#'model <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
#'           s =~ 0*t1 + s1*t2 + s2*t3 + 3*t4 '
#'fit <- growth(model, data=Demo.growth)
#'summary(fit)
#'fitmeasures(fit)

#'fit3 <- multi_optim(fit, lambda=0.2, type="lasso")
#'summary(fit3)
#'}


multi_optim <- function(model,max.try=10,lambda=0,
                        alpha=.5,
                        gamma=3.7,
                        random.alpha=0.5,
                         LB=-Inf,
                        UB=Inf,
                        par.lim=c(-Inf,Inf),
                        block=TRUE,
                        full=TRUE,
                        type="lasso",
                        optMethod="rsolnp",
                        gradFun="ram",
                        pars_pen="regressions",
                        diff_par=NULL,
                        hessFun="none",
                        tol=1e-5,
                        round=3,
                        solver=FALSE,
                        quasi=FALSE,
                        solver.maxit=50000,
                        alpha.inc=FALSE,
                        line.search=FALSE,
                        prerun=FALSE,
                        step=.1,
                        momentum=FALSE,
                        step.ratio=FALSE,
                        verbose=FALSE,warm.start=FALSE,Start2=NULL,
                        nlminb.control=NULL,
                        max.iter=500){


#  if(optMethod=="default" & type=="lasso"){
#    optMethod<-"coord_desc"
#  }

#  if(optMethod=="default" & type!="lasso"){
#    optMethod <- "nlminb"
#  }

 # warning("Note it is not currently recommended to use multi_optim")

#  if(gradFun=="norm"){
#    stop("Only recommended grad function is ram or none at this time")
#  }

#  if(type=="ridge" & gradFun != "none"){
#    warning("At this time, only gradFun=none recommended with ridge penalties")
#  }

#  if(type=="lasso" & gradFun != "ram"){
#    warning("At this time, only gradFun=ram recommended with lasso penalties")
#  }

  if(type=="ridge"){
    type="enet";alpha=1
  }

  if(quasi==TRUE){
    warnings("The quasi-Newton method is currently not recommended")
  }
#if(warm.start==TRUE){
#  stop("warm start not currently functioning")
#}
    if(warm.start==TRUE & is.null(Start2)){
     # fit99 <- suppressWarnings(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
      #                                 Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
      #                                 LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol))
      Start2 = rep(0.5,length(extractMatrices(model)$parameters)) +
        rnorm(length(extractMatrices(model)$parameters),0,0.05)
#
 #     if(length(start.vals) != length(extractMatrices(model)$parameters)){
 #       start.vals = rep(0.5,length(extractMatrices(model)$parameters)) +
 #         rnorm(length(extractMatrices(model)$parameters),0,0.05)
  #    }

    }

  mats <- extractMatrices(model)

  val1 = max(mats$A)
  val2 = max(mats$S) - max(mats$A)

  sss = seq(1,max.try)

    mult_run <- function(model,n.try=1,lambda,LB=-Inf,tol,alpha,gamma,
                         UB=Inf,
                         par.lim,
                         random.alpha,
                         block,
                         full,
                         type,optMethod,warm.start,
                         solver,
                         quasi,
                         round,
                         solver.maxit,
                         alpha.inc,
                         step,
                         momentum,prerun,
                         max.iter,
                         line.search,
                         step.ratio,
                         gradFun,n.optim,pars_pen,nlminb.control,
                         diff_par,hessFun,Start2){
      mtt = matrix(NA,n.try,3)
      mtt[,3] = round(runif(n.try,1,99999))
      start.optim=NULL
      n.optim = 0
      while(n.optim < n.try){
        n.optim=n.optim+1
        set.seed(mtt[n.optim,3])

        if(warm.start==FALSE){
          if(is.null(Start2)){



            #start.optim = c(rep(0,val1) + rnorm(val1,0,0.1),abs(rep(0.5,val2) + rnorm(val2,0,0.1)))
            start.optim=mats$parameters

          }else{
            start.optim = Start2
          }
        }else if(warm.start==TRUE){
          start.optim= mats$parameters
        }

        #else if(warm.start==TRUE){
         # start.optim <- rep(0,length(start.vals))
         # for(i in 1:length(start.vals)){
         #   start.optim[i] <- as.numeric(start.vals[i]) + rnorm(1,0,0.02)
         # }
       # }


        fit1 <- suppressWarnings(try(regsem(model,lambda=lambda,alpha=alpha,gamma=gamma,
                                            type=type,optMethod=optMethod,
                                            Start=start.optim,gradFun=gradFun,hessFun=hessFun,
                                            nlminb.control=nlminb.control,tol=tol,
                                            solver=solver,
                                            random.alpha=random.alpha,
                                            quasi=quasi,
                                            block=block,
                                            round=round,
                                            par.lim=par.lim,
                                            full=full,prerun=prerun,
                                            solver.maxit=solver.maxit,
                                            alpha.inc=alpha.inc,
                                            max.iter=max.iter,
                                            line.search=line.search,
                                            step=step,
                                            momentum=momentum,
                                            step.ratio=step.ratio,
                                            LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par),silent=T))



        if(inherits(fit1, "try-error")) {
          mtt[n.optim,1] = NA
          mtt[n.optim,2] = NA
          start.vals = NULL

          if(warm.start == TRUE){
            start.vals <- mats$parameters +
              rnorm(length(mats$parameters),0,0.05)
          }
        }else{
          start.vals = NULL

          if(warm.start==TRUE){
            start.vals <- as.numeric(fit1$coefficients) + rnorm(length(mats$parameters),0,0.0001)
          }

          if(optMethod=="nlminb"){
            mtt[n.optim,1] = fit1$out$objective
            mtt[n.optim,2] = fit1$out$convergence
          }else{
            #print(fit1$out$value)
            mtt[n.optim,1] = fit1$optim_fit
            mtt[n.optim,2] = fit1$convergence
          }
        }
    }
    ret =  list(mtt=mtt,start.vals=start.vals,fit1=fit1)
    ret
  }


iter.optim = 0
while(iter.optim < max.try){
iter.optim = iter.optim + 1


    ret.mult = mult_run(model,n.try=1,lambda=lambda,alpha=alpha,gamma=gamma,
                        LB,UB,type,warm.start=warm.start,
                        nlminb.control=nlminb.control,tol=tol,
                        solver=solver,
                        random.alpha=random.alpha,
                        quasi=quasi,
                        block=block,
                        round=round,
                        par.lim=par.lim,
                        full=full,prerun=prerun,
                        solver.maxit=solver.maxit,
                        max.iter=max.iter,
                        alpha.inc=alpha.inc,
                        line.search=line.search,
                        step=step,
                        momentum=momentum,
                        step.ratio=step.ratio,
                    optMethod=optMethod,gradFun=gradFun,n.optim=iter.optim,Start2=Start2,
                    pars_pen=pars_pen,diff_par=diff_par,hessFun=hessFun)




    if(warm.start == TRUE & ret.mult$mtt[1,1] > 0){

      Start2 = ret.mult$start.vals +
        c(rep(0,val1) + rnorm(val1,0,0.01),abs(rep(0,val2) + rnorm(val2,0,0.001)))
    }else{

      Start2 = c(rep(0,val1) + rnorm(val1,0,0.1),abs(rep(0.5,val2) + rnorm(val2,0,0.1)))
    }

    outt = ret.mult$mtt


    if(verbose==TRUE) print(c(iter.optim,outt[,1],outt[,2]))

 # if(all(is.na(outt[,2])==TRUE)){
   # print(19)
   #     return(NA)
#
   #}else
     if(any(na.omit(outt[,2]) == 0)){

    #if(any(is.na(outt[,1]) == FALSE & outt[,1] < 999999 & outt[,1] > 0)){
       # row = which(outt[,1] == min(outt[which(is.na(outt[,1]) == FALSE & outt[,1] > 0 & outt[,2] == 0)]))[1]
       # set.seed(outt[row,3])
      #  if(warm.start==FALSE){
       #   start.optim = rep(0.5,length(extractMatrices(model)$parameters)) +
       #     rnorm(length(extractMatrices(model)$parameters),0,0.05)
       # }else if(warm.start==TRUE){
       #   start.optim = Start2
       # }

        #fit1 <- suppressWarnings(regsem(model,lambda=lambda,type=type,optMethod=optMethod,
        #               Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
         #              LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol))

        return(ret.mult$fit1)

        break
      # }else{
      #    return(NA)
      # }
    }

}



    if(exists("fit1")==FALSE){

      if(warm.start==TRUE){
        Start=Start2
      }else{
        Start="default"
      }

      fit1 <- suppressWarnings(regsem(model,lambda=lambda,
                     alpha=alpha,gamma=gamma,
                     random.alpha=random.alpha,
                     Start=Start,gradFun=gradFun,hessFun=hessFun,
                     nlminb.control=nlminb.control,tol=tol,
                     solver=solver,
                     quasi=quasi,
                     block=block,
                     full=full,
                     round=round,
                     par.lim=par.lim,prerun=prerun,
                     max.iter=max.iter,
                     solver.maxit=solver.maxit,
                     alpha.inc=alpha.inc,
                     step=step,
                     line.search=line.search,
                     momentum=momentum,
                     step.ratio=step.ratio,
                     LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par))

        fit1$convergence <- 99
        return(fit1)
    }



    if(fit1$convergence != 0){
      warning("WARNING: Model did not converge! It is recommended to increase max.try")
    }

}



