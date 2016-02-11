#'
#'
#' uses multiple starts to try and get convergence, running regsem.
#' @param lavObject lavaan output object.
#' @param max.try number of starts to try before convergence.
#' @param lambda shrinkage parameter.
#' @param LB LB bound vector.
#' @param UB UB bound vector.
#' @param type type of penalty to use.
#' @param optMethod solver to use.
#' @param gradFun gradient function to use.
#' @param pars_pen parameter indicators to penalize.
#' @param diff_par parameter values to deviate from.
#' @param hessFun hessian function to use.
#' @param verbose Whether to print iteration number.
#' @param warm.start Whether start values are based on previous iteration.
#' @param tol absolute tolerance for convergence.
#' @param max.iter max iterations for optimization.
#' @keywords multiple optim
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
#' fit1 <- multi_optim(outt,max.try=40,
#'                    lambda=0.1,type="lasso",
#'                    optMethod="nlminb", gradFun="timo")
#'}


multi_optim <- function(lavObject,max.try=10,lambda,
                         LB=-Inf,UB=Inf,type,optMethod="nlminb",gradFun="timo",
                         pars_pen=NULL,diff_par=NULL,hessFun="none",
                        verbose=TRUE,warm.start=TRUE,
                        tol=1e-6,max.iter=150){



    if(warm.start==TRUE){
      fit99 <- suppressWarnings(regsem(lavObject,lambda=lambda,type=type,optMethod=optMethod,
                                      Start="default",gradFun=gradFun,hessFun=hessFun,tol=tol,
                                      LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par))
      start.vals <- as.vector(fit99$par.ests)
    }

  sss = seq(1,max.try)

    mult_run <- function(lavObject,n.try=1,lambda,LB=-Inf,
                         UB=Inf,type,optMethod,
                         gradFun,n.optim,pars_pen,
                         diff_par,hessFun){
      mtt = matrix(NA,n.try,3)
      mtt[,3] = round(runif(n.try,1,99999))
      start.optim=NULL
      n.optim = 0
      while(n.optim < n.try){
        n.optim=n.optim+1
        set.seed(mtt[n.optim,3])

        if(warm.start==FALSE){
          start.optim = runif(max(parTable(lavObject)$free))
        }else if(warm.start==TRUE){
          start.optim <- rep(0,length(start.vals))
          for(i in 1:length(start.vals)){
            start.optim[i] <- as.numeric(start.vals[i]) + rnorm(1,0,0.02)
          }
        }


        fit1 <- suppressWarnings(try(regsem(lavObject,lambda=lambda,type=type,optMethod=optMethod,
                           Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
                           LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol),silent=T))

        if(inherits(fit1, "try-error")) {
          mtt[n.optim,1] = NA
          mtt[n.optim,2] = NA
        }else{
          mtt[n.optim,1] = fit1$optim_fit
          mtt[n.optim,2] = fit1$convergence
        }
      }
      mtt
    }

    iter.optim = 0
    while(iter.optim < max.try){
    iter.optim = iter.optim + 1

    outt = mult_run(lavObject,n.try=1,lambda=lambda,LB,UB,type,
                    optMethod=optMethod,gradFun=gradFun,n.optim=iter.optim,
                    pars_pen=pars_pen,diff_par=diff_par,hessFun=hessFun)

    if(verbose==TRUE) print(c(iter.optim,outt[,1],outt[,2]))

      if(all(is.na(outt[,2])==TRUE)){
        return
      }else if(any(na.omit(outt[,2]) == 0)){
        if(any(is.na(outt[,1]) == FALSE & outt[,1] > 0)){
        row = which(outt[,1] == min(outt[which(is.na(outt[,1]) == FALSE & outt[,1] > 0 & outt[,2] == 0)]))[1]
        set.seed(outt[row,3])
        start.optim = runif(max(parTable(lavObject)$free))

        fit1 <- suppressWarnings(regsem(lavObject,lambda=lambda,type=type,optMethod=optMethod,
                       Start=start.optim,gradFun=gradFun,hessFun=hessFun,max.iter=max.iter,
                       LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par,tol=tol))
        fit1
        break
        }else{
          return
        }
      }else{
        return
      }
    }
    if(exists("fit1")==TRUE){
      fit1
    }else{
      fit1 <- suppressWarnings(regsem(lavObject,lambda=lambda,type=type,optMethod=optMethod,
                     Start="default",gradFun=gradFun,hessFun=hessFun,tol=tol,
                     LB=LB,UB=UB,pars_pen=pars_pen,diff_par=diff_par))
      fit1$convergence <- 99
      fit1
    }



    if(fit1$convergence != 0){
      warning("WARNING: Model did not converge! It is recommended to increase max.try")
    }

}



