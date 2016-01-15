#'
#'
#' uses multiple starts to try and get convergence, running regsem.
#' @param lavObject lavaan output object.
#' @param max.try number of starts to try before convergence.
#' @param SHRINK shrinkage parameter.
#' @param lower lower bound vector.
#' @param upper upper bound vector.
#' @param type type of penalty to use.
#' @param optMethod solver to use.
#' @param gradFun gradient function to use.
#' @param pars_pen parameter indicators to penalize.
#' @param diff_par parameter values to deviate from.
#' @param hessFun hessian function to use.
#' @keywords multiple optim
#' @export
#' @examples
#' \dontrun{
#' multi_optim()
#' }

multi_optim <- function(lavObject,max.try=10,SHRINK,
                         lower=-Inf,upper=Inf,type,optMethod="nlminb",gradFun="timo",
                         pars_pen,diff_par=NULL,hessFun="timo"){

  sss = seq(1,max.try)

    mult_run <- function(lavObject,n.try=5,SHRINK,lower=-Inf,
                         upper=Inf,type,optMethod=optMethod,
                         gradFun=gradFun,n.optim,pars_pen=pars_pen,
                         diff_par=diff_par,hessFun=hessFun){
      mtt = matrix(NA,n.try,3)
      mtt[,3] = round(runif(5,1,99999))
      start.optim=NULL
      n.optim = 0
      while(n.optim < n.try){
        n.optim=n.optim+1
        set.seed(mtt[n.optim,3])
        start.optim = runif(max(parTable(lavObject)$free))
        fit1 <- try(regsem(lavObject,lambda=SHRINK,type=type,optMethod=optMethod,
                           Start=start.optim,gradFun=gradFun,hessFun=hessFun,
                           LB=lower,UB=upper,pars_pen=pars_pen,diff_par=diff_par),silent=T)

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
    iter.optim = iter.optim + 5

    outt = mult_run(lavObject,n.try=10,SHRINK,lower,upper,type,
                    optMethod=optMethod,gradFun,n.optim=iter.optim,
                    pars_pen=pars_pen,diff_par=diff_par,hessFun=hessFun)

      if(all(is.na(outt[,2])==TRUE)){
        return
      }else if(any(na.omit(outt[,2]) == 0)){
        if(any(is.na(outt[,1]) == FALSE & outt[,1] > 0)){
        row = which(outt[,1] == min(outt[which(is.na(outt[,1]) == FALSE & outt[,1] > 0 & outt[,2] == 0)]))[1]
        set.seed(outt[row,3])
        start.optim = runif(max(parTable(lavObject)$free))

        fit1 <- regsem(lavObject,lambda=SHRINK,type=type,optMethod=optMethod,
                       Start=start.optim,gradFun=gradFun,hessFun=hessFun,
                       LB=lower,UB=upper,pars_pen=pars_pen,diff_par=diff_par)
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
      fit1 <- regsem(lavObject,lambda=SHRINK,type=type,optMethod=optMethod,
                     Start="default",gradFun=gradFun,hessFun=hessFun,
                     LB=lower,UB=upper,pars_pen=pars_pen,diff_par=diff_par)
      fit1$convergence <- 99
      fit1
    }

}



