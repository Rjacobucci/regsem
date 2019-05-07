#'
#' print information about cvregsem object
#' @param object cv_regsem object
#' @param ... Additional arguments
#' @method summary cvregsem
#' @export
summary.cvregsem <- function(object,...)
{


    fit=object$metric

  fitt = object$fits[,fit]
  conv = object$fits[,"conv"]
  lowest.id <- which(fitt==min(fitt[conv==0 & is.nan(fitt)==FALSE & is.na(conv)==FALSE]))
  lenpar <- length(object$pars_pen)
  final_pars <- object$final_pars

  sm <- list(lowest.id=lowest.id, lenpar=lenpar, min.lambda=min(object$fits[,1],na.rm=TRUE),
             max.lambda=max(object$fits[,1],na.rm=TRUE), lowest.lambda=object$fits[lowest.id,1],
             metric=fit,
             final_pars=final_pars,
             num.conv = sum(conv==0,na.rm=TRUE))

  class(sm) <- "summary.cvregsem"


  print.summary.cvregsem <- function(x, ...)
  {
    string <- paste0("CV regsem Object\n",
                     " Number of parameters regularized: ",x$lenpar,"\n",
                     " Lambda ranging from ",x$min.lambda," to ",x$max.lambda,"\n",
                     " Lowest Fit Lambda: ", x$lowest.lambda,"\n",
                     " Metric: ", x$metric,"\n",
                     " Number Converged: ", x$num.conv,"\n\n",collapse=", ")
    cat(string)
   # string2 <- paste0(c("Final Parameters: ", names(x$final_pars)))
   # cat(string2)
   # string3 <- paste0(c("Final Parameters: ", x$final_pars))
   # cat(string3)
  }



  print(sm)
  #return(sm)

}




