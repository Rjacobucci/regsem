#' Summary results from regsem.
#'
#' @param object An object from regsem.
#' @param ... Other arguments.
#' @return Details regarding convergence and fit
#' @method summary regsem
#' @export


summary.regsem <- function(object,...){
  fits = fit_indices(object,CV=FALSE)

  if(any(is.null(object$mediation))==FALSE){
    TAB <- cbind(convergence = object$convergence,
                 df = object$df,
                 fit=object$fit,
                 rmsea = fits$fit["rmsea"],
                 BIC = fits$fit["BIC"])

    ret <- list(call=object$call,
                estimates = object$coefficients,
                returnVals=TAB,
                mediation=object$mediation)
  }else{
    TAB <- cbind(convergence = object$convergence,
                 df = object$df,
                 fit=object$fit,
                 rmsea = fits$fit["rmsea"],
                 BIC = fits$fit["BIC"])

    ret <- list(call=object$call,
                estimates = object$coefficients,
                returnVals=TAB)
  }



  class(ret) <- "summary.regsem"
  print(ret)
}
