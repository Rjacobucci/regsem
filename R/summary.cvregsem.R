#'
#' print information about cvregsem object
#' @param object cv_regsem object
#' @param ... Additional arguments, default fit is fit="BIC"
#' @export
summary.cvregsem <- function(object,...)
{

  if(!exists("fit")){
    fit="BIC"
  }
  fitt = object[[2]][,fit]
  conv = object[[2]][,"conv"]
  lowest.id <- which(fitt==min(fitt[conv!=99 & is.na(conv)==FALSE]))
  lenpar <- length(object$pars_pen)

  sm <- list(lowest.id=lowest.id, lenpar=lenpar, min.lambda=min(object[[2]][,1]),
             max.lambda=max(object[[2]][,1]), lowest.lambda=object[[2]][lowest.id,1])

  class(sm) <- "summary.cvregsem"


  print.summary.cvregsem <- function(x, ...)
  {
    string <- paste0("CV regsem Object\n",
                     " Number of parameters regularized: ",x$lenpar,"\n",
                     " Lambda ranging from ",x$min.lambda," to ",x$max.lambda,"\n",
                     " Lowest Fit Lambda: ", x$lowest.lambda,"\n\n")
    cat(string)
  }



  print(sm)

}




