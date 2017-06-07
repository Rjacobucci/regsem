#'
#' print information about cvregsem object
#' @param x cv_regsem
#' @param fit Which fit index to use
#' @export
summary.cvregsem <- function(x,fit="BIC",...)
{
  lowest.id <- which.min(x[[2]][,fit])
  lenpar <- length(x$pars_pen)

  sm <- list(lowest.id=lowest.id, lenpar=lenpar, min.lambda=min(x[[2]][,1]),
             max.lambda=max(x[[2]][,1]), lowest.lambda=x[[2]][lowest.id,1])

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




