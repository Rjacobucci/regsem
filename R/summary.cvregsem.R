#'
#' print information about cvregsem object
summary.cvregsem <- function(x, ...) 
{
  lowest.id <- which.min(r1[[2]][,4])
  lenpar <- length(r1$pars_pen)

  sm <- list(lowest.id=lowest.id, lenpar=lenpar, min.lambda=min(r1[[2]][,1]),
             max.lambda=max(r1[[2]][,1]), lowest.lambda=r1[[2]][lowest.id,1])
  
  class(sm) <- "summary.cvregsem"
  return(sm)
}


print.summary.cvregsem <- function(x, ...)
{
  string <- paste0("CV regsem Object\n",
                   " Number of parameters regularized: ",x$lenpar,"\n",
                   " Lambda ranging from ",x$min.lambda," to ",x$max.lambda,"\n",
                   " Lowest BIC Lambda: ", x$lowest.lambda,"\n\n")
  cat(string)
}