#'
#' print information about cvregsem object
summary.cvregsem <- function(x, ...) 
{
  lowest.id <- which.min(r1[[2]][,4])
  lenpar <- length(r1$pars_pen)
  string <- paste0("CV regsem Object\n",
    " Number of parameters regularized: ",lenpar,"\n",
    " Lambda ranging from ",min(r1[[2]][,1])," to ",max(r1[[2]][,1]),"\n",
    " Lowest BIC Lambda: ", r1[[2]][lowest.id,1],"\n\n")
  cat(string)
}