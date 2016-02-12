#' Calculates the trace of a matrix.
#'
#' @param mat A matrix.
#' @return The sum of the diagonal entries in \code{mat}.
#' @keywords trace
#' @examples
#' \dontrun{
#' mm = matrix(c(1,2,3,4),2,2)
#' trace(mm)
#' }
#'

trace <- function(mat){
  return(sum(diag(mat)))
}
