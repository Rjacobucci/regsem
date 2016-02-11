#'
#'
#' returns final hessian matrix
#' @param final_par parameter estimates.
#' @param A A matrix with parameter labels.
#' @param S S matrix with parameter labels.
#' @param F F matrix.
#' @param A_fixed A matrix with fixed indicators.
#' @param A_est A matrix with parameter estimates.
#' @param S_fixed S matrix with fixed indicators.
#' @param S_est S matrix with parameter estimates.



ret_hess <- function(final_par,A,S,F,A_fixed,A_est,S_fixed,S_est){
  mult = RAMmult(par=final_par,A,S,F,A_fixed,A_est,S_fixed,S_est)
  retH = hessian(par=final_par,ImpCov=mult$ImpCov,A,A_fixed,A_est,
                          S,S_fixed,S_est,F)
  retH
}

