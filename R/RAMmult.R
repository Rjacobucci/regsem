#'
#'
#' Take RAM matrices, multiplies, and returns Implied Covariance matrix.
#' @param par parameter estimates.
#' @param A A matrix with parameter labels.
#' @param S S matrix with parameter labels.
#' @param F F matrix.
#' @param A_fixed A matrix with fixed indicators.
#' @param A_est A matrix with parameter estimates.
#' @param S_fixed S matrix with fixed indicators.
#' @param S_est S matrix with parameter estimates.


RAMmult <- function(par,A,S,F,A_fixed,A_est,S_fixed,S_est){

  A2 <- A
  S2 <- S
  # doesn't work for
  for(i in 1:length(par)){
    A2[A2== i] <- par[i]
    S2[S2== i] <- par[i]
  }

  A2[A_fixed] <- A_est[A_fixed]
  S2[S_fixed] <- S_est[S_fixed]


  I = diag(nrow(A))

  ImpCov = F %*% solve(I-A2) %*% S2 %*% t(solve(I-A2)) %*% t(F)

  res <- list()
  res$ImpCov <- ImpCov;
  res$A_est22 <- A2;
  res$S_est22 <- S2
 # res$S <- S; res$S_fixed <- S_fixed;
 # res$A_fixed <- A_fixed; res$F <- F
  res
}

