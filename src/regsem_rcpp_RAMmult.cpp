
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


//' Take RAM matrices, multiplies, and returns Implied Covariance matrix.
//'
//' @param par parameter estimates.
//' @param A A matrix with parameter labels.
//' @param S S matrix with parameter labels.
//' @param S_fixed S matrix with fixed indicators.
//' @param A_fixed A matrix with fixed indicators.
//' @param A_est A matrix with parameter estimates.
//' @param S_est S matrix with parameter estimates.
//' @param Fmat Fmat matrix.
//' @param I Diagonal matrix of ones.
//'
//'
// [[Rcpp::export]]
List rcpp_RAMmult(NumericVector par,
          NumericMatrix A,
          NumericMatrix S,
          LogicalMatrix S_fixed,
          LogicalMatrix A_fixed,
          NumericMatrix A_est,
          NumericMatrix S_est,
          IntegerMatrix Fmat,
          IntegerMatrix I) {
//  mat A2;
//  mat S2;
  NumericMatrix A2;
  A2 = clone(A);
  NumericMatrix S2;
  S2 = clone(S);

//  arma::vec i = 1;
double Asize = A2.nrow() * A2.ncol();
for (double i = 0; i < Asize; i++) {
  if (A2[i] > 0) {
    A2[i] = par[A2[i]-1];
  }
}
double Ssize = S2.nrow() * S2.ncol();
for (double i = 0; i < Ssize; i++) {
  if (S2[i] > 0) {
    S2[i] = par[S2[i]-1];
  }
}


//double Asize = A2.nrow() * A2.ncol();
for (double i = 0; i < Asize; i++) {
  if (A_fixed[i] == TRUE) {
    A2[i] = A_est[i];
  }
}
//double Ssize = S2.nrow() * S2.ncol();
for (double i = 0; i < Ssize; i++) {
  if (S_fixed[i] == TRUE) {
    S2[i] = S_est[i];
  }
}



//NumericMatrix ImpCov;
arma::mat A3 = Rcpp::as <arma::mat>(A2);
arma::mat I3 = Rcpp::as<arma::mat>(I);
arma::mat Fmat3 = Rcpp::as <arma::mat>(Fmat);
arma::mat S3 = Rcpp::as<arma::mat>(S2);




arma::mat ImpCov = Fmat3 * pinv(I3-A3) * S3 * (pinv(I3-A3)).t() * Fmat3.t();


return Rcpp::List::create(
  Rcpp::Named("ImpCov") = ImpCov,
  Rcpp::Named("A_est22") = A2,
  Rcpp::Named("S_est22") = S2,
  Rcpp::Named("S") = S,
  Rcpp::Named("A") = A);

}


