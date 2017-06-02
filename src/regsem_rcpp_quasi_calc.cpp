
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Compute quasi Hessian
//'
//' @param I identity matrix.
//' @param s s vector.
//' @param p p vector.
//' @param y y vector.
//' @param H previous Hessian.
// [[Rcpp::export]]
List rcpp_quasi_calc(NumericMatrix I,
                NumericVector s,
                NumericVector p,
                NumericVector y,
                NumericMatrix H) {

  arma::vec p2 = Rcpp::as <arma::vec>(p);
  arma::vec s2 = Rcpp::as<arma::vec>(s);
  arma::vec y2 = Rcpp::as <arma::vec>(y);
  arma::mat HH = Rcpp::as <arma::mat>(H);
  arma::mat I2 = Rcpp::as<arma::mat>(I);


  arma::mat H2 = I2 - (p2 = s2.each_row() *y2.t())*HH*I2 - p2 % y2*s2.t() + p2 % s2*s2.t();

    return Rcpp::List::create(
      Rcpp::Named("H2") = H2);

}
