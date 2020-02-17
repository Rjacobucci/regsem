
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


//' Calculates the objective function values.
//'
//' @param ImpCov expected covariance matrix.
//' @param SampCov Sample covariance matrix.
//' @param type2 penalty type.
//' @param lambda penalty value.
//' @param gamma additional penalty for mcp and scad
//' @param pen_vec vector of penalized parameters.
//' @param pen_diff Vector of values to take deviation from.
//' @param e_alpha Alpha for elastic net
//' @param rlasso_pen Alpha for rlasso2
//'
// [[Rcpp::export]]
double rcpp_fit_fun(Rcpp::NumericMatrix ImpCov,
                 Rcpp::NumericMatrix SampCov,
                 int type2,
                 double lambda,
                 double gamma,
                 arma::vec pen_vec,
                 arma::vec pen_diff,
                 double e_alpha,
                 double rlasso_pen){//,
                 //int estimator2,
                // arma::vec poly_vec,
                // arma::vec imp_vec){
           //      int alpha,


double m;
m = ImpCov.nrow();
double fit;
double add;
double fit_base;
arma::mat ImpCov2 = Rcpp::as <arma::mat>(ImpCov);
arma::mat SampCov2 = Rcpp::as<arma::mat>(SampCov);
//arma::mat Areg2 = Rcpp::as<arma::mat>(Areg);
//arma::double lambda2 = Rcpp::as<arma::int>(lambda);
//arma::vec pen_vec2 = Rcpp::as<arma::vec>(pen_vec);


//NumericMatrix A2;
//A2 = clone(A);
//NumericMatrix S2;
//S2 = clone(S);

//if (estimator2 == 1) { // ML
  fit_base = log(det(ImpCov2)) + trace(SampCov2 * (pinv(ImpCov2))) - log(det(SampCov2))  - m;
//}
//else if (estimator2 == 2) {
//   fit_base = as_scalar((poly_vec - imp_vec).t() * (poly_vec - imp_vec));
//}



  if (type2 == 0) {
    fit = fit_base;
  }
  else if (type2 == 1) {
    fit = fit_base  + lambda * norm(pen_vec,1);
  }
  else if (type2 == 2) {
    fit = fit_base  + lambda * norm(pen_vec,2);
  }
  else if (type2 == 3) {
    fit = fit_base  + lambda * norm(pen_diff,1);
  }
  else if (type2 == 4) {
    fit = fit_base  + lambda * ((1-e_alpha)*norm(pen_vec,1) + e_alpha*norm(pen_vec,2));
  }
  else if (type2 == 5) {
    fit = fit_base  + rlasso_pen;
  }
  else if (type2 == 6) {
    add = 0;
    for (double i = 0; i < pen_vec.n_elem; i++) {
      if (std::abs(pen_vec[i]) <= lambda){
        add = add + lambda*(std::abs(pen_vec[i]));
      }
      else if(std::abs(pen_vec[i]) > lambda && std::abs(pen_vec[i]) <= lambda*gamma){
        add = add - (lambda*lambda + std::abs(pen_vec[i])*std::abs(pen_vec[i]) + 2*gamma*lambda*(std::abs(pen_vec[i])))/(2*(gamma-1));
      }
      else if(std::abs(pen_vec[i]) > lambda*gamma){
// add = add + ((lambda*lambda)*(gamma*gamma - 1))/(2*(gamma-1));
        add = add + ((gamma+1)*(lambda*lambda))/2;
        }
        else {
          add = add + 0;
        }
      }
    fit = fit_base + add;
  }
  else if (type2 == 7) {
    add = 0;
    for (double i = 0; i < pen_vec.n_elem; i++) {
      if (std::abs(pen_vec[i]) <= lambda*gamma){
        add = add + (lambda * (std::abs(pen_vec[i]) - ((pen_vec[i]*pen_vec[i])/(2*lambda*gamma))));
      }
      else if(std::abs(pen_vec[i]) > lambda * gamma){
        add = add + ((lambda*lambda)*gamma)/2;
      }
    }
    fit = fit_base + add;
  }
  else{
    fit = -9999999;
  }

//  List ret;
//  ret["fit"] = fit;
 // double fitt = fit;
 // return fit;
  //return fit;

  //return(std::cout << std::setprecision(5) << fit);
  return 0.5*fit;
}

