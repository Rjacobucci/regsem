
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;



//' Calculates the gradient vector based on Von Oertzen \& Brick, 2014
//'
//' @param par vector with parameters.
//' @param ImpCov expected covariance matrix.
//' @param SampCov Sample covariance matrix.
//' @param Areg A matrix with current parameter estimates.
//' @param Sreg S matrix with current parameter estimates.
//' @param A A matrix with parameter labels.
//' @param S S matrix with parameter labels.
//' @param F F matrix.
//' @param lambda penalty value.
//' @param type2 penalty type.
//' @param pen_vec parameter indicators to be penalized.
//' @param diff_par parameter values to take deviations from.
//'
// [[Rcpp::export]]

arma::vec rcpp_grad_ram(arma::vec par,
                  arma::mat ImpCov,
                  arma::mat SampCov,
                  arma::mat Areg,
                  arma::mat Sreg,
                  arma::mat A,
                  arma::mat S,
                  arma::mat F,
                  double lambda,
                  int type2,
                  arma::vec pen_vec,
                  arma::vec diff_par) {

    double add = 0;
   // double m;
    //m = ImpCov.n_rows;
    arma::vec grad_out; grad_out.zeros(par.n_elem);
    //double deriv15;
    arma::mat deriv15;
    arma::mat A2; A2.zeros(size(A));
    arma::mat S2; S2.zeros(size(S));
    arma::mat I_A = eye(size(A));
    arma::mat B = pinv(I_A - Areg);
    arma::mat I_ImpCov = eye(size(ImpCov));
    arma::mat C = I_ImpCov - pinv(ImpCov) * SampCov;
    arma::mat E = B * Sreg * B.t();

    double Asize = A.n_rows * A.n_cols;
    double Ssize = S.n_rows * S.n_cols;

// ml
    if(type2==0){

      for (double i = 0; i < grad_out.n_elem; i++) {

        arma::mat A2; A2.zeros(size(A));
        arma::mat S2; S2.zeros(size(S));

        for (double j = 0; j < Asize; j++) {
          if (A[j]==i+1) {
            A2[j] = 1;
          }
        }
        for (double j = 0; j < Ssize; j++) {
          if (S[j]==i+1) {
            S2[j] = 1;
          }
        }


        deriv15 = F * B * A2 * E * F.t() + F * B * S2 * B.t() * F.t();
          // left out mean part
        grad_out[i]  = trace(pinv(ImpCov) * deriv15 * C);


      }

    }
// lasso
    else if(type2==1){
      //int add = 0;
      for (double i = 0; i < grad_out.n_elem; i++) {

        arma::mat A2; A2.zeros(size(A));
        arma::mat S2; S2.zeros(size(S));

        for (double j = 0; j < Asize; j++) {
          if (A[j]==i+1) {
            A2[j] = 1;

           // if(any((i+1)==pen_vec)==true) {
          //  add = lambda * ((Areg[j] > 0) - (Areg[j] < 0));
           //add = ((Areg[j]- lambda > 0) - (Areg[j] - lambda < 0));
            }
          //  else{
              add = 0;
           // }


        //  }
        }
        for (double j = 0; j < Ssize; j++) {
          if (S[j]==i+1) {
            S2[j] = 1;
          }
        }


        deriv15 = F * B * A2 * E * F.t() + F * B * S2 * B.t() * F.t();
        // left out mean part


        grad_out[i]  = trace(pinv(ImpCov) * deriv15 * C)  + add;
       // add = 0;

      }

    }
// ridge
    else if(type2==2){
      //int add = 0;
      for (double i = 0; i < grad_out.n_elem; i++) {

        arma::mat A2; A2.zeros(size(A));
        arma::mat S2; S2.zeros(size(S));

        for (double j = 0; j < Asize; j++) {
          if (A[j]==i+1) {
            A2[j] = 1;

            if(any((i+1)==pen_vec)==true) {
              add = 2 * lambda * Areg[j];
            }
            //  else{
            //    add = 0;
            // }


          }
        }
        for (double j = 0; j < Ssize; j++) {
          if (S[j]==i+1) {
            S2[j] = 1;
          }
        }


        deriv15 = F * B * A2 * E * F.t() + F * B * S2 * B.t() * F.t();
        // left out mean part


        grad_out[i]  = trace(pinv(ImpCov) * deriv15 * C)  + add;
        add = 0;

      }

    }
    // diff lasso
    else if(type2==3){
      //int add = 0;
      int count = -1;
      for (double i = 0; i < grad_out.n_elem; i++) {

        arma::mat A2; A2.zeros(size(A));
        arma::mat S2; S2.zeros(size(S));

        for (double j = 0; j < Asize; j++) {
          if (A[j]==i+1) {
            A2[j] = 1;

            if(any((i+1)==pen_vec)==true) {
              count = count + 1;
              double diff = Areg[j] - diff_par[count];
              add = lambda * ((diff > 0) - (diff < 0));
            }
            //  else{
            //    add = 0;
            // }


          }
        }
        for (double j = 0; j < Ssize; j++) {
          if (S[j]==i+1) {
            S2[j] = 1;
          }
        }


        deriv15 = F * B * A2 * E * F.t() + F * B * S2 * B.t() * F.t();
        // left out mean part


        grad_out[i]  = trace(pinv(ImpCov) * deriv15 * C)  + add;
        add = 0;

      }

    }


      int idx1 = max(max(A));
      int idx2 = max(max(S))-1;

       arma::vec grad_out2 = join_cols<mat>(grad_out(span(0,idx1-1)),grad_out(span(idx1,idx2)) * 0.5);
      //grad_out(span(idx1,idx2)) = grad_out(span(idx1,idx2)) * 0.5;

      //return Rcpp::List::create(
       // Rcpp::Named("span") = grad_out);
    //    Rcpp::Named("A")= A);
    //    Rcpp::Named("ImpCov")= ImpCov,
    //    Rcpp::Named("C")= C,
    //    Rcpp::Named("SampCov")= SampCov);
    return grad_out2;

}
