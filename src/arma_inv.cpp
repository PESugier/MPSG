#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
//using namespace arma;

// [[Rcpp::export]]
arma::mat arma_inv (arma::mat M){
  arma::mat  c_inv = inv_sympd(M);
  return c_inv;
}



