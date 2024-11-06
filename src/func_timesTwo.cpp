// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector& x) {

  NumericVector typeA = Rcpp::RcppArmadillo::sample(x, 1, false,
                                          NumericVector::create(0.9, 0.1));
  std::cout<< typeA<<std::endl;
  return x * 2;

}

