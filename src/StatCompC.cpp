#include <Rcpp.h>
using namespace Rcpp;

//' @title Metropolis samples using Rcpp
//' @description A Metropolis sampler using Rcpp
//' @param sigma the sample variance
//' @param x0 the given initial value
//' @param N the number of the samples
//' @return a random sample of size \code{n}
//' @useDynLib StatComp20039
//' @examples
//' \dontrun{
//' sigma = 2
//' x0 = 25
//' N = 2000
//' rwMetropolis(sigma,x0,N)
//' }
//' @export
//[[Rcpp::export]]
NumericVector rwMetropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-abs(y[0])) / exp(-abs(x[i-1])))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1]; 
    }
  }
  return(x);
} 
