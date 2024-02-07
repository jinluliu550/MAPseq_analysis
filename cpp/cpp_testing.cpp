#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]


double ddirmultinomial(const IntegerVector x,
                       int n,
                       const NumericVector alpha,
                       bool log){
  
  double alpha_sum = Rcpp::sum(alpha);
  double log_prob = lgamma(alpha_sum) + lgamma(n+1) - lgamma(n+alpha_sum) + Rcpp::sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x+1));
  
  if(log == TRUE){
    
    return log_prob;
    
  }else{
    
    return exp(log_prob);
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
