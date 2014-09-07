
## > sessionInfo()
## R version 3.0.2 (2013-09-25)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)

## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     

## other attached packages:
## [1] RcppArmadillo_0.4.320.0 Rcpp_0.11.2             devtools_1.5           

## loaded via a namespace (and not attached):
##  [1] compiler_3.0.2 digest_0.6.3   evaluate_0.4.3 httr_0.2       memoise_0.1    parallel_3.0.2 RCurl_1.95-4.1 stringr_0.6.2  tools_3.0.2    whisker_0.3-2 




### working examples of how to use attributes with rcpparmadillo
library(Rcpp)
library(RcppArmadillo)
cppFunction('double testarma(arma::vec x, arma::rowvec y) {
      return as_scalar(y * x);
  }',depends='RcppArmadillo')

## I get one warning:
## ld: warning: directory '/usr/local/lib/x86_64' following -L not found
testarma(1:2, 2:3)


## below also works with  sessionInfo() as above (this code is example from ?Rcpp::cppFunction)
cppFunction(depends = "RcppArmadillo",
            'List fastLm(NumericVector yr, NumericMatrix Xr) {
             
             int n = Xr.nrow(), k = Xr.ncol();
             
             arma::mat X(Xr.begin(), n, k, false); 
             arma::colvec y(yr.begin(), yr.size(), false);
             
             arma::colvec coef = arma::solve(X, y);
             arma::colvec resid = y - X*coef;
             
             double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k) );
             arma::colvec stderrest = arma::sqrt(
                 sig2 * arma::diagvec(arma::inv(arma::trans(X)*X)));
             
             return List::create(Named("coefficients") = coef,
                 Named("stderr")       = stderrest
             );
         }')


## alternatively, you can use sourceCpp

sourceCpp(code='
#include "RcppArmadillo.h" // needed after version 4.3
// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp; // not necessary with sourceCpp but include this if using in a package

// [[Rcpp::export]]
double testarma(arma::vec x, arma::rowvec y) {
      return as_scalar(y * x);
  }')

## package: 
## to get started see
## Rcpp::compileAttributes
## and
## Rcpp::Rcpp.package.skeleton (or  Rcpp::RcppArmadillo.package.skeleton if using arma)


## more examples 
library(Rcpp)
sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector mypnorm(Rcpp::NumericVector x) {
    
   int n = x.size();
   Rcpp::NumericVector y(n);

   for (int i=0; i<n; i++) 
      y[i] = R::pnorm(x[i], 0.0, 1.0, 1, 0);

   return y;
}')


sourceCpp(code='
// [[Rcpp::depends(RcppGSL)]]

#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// [[Rcpp::export]]
Rcpp::NumericVector colNorm(Rcpp::NumericMatrix sM) {

    RcppGSL::matrix<double> M(sM); 	// create gsl data structures from SEXP
    int k = M.ncol();
    Rcpp::NumericVector n(k); 		// to store results 

    for (int j = 0; j < k; j++) {
        RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
	n[j] = gsl_blas_dnrm2(colview);
    }
    M.free() ;                          // important: GSL wrappers use C structure
    return n;				// return vector  
}')




 code <- '
    arma::mat GG = Rcpp::as<arma::mat>(G);
    arma::colvec env = Rcpp::as<arma::colvec>(e);
     // env should be column vector
    return Rcpp::wrap(trans(env) * GG * env);
  '

Va <- cppFunction(signature(G="numeric",e="numeric"),
                         code,plugin="RcppArmadillo")
