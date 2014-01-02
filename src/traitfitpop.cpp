#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; // done by sourceCpp?


// ** below at least compile when used in 
// cppFunction(code="//code here", depends="RcppArmadillo")

// [[Rcpp::export]]
double Log_W_bar (double zbar, double theta, double Oz2, double gamma){
    return 0.5 * log(gamma * Oz2) - gamma /2 * pow(zbar - theta, 2); //## compare with eqn 2c lande 2009
}

// [[Rcpp::export]]
double W_bar (double zbar, double theta, double Oz2, double gamma, bool LOG){
  if(LOG){
    return Log_W_bar(zbar, theta, Oz2, gamma); 
      }
  else {
    return  sqrt(gamma * Oz2) * exp(- gamma /2 * pow(zbar - theta,2)); 
  }
}

// [[Rcpp::export]]
double  R_bar(double N, double R0, double Wbar, double K, double thetaL){
    return pow(R0, (1 -pow(N/K, thetaL))) * Wbar;
      }

// [[Rcpp::export]]
double  Log_R_bar(double N, double R0, double logWbar, double K, double thetaL){
  return log(R0) *(1 -pow(N/K, thetaL)) + logWbar;
      }

// [[Rcpp::export]]
arma::vec Beta(double gamma, double A, double B, double a, double b, double e_t, double e_plast){
  //## eqn 3b
  double beta1;
  beta1 = -gamma * (a - A + b * e_plast - B * e_t);
  arma::vec beta(2);
  beta(0) =  beta1;
  beta(1) = beta1 * e_plast;
  return beta;
}

// [[Rcpp::export]]
arma::mat Va (arma::vec env, arma::mat GG){
    // env should be column vector
    return trans(env) * GG * env;
}

//- Compute a  white noise environment with a shift
//- @param env previous years environment NOT used in current formulation
//- @param t the time point
//- @param rho the correlation between environment of selection and the cue
//- @param env.args other args:
//- t.jump the location for the jump,
//- delta env change that takes place at t.jump
//- sd noise in env
//- sdc noise in cue
//- @details nada
//- @import MASS
//- @export
// [[Rcpp::export]]
arma::vec Env_shift(int t, List env_args) {
  // returns (selection env, cue env)
  double t_jump = as<double>(env_args["t.jump"]);
  arma::mat sigma = as<arma::mat>(env_args["sigma"]);
  //double rho  = as<double>(env_args["rho"]);
  //double k = as<double>(env_args["k"]);
  arma::vec out(2);

  if( t < t_jump){
    out.zeros();
  } else {
    out.fill(Rcpp::as<double>(env_args["delta"]));
  }
  // Below equivalent to R code:  
  // > out <- mvrnorm(1, mu = c(0,0) , Sigma = matrix(c(sd, sd*sdc*rho, sd*sdc*rho, sdc), nrow=2)
  // generate mvrandom w/ chol,
  // [todo] - break into sep mvrnormarma see http://gallery.rcpp.org/articles/simulate-multivariate-normal
  int ncols = sigma.n_cols;
  arma::mat Z = arma::randn(1, ncols);    
  return out + Z * arma::chol(sigma); // wld use arma::repmat(out, 1, n) for n>1
  // for red noise do somehting like k * env  + sqrt(1 - k^2) * rnorm (n=1, mean=0, sd=sd)) but hard
  // to figure this out with the cue / correlation thing

}


//- Compute trait + demographic change under stabilizing selection as function of environment via Lande Chevin
//- @param t timestep
//- @param X parameters (a, b, env, N)
//- @param p (rho, gamma, A, B, R0, K, thetaL)
//- @param G the (constant) G matrix
//- @param env.fn function to compute environment with signature function(environment, time, ...)
//- @param env.args extra args for env.fn
//- @details function signature for use with _deSolve_ and _iterate_, imposes fitness based on Lande
//- and demography after CL 2010
//- @import deSolve
//- @export

// [[Rcpp::export]]
arma::mat pdLande (int t, arma::rowvec X, List params, 
		   arma::mat G, List env_args) {
  // state vars in X
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  //// env <- X[3]
  // logN = X[4]

  double gamma  = as<double>(params["gamma"]);
  double Oz2  = as<double>(params["Oz2"]);
  arma::rowvec AB  = as<arma::rowvec>(params["AB"]);
  double R0 = as<double>(params["R0"]);
  double K = as<double>(params["K"]);
  double thetaL = as<double>(params["thetaL"]);
  
  arma::vec ee = Env_shift(t, env_args);
  arma::vec env(2);
  env.fill(1.0);
  env(1)= ee(1);
  
  arma::vec z_bar = X.cols(0,1) * env; 
  //double ;
  double log_w_bar = Log_W_bar(z_bar, as<arma::vec>(AB * env), Oz2, gamma); // todo argument ab*env should be double
  
  arma::vec beta = Beta(gamma, AB(0), AB(1), X(0), X(1), ee(0), ee(1));

  // construct the output
  arma::rowvec out(X.size());
  out.cols(0,1) = X.cols(0,1) + G * beta;
  out(2) = env(1);
  out(3) = Log_R_bar(exp(X(4)), R0, log_w_bar, K, thetaL)  + X(4); // demo update based on mean fitness
  return out;
}


