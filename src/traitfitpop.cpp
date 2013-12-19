// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp; // done by sourceCpp

//' Compute trait + demographic change under stabilizing selection as function of environment via Lande Chevin
//' @param t timestep
//' @param X parameters (a, b, env, N)
//' @param p (rho, gamma, A, B, R0, K, thetaL)
//' @param G the (constant) G matrix
//' @param env.fn function to compute environment with signature function(environment, time, ...)
//' @param env.args extra args for env.fn
//' @details function signature for use with deSolve and 'iterate', imposes fitness based on Lande
//' and demography after CL 2010
//' @import deSolve
//' @export

// [[Rcpp::export]]

arma::mat pdLande (int t, Arma::Vec X, List params, 
		       NumericMatrix G, List env_args) {
  // state vars in X
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  // env <- X[3]
  // N = X[4]

  double t_jump = as<double>(env_args["t.jump"]);
  double k = as<double>(env_args["k"]);
  arma::mat sigma = as<arma::mat>(env_args["sigma"])

  double rho  = as<double>(params['rho']);
  double gamma  = as<double>(params['gamma']);
  double Oz2  = as<double>(params['Oz2']);
  double A  = as<double>(params['A']);
  double B  = as<double>(params['B']);
  double R0 = as<double>(params['R0']);
  double K = as<double>(params['K']);
  double thetaL = as<double>(params['thetaL']);
  
  arma::colvec EE = Env_shift(t, rho, sigma, t_jump, k);

  env <- EE[1]
  e_plast <- EE[2]
  
  z_bar <- a + b * e_plast
  log_W_bar <- W_bar(z_bar, A + B*env, Oz2, gamma)
  N <- R_bar(N, R0, exp(log_W_bar), K, thetaL) * N // demo update based on mean fitness
  
  beta <- Beta(gamma, A, B, a, b, env, e_plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e_plast, G) 
  list(c(X, env, N, va))
}



// ** below at least compile when used in 
// cppFunction(code='//code here ', depends='RcppArmadillo')

// [[Rcpp::export]]
double W_bar (double zbar, double theta, double Oz2, double gamma, bool LOG){
  if(LOG){
    return 0.5 * log(gamma * Oz2) - gamma /2 * pow(zbar - theta, 2); //## compare with eqn 2c lande 2009
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
arma::colvec Beta(double gamma, double A, double B, double a, double b, double e_t, double e_plast){
  //## eqn 3b
  double beta1;
  beta1 = -gamma * (a - A + b * e_plast - B * e_t);
  arma::colvec beta(2);
  beta(1) =  beta1;
  beta(2) = beta1 * e_plast;
  return beta;
}

// [[Rcpp::export]]
arma::mat Va (arma::colvec env, arma::mat GG){
    // env should be column vector
    return trans(env) * GG * env;
}

//' Compute a  white noise environment with a shift
//' @param env previous year's environment NOT used in current formulation
//' @param t the time point
//' @param rho the correlation between environment of selection and the cue
//' @param env.args other args:
//' t.jump the location for the jump,
//' delta env change that takes place at t.jump
//' sd noise in env
//' sdc noise in cue
//' @details nada
//' @import MASS
//' @export
// [[Rcpp::export]]
arma::colvec Env_shift(int t, double rho, arma::mat sigma, double t_jump, double k) {
  // Rcpp::List rparam(env_args); implicit? 
  //  double t_jump = as<double>(env_args["t.jump"]); // no Rcpp:: needed?
  //  double k = as<double>(env_args["k"]);
  arma::colvec out;
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




