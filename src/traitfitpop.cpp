#include "RcppArmadillo.h" // after version 4.3
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; // done by sourceCpp?


// ** below at least compile when used in 
// cppFunction(code="//code here", depends="RcppArmadillo")

//' Compute mean fitness under stabilizing selection
//' @param zbar mean trait prior to selection
//' @param theta enironmental optimum
//' @param Oz2 strength of stabilizing selection
//' @param gamma 1/(Oz2 + Vz2), with Vz2 phenotypic variance 
//' @details no details
//' @export
// [[Rcpp::export]]
double Log_W_bar (double zbar, double theta, double Oz2, double gamma){
    return 0.5 * log(gamma * Oz2) - gamma /2 * pow(zbar - theta, 2); //## compare with eqn 2c lande 2009
}

//' Compute mean fitness under stabilizing selection
//' @inheritParams Log_W_bar
//' @param LOG (=TRUE) whether to return the log of fitness  
//' @details no details
//' @export
// [[Rcpp::export]]
double W_bar (double zbar, double theta, double Oz2, double gamma, bool LOG){
  if(LOG){
    return Log_W_bar(zbar, theta, Oz2, gamma); 
      }
  else {
    return  sqrt(gamma * Oz2) * exp(- gamma /2 * pow(zbar - theta,2)); 
  }
}

//' Compute population growth rate under stabilizing selection
//' @param R0 basic reproductive number
//' @param Wbar average fitness
//' @details Assumes density-independent
//' @export
// [[Rcpp::export]]
double  R_bar (double R0, double Wbar){
    return R0 * Wbar;
      }

//' Compute LOG population growth rate under stabilizing selection
//' @inheritParams R_bar
//' @param logWbar log average fitness
//' @details Assumes theta-logistic population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  Log_R_bar (double R0, double logWbar){
  return log(R0) + logWbar;
      }

//' Compute population growth rate under stabilizing selection
//' @inheritParams R_bar
//' @param N number of individuals in this generation
//' @param K carrying capacity
//' @param thetaL theta-logistic parameter for density dependence
//' @details Assumes theta-logistic population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  R_bar_dd(double R0, double Wbar, double N, double K, double thetaL){
    return pow(R0, (1 -pow(N/K, thetaL))) * Wbar;
      }

//' Selection on plasticity as as function of environment assuming stabilizing selection
//' @param gamma = 1/(Oz2 + Vz2), with Vz2 phenotypic variance
//' @param A the environmental optimum RN int
//' @param B the environmental optimum RN slope
//' @param a current value of RN int
//' @param b current value of RN slope
//' @param e_t environment now
//' @param e_plast env that cues plasticity
//' @export
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

//' Va additive genetic variance in the phenotype as a function of the environment
//' @param env environment in which plastiicty is cued
//' @param GG additive genetic variance-covariance matrix
//' @export
// [[Rcpp::export]]
double Va (arma::vec env, arma::mat GG){
    // env should be column vector
  return as_scalar(trans(env) * GG * env);
}

//' Compute a white noise environment with a shift
//' @param t the time point
//' @param env_args other args:
//' t.jump the location for the jump,
//' delta env change that takes place at t.jump
//' sd noise in env
//' sdc noise in cue
//' @details nada
//' @export
// [[Rcpp::export]]
arma::vec Env_shift_cpp(int t, List env_args) {
  // returns (selection env, cue env)
  double t_jump = as<double>(env_args["t.jump"]);
  arma::mat sigma = as<arma::mat>(env_args["sigma"]);
  //double rho  = as<double>(env_args["rho"]);
  //double k = as<double>(env_args["k"]);
  arma::vec out(2); // maybe make a row vec to avoid trans() in output below?

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
  return out + arma::trans(Z * arma::chol(sigma)); // wld use arma::repmat(out, 1, n) for n>1
  // for red noise do somehting like k * env  + sqrt(1 - k^2) * rnorm (n=1, mean=0, sd=sd)) but hard
  // to figure this out with the cue / correlation thing

}


//' Compute phenotypic dynamic Time Series of trait + demographic change under stabilizing selection as 
//' function of environment (after Lande Chevin)
//' @param T end time, assuming start time of 1
//' @param X parameters (z, a, b, wbar, logN, theta)
//' @param params a list with (gamma_sh, omegaz, A, B, R0, Va, Vb, delta, sigma_xi, rho_tau, fractgen)
//' @param env_args extra args for env.fn
//' @details NB - for now assume Tchange = 0 and demography after CL 2010
//' @export
// [[Rcpp::export]]

arma::mat pdTS (int T, arma::rowvec X, List params, 
		       List env_args) {
  // TODO - make above functions consistent with this approach
  // NB - for now assume Tchange = 0
  // state vars in X (R indexing)
  // zbar <- X[1]
  // abar <- X[2] // elevation
  // bbar <- X[3]  // slope
  // Wbar <- X[4]
  // Npop <- X[5]
  // Theta <- X[6]

  //tmp vars for bio
  double zbart;
  double abart; 
  double bbart;
  double wbart;
  double thetat;
  double betat;


  double gamma_sh = as<double>(params["gamma_sh"]);
  double omegaz = as<double>(params["omegaz"]);
  double A =  as<double>(params["A"]);
  double B = as<double>(params["B"]);
  double R0 = as<double>(params["R0"]);
  double Va = as<double>(params["Va"]);
  double Vb = as<double>(params["Vb"]);

  double delta = as<double>(env_args["delta"]);
  double sigma_xi = as<double>(env_args["sigma_xi"]);
  double rho_tau = as<double>(env_args["rho_tau"]);
  double fractgen = as<double>(env_args["fractgen"]);

  double gamma = 0;
  gamma  = gamma_sh; 
  
  int len = X.size();
  arma::mat out(len, T + 1);
  out.col(0) = trans(X); 

  // tmp vars for env
  double xi_dev_temp; 
  double xi_sel_temp; 
  double eps_sel; 
  double eps_dev; 

   xi_sel_temp = 0;

  arma::vec all_env = arma::randn((T + 1) * fractgen);    
  int env_ctr = 0; 

  all_env = all_env * sigma_xi;  // scale noise

  for(int t=0; t <= T; t++) {

    if (t == 0)
      xi_dev_temp = 0;
    else
      xi_dev_temp = xi_sel_temp;

    for(int i=1; i < fractgen; i++) {
      xi_dev_temp = rho_tau * xi_dev_temp + sqrt(1 - pow(rho_tau, 2)) * all_env(env_ctr); // make drawenv
      env_ctr = env_ctr + 1; 
    }
    

    xi_sel_temp =  rho_tau * xi_dev_temp + sqrt(1 - pow(rho_tau, 2)) * all_env(env_ctr); // make drawenv
    env_ctr = env_ctr + 1; 

    //##    computing the corresponding mean/variance genetic values, and optimal environmen  
    eps_dev =  delta + xi_dev_temp; // ##xi_temp = xi_dev
    eps_sel =  delta +  xi_sel_temp; 
    
    abart = out(1, t);
    bbart = out(2, t);
    zbart =     abart + bbart * eps_dev; // zbar = abar + bbar * eps_dev
    out(0, t) = zbart;
    //    varz <- Va + 2 * Vab * eps_dev + Vb * pow(eps_dev, 2) + Ve;
    thetat = A + B * eps_sel; //Theta
    out(5, t) = thetat; 
    //        ## evolutionary change with and population growth
    wbart = R0 *  sqrt(gamma * pow(omegaz, 2)) * exp(-gamma / 2 *  pow(zbart - thetat, 2)); 
    out(3, t) = wbart;
    // Wbar= f( zbar - theta)
    
    //##  in R version increment here to update the vars with init conditions already filled in
    //    t = t + 1;

    if(t < T) {
      betat = gamma * (zbart - thetat);
      // abar = abar - gamma(zbar - theta) * Va
      out(1, t + 1) = abart - betat * Va; 
      // bbar = bbar - gamma (zbar - theta) eps_dev Vb
      out(2, t + 1) =    bbart - betat * eps_dev * Vb ; 
      out(4, t + 1) = out(4, t) * wbart; // N = N * wbar
    }
  }
  return out;
}
