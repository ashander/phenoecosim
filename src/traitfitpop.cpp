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

//' Compute LOG population growth rate under stabilizing selection
//' @inheritParams R_bar_dd
//' @param logWbar log average fitness
//' @details Assumes theta-logistic population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  Log_R_bar_dd(double R0, double logWbar, double N, double K, double thetaL){
  return log(R0) *(1 -pow(N/K, thetaL)) + logWbar;
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

// Internal function
arma::rowvec pdIter(int len, double z_bar, double theta_t, double Oz2, double gamma, double A, double B, 
		    double a, double b, double e_t, double e_plast, 
		    arma::mat G,
		    double logN, double R0, double K, double thetaL){
  double log_w_bar = Log_W_bar(z_bar, theta_t , Oz2, gamma); 

  arma::vec beta = Beta(gamma, A, B, a, b, e_t, e_plast);

  // construct the output
  arma::rowvec out(len);
  out.cols(0,1) = trans(G * beta);
  out(0) += a;
  out(1) += b;
  out(2) = e_t;
  out(3) = Log_R_bar(R0, log_w_bar)  + logN; // demo update based on mean fitness
    //Log_R_bar_dd( R0, log_w_bar, exp(logN), K, thetaL)  + logN; 
  return out;
}


//' Compute change in fitness over one generation after stabilizing selection
//' function of environment (after Lande Chevin)
//' @param t timestep
//' @param X parameters (a, b, env, logN)
//' @param params a list with (Oz2, AB, R0, K, theta)
//' @param G the (constant) G matrix
//' @param env_args extra args for env.fn
//' @details function signature as for use with deSolve
//' imposes fitness based on Lande
//' and demography after CL 2010
//' @return rowvector: log R bar, delta R bar, variance load
//' @export
// [[Rcpp::export]]
arma::rowvec pdRespVarLoad (int t, arma::rowvec X, List params, 
		   arma::mat G, List env_args) {
  // state vars in X (R indexing)
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  //// env <- X[3]
  // logN = X[4]

  double Oz2  = as<double>(params["Oz2"]);
  arma::rowvec AB = as<arma::rowvec>(params["AB"]);// or set size first? AB(2);
  double R0 = as<double>(params["R0"]);
  double K = as<double>(params["K"]);
  double thetaL = as<double>(params["thetaL"]);
  double Ve = as<double>(params["Ve"]);
  
  arma::vec ee = Env_shift_cpp(t, env_args);
  arma::vec env(2);
  env.fill(1.0);
  env(1)= ee(1);
  
  double z_bar = arma::as_scalar(X.cols(0,1) * env); 
  double theta_t = arma::as_scalar(AB * env);
  double gamma = 1 / (Oz2 + Va(env, G) + Ve);
  arma::rowvec out(3);
  out(0) = Log_R_bar( R0, Log_W_bar(z_bar, theta_t, Oz2, gamma)); // mean fitness pre selection
  //Log_R_bar_dd( R0, Log_W_bar(z_bar, theta_t, Oz2, gamma), exp(X(3)),K, thetaL); 

  // selection
  arma::rowvec oneiter =  pdIter(X.size(), z_bar, theta_t, Oz2, gamma, AB(0), AB(1), 
		X(0), X(1), ee(0), ee(1),
		G,
		X(3), R0, K, thetaL);
  z_bar = as_scalar(oneiter.cols(0,1) * env); // update mean z based on selection

  out(1) = Log_R_bar(R0, Log_W_bar(z_bar, theta_t, Oz2, gamma)) - out(0);
  //Log_R_bar_dd( R0, Log_W_bar(z_bar, theta_t, Oz2, gamma), exp(oneiter(3)),K, thetaL) - out(0);
  out(2) = 0.5 * log( 1 + (Va(env, G) + Ve) / Oz2);
  return out;
}


//' Compute trait + demographic change over one generation under stabilizing selection as 
//' function of environment (after Lande Chevin)
//' @param t timestep
//' @param X parameters (a, b, env, logN)
//' @param params a list with (Oz2, AB, R0, K, theta)
//' @param G the (constant) G matrix
//' @param env_args extra args for env.fn
//' @details function signature as for use with deSolve
//' imposes fitness based on Lande
//' and demography after CL 2010
//' @export
// [[Rcpp::export]]
arma::rowvec pdLande (int t, arma::rowvec X, List params, 
		   arma::mat G, List env_args) {
  // state vars in X (R indexing)
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  //// env <- X[3]
  // logN = X[4]

  double Oz2  = as<double>(params["Oz2"]);
  arma::rowvec AB = as<arma::rowvec>(params["AB"]);// or set size first? AB(2);
  double R0 = as<double>(params["R0"]);
  double K = as<double>(params["K"]);
  double thetaL = as<double>(params["thetaL"]);
  double Ve = as<double>(params["Ve"]);

  arma::vec ee = Env_shift_cpp(t, env_args);
  arma::vec env(2);
  env.fill(1.0);
  env(1)= ee(1);
  
  double z_bar = arma::as_scalar(X.cols(0,1) * env); 
  double theta_t = arma::as_scalar(AB * env);
  double gamma = 1 / (Oz2 + Va(env, G) + Ve);
  return pdIter(X.size(), z_bar, theta_t, Oz2, gamma, AB(0), AB(1), 
		X(0), X(1), ee(0), ee(1),
		G,
		X(3), R0, K, thetaL);
}


//' Compute trait + demographic change over T generations under stabilizing selection as 
//' function of environment (after Lande Chevin)
//' @param T end time, assuming start time of 1
//' @param X parameters (a, b, env, logN)
//' @param params a list with (Oz2, AB, R0, K, theta)
//' @param G the (constant) G matrix
//' @param env_args extra args for env.fn
//' @details function signature as for use with deSolve
//' imposes fitness based on Lande
//' and demography after CL 2010
//' @export
// [[Rcpp::export]]
arma::rowvec pdLandeT (int T, arma::rowvec X, List params, 
		   arma::mat G, List env_args) {
  // state vars in X (R indexing)
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  //// env <- X[3]
  // logN = X[4]

  double Oz2  = as<double>(params["Oz2"]);
  arma::rowvec AB = as<arma::rowvec>(params["AB"]);// or set size first? AB(2);
  double R0 = as<double>(params["R0"]);
  double K = as<double>(params["K"]);
  double thetaL = as<double>(params["thetaL"]);
  double Ve = as<double>(params["Ve"]);

  double gamma = 0;
  int len = X.size();
  arma::rowvec out(len);
  out = X; 
  for(int t=1; t <= T; t++){
    arma::vec ee = Env_shift_cpp(t, env_args);
    arma::vec env(2);
    env.fill(1.0);
    env(1)= ee(1);

    double z_bar = arma::as_scalar(X.cols(0,1) * env); 
    double theta_t = arma::as_scalar(AB * env);
    gamma = 1 / (Oz2 + Va(env, G) + Ve);

    out = pdIter(len, z_bar, theta_t, Oz2, gamma, AB(0), AB(1), 
		  out(0), out(1), ee(0), ee(1),
		  G,
		  out(3), R0, K, thetaL);
  }
  return out;
}


//' Compute trait + demographic change over over T generations and
//' long run growth rate over last N_lam generations
//' under stabilizing selection as 
//' function of environment (after Lande Chevin)
//' @param T end time, assuming start time of 1
//' @param N_lam number of gens over which to compute long-run growth  -- must be less than T
//' @param X parameters (a, b, env, logN)
//' @param params a list with (Oz2, AB, R0, K, theta)
//' @param G the (constant) G matrix
//' @param env_args extra args for env.fn
//' @details function signature as for use with deSolve
//' imposes fitness based on Lande
//' and demography after CL 2010
//' @export
// [[Rcpp::export]]
double pdLandeLRG (int T, int N_lam, arma::rowvec X, List params, 
		   arma::mat G, List env_args) {
  // state vars in X (R indexing)
  // a <- X[1] // elevation
  // b <- X[2]  // slope
  //// env <- X[3]
  // logN = X[4]

  // TODO - error check  N_lam < T
  double Oz2  = as<double>(params["Oz2"]);
  arma::rowvec AB = as<arma::rowvec>(params["AB"]);// or set size first? AB(2);
  double R0 = as<double>(params["R0"]);
  double K = as<double>(params["K"]);
  double thetaL = as<double>(params["thetaL"]);
  double Ve = as<double>(params["Ve"]);

  double gamma = 0;
  int len = X.size();
  arma::rowvec out(len);
  out = X; 
  for(int t=1; t <= T - N_lam; t++){
    arma::vec ee = Env_shift_cpp(t, env_args);
    arma::vec env(2);
    env.fill(1.0);
    env(1)= ee(1);

    double z_bar = arma::as_scalar(X.cols(0,1) * env); 
    double theta_t = arma::as_scalar(AB * env);
    gamma = 1 / (Oz2 + Va(env, G) + Ve);

    out = pdIter(len, z_bar, theta_t, Oz2, gamma, AB(0), AB(1), 
		  out(0), out(1), ee(0), ee(1),
		  G,
		  out(3), R0, K, thetaL);
  }



  double logN_0 = out(3); // 0 index  -- store logN to calculate lrg
  for(int t=T - N_lam + 1; t <= T; t++){
    arma::vec ee = Env_shift_cpp(t, env_args);
    arma::vec env(2);
    env.fill(1.0);
    env(1)= ee(1);

    double z_bar = arma::as_scalar(X.cols(0,1) * env); 
    double theta_t = arma::as_scalar(AB * env);
    gamma = 1 / (Oz2 + Va(env, G) + Ve);

    out = pdIter(len, z_bar, theta_t, Oz2, gamma, AB(0), AB(1), 
		  out(0), out(1), ee(0), ee(1),
		  G,
		  out(3), R0, K, thetaL);
  }

  
  return (out(3) - logN_0) / N_lam ;  // calculate long run growth rate (logN(T) - logN(0)) / T
}


//' Compute phenotypic dynamic Time Series of trait + demographic change under stabilizing selection as 
//' function of environment (after Lande Chevin)
//' @param T end time, assuming start time of 1
//' @param X parameters (z, a, b, wbar, logN, theta)
//' @param params a list with (Oz2, AB, R0, K, theta)
//' @param env_args extra args for env.fn
//' @details NB - for now assume Tchange = 0
//' and demography after CL 2010
//' @export
// [[Rcpp::export]]

// TODO - make above functions consistent with this approach
arma::mat pdTS (int T, arma::rowvec X, List params, 
		       List env_args) {
  // NB - for now assume Tchange = 0
  // state vars in X (R indexing)
  // zbar <- X[1]
  // abar <- X[2] // elevation
  // bbar <- X[3]  // slope
  // Wbar <- X[4]
  // Npop <- X[5]
  // Theta <- X[6]

  //  double Vz = as<double>(params["Vz"]);
  double gamma_sh = as<double>(params["gamma_sh"]);
  //  double rmax = as<double>(params["rmax"]);
  double omegaz = as<double>(params["omegaz"]);
  //  arma::rowvec AB = as<arma::rowvec>(params["AB"]);// or set size first? AB(2);
  double A =  as<double>(params["A"]);
  double B = as<double>(params["B"]);
  double R0 = as<double>(params["R0"]);
  double Va = as<double>(params["Va"]);
  double Vb = as<double>(params["Vb"]);
  //  double Ve = as<double>(params["Ve"]); // only to compute   double varz; 

  double delta = as<double>(env_args["delta"]);
  double sigma_xi = as<double>(env_args["sigma_xi"]);
  double rho_tau = as<double>(env_args["rho_tau"]);
  double fractgen = as<double>(env_args["fractgen"]);

  double gamma = 0;
  gamma  = gamma_sh; 


  int len = X.size();
  arma::mat out(len, T + 1);
  out.col(0) = trans(X); 

  // tmp vars
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
    
    
    out(0 , t) =     out(1, t) + out(2, t) * eps_dev; // zbar = abar + bbar * eps_dev
    //    varz <- Va + 2 * Vab * eps_dev + Vb * pow(eps_dev, 2) + Ve;
    out(5, t) = A + B * eps_sel; //Theta
    //        ## evolutionary change with and population growth
    out(3, t) = R0 *  sqrt(gamma * pow(omegaz, 2)) * exp(-gamma / 2 *  pow(out(0, t) - out(5, t), 2)); // Wbar= f( zbar - theta)
    
    //##  in R version increment here to update the vars with init conditions already filled in
    //    t = t + 1;

    if(t < T) {
      
      out(1, t + 1) = out(1, t) - gamma * (out(0, t) - out(5, t)) * Va; // abar = abar - gamma(zbar - theta) * Va
      out(2, t + 1) =    out(2, t) - gamma * (out(0, t) - out(5, t)) * eps_dev * Vb ; // bbar = bbar - gamma (zbar - theta) eps_dev Vb
      out(4, t + 1) = out(4, t) * out(3, t); // N = N * wbar
    }
  }
  return out;
}
