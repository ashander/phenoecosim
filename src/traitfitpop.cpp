#include <RcppArmadillo.h> // after version 4.3
#include <math.h> // to suppress linter warnings
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; // done by sourceCpp?

//' Compute mean fitness under stabilizing selection
//' @param zbar mean trait prior to selection
//' @param theta enironmental optimum
//' @param oz2 strength of stabilizing selection
//' @param gamma 1/(oz2 + Vz2), with Vz2 phenotypic variance
//' @details IN USE
//' @export
// [[Rcpp::export]]
double Log_W_bar (double zbar, double theta, double oz2, double gamma){
	return 0.5 * log(gamma * oz2) - gamma / 2 * pow(zbar - theta, 2); //## compare with eqn 2c lande 2009
}

//' Compute mean fitness under stabilizing selection
//' @inheritParams Log_W_bar
//' @param LOG (=TRUE) whether to return the log of fitness
//' @details IN USE
//' @export
// [[Rcpp::export]]
double W_bar (double zbar, double theta, double oz2, double gamma, bool LOG){
	if(LOG){
		return Log_W_bar(zbar, theta, oz2, gamma);
	}
	else {
		return  sqrt(gamma * oz2) * exp(- gamma / 2 * pow(zbar - theta, 2));
	}
}

//' Compute population growth rate under stabilizing selection and no regulation
//' @param R0 basic reproductive number
//' @param Wbar average fitness
//' @param N number of individuals in this generation
//' @details Assumes ceiling population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  R_bar(double R0, double Wbar, double N) {
	return R0 * Wbar * N;
}

// @param extra, extra2 unused parameter so function signature is consistent
double  R_bar_wrap(double R0, double Wbar, double N, double extra, double extra2) {
	return R_bar(R0, Wbar, N);
}

//' Compute population growth rate under stabilizing selection and ceiling regulation
//' @inheritParams R_bar
//' @param K carrying capacity
//' @details Assumes ceiling population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  R_bar_ceiling(double R0, double Wbar, double N, double K) {
	N = R0 * Wbar * N;
	if(N > K) {
		return K;
	} else {
		return N;
	}
}
double  R_bar_ceiling_wrap (double R0, double Wbar, double N, double K, double extra) {
	return R_bar_ceiling(R0, Wbar, N, K);
}

//' Compute population growth rate under stabilizing selection and theta-logistic regulation
//' @inheritParams R_bar
//' @param K0 carrying capacity with optimal trait
//' @param thetaL theta-logistic parameter for density dependence
//' @details Assumes theta-logistic population regulation
//' would be good to have separate DD function specified after
//' chevin and lande 2010
//' @export
// [[Rcpp::export]]
double  R_bar_thetalog(double R0, double Wbar, double N, double K0, double thetaL) {
	return N * pow(R0, 1 - pow(N, thetaL) / pow(K0, thetaL)) * Wbar;
}

//' Compute population growth rate under stabilizing selection and Gompertz regulation
//' @inheritParams R_bar
//' @param K0 carrying capacity at optimum trait
//' @details Assumes Gompertz population regulation
//' @export
// [[Rcpp::export]]
double  R_bar_gompertz (double R0, double Wbar, double N, double K0) {
	return N * pow(R0, 1 - log(N) / log(K0)) * Wbar;
}
double  R_bar_gompertz_wrap (double R0, double Wbar, double N, double K0, double extra) {
	return R_bar_gompertz(R0, Wbar, N, K0);
}

//' Va additive genetic variance in the phenotype as a function of the environment
//' @param env environment in which plastiicty is cued
//' @param GG additive genetic variance-covariance matrix
//' @export
// [[Rcpp::export]]
double Va(arma::vec env, arma::mat GG){
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

//' Compute correctly-scaled environment of development and selection
//' @param env_args environment args
//' @param T how many gens
//' @details given a pre-generated vector of random environments
//' @export
// [[Rcpp::export]]
arma::mat make_env(int T, List env_args) {
	double sigma_xi = as<double>(env_args["sigma_xi"]);
	double rho_tau = as<double>(env_args["rho_tau"]);
	double fractgen = as<double>(env_args["fractgen"]);
	double delta = as<double>(env_args["delta"]);

	//temp vars
	double rho2 = pow(rho_tau, 2);
	double xi_dev;
	double xi_sel = 0;
	arma::vec all_env = arma::randn((T + 1) * fractgen); // random normal for every timestep
	int env_ctr = 0;
	all_env = all_env * sigma_xi;  // scale noise

	arma::mat xi_all(2, T + 1);
	for(int t=0; t <= T; t++) {
		xi_dev = xi_sel;

		for(int i=1; i < fractgen; i++) {
			xi_dev = rho_tau * xi_dev + sqrt(1 - rho2) * all_env(env_ctr);
			env_ctr = env_ctr + 1;
		}
		xi_sel =  rho_tau * xi_dev + sqrt(1 - rho2) * all_env(env_ctr);

		xi_all(0, t) = xi_dev + delta;
		xi_all(1, t) = xi_sel + delta;
	}
	return xi_all;
}


// helpers for simulate_pheno_ts
enum growth_fn_code {
	c_default,
	c_ceiling,
	c_thetalogistic,
	c_gompertz
};

growth_fn_code hash_str (std::string const& in_string) {
	if (in_string ==  "default") return c_default;
	if (in_string ==  "ceiling") return c_ceiling;
	if (in_string ==  "thetalogistic") return c_thetalogistic;
	if (in_string ==  "gompertz") return c_gompertz;
	return c_default;
}

//' Compute phenotypic dynamic Time Series of trait + demographic change under stabilizing selection as
//' function of environment (after Lande Chevin)
//' @param T end time, assuming start time of 1
//' @param X parameters (z, a, b, wbar, logN, theta)
//' @param params a list with (gamma_sh, omegaz, A, B, R0, var_a, Vb, delta, sigma_xi, rho_tau, fractgen)
//' @param env_args extra args for env.fn
//' @param growth_fun should be one of "default" (no dd), "gompertz" or "ceiling"
//' @details NB - for now assume Tchange = 0 and demography after CL 2010
//' @return a long matrix with columns zbar, abar, bbar, Wbar, Npop, theta
//' @export
// [[Rcpp::export]]
arma::mat simulate_pheno_ts(int T, arma::rowvec X, List params,
		List env_args, std::string growth_fun = "default") {
	// TODO - make above functions consistent with this approach
	// NB - for now assume Tchange = 0
	// state vars in X (R indexing)
	// zbar <- X[1]
	// abar <- X[2] // elevation
	// bbar <- X[3]  // slope
	// Wbar <- X[4]
	// Npop <- X[5]
	// Theta <- X[6]
	int len = X.size();
	if(len != 6) {
		throw std::invalid_argument("Initial conditions must be a vector of length 6!");
	}

	//tmp vars for bio
	double zbart;
	double abart;
	double bbart;
	double wbart;
	double thetat;
	double betat;

	double gamma = as<double>(params["gamma_sh"]);
	double omegaz = as<double>(params["omegaz"]);
	double A =  as<double>(params["A"]);
	double B = as<double>(params["B"]);
	double R0 = as<double>(params["R0"]);
	double var_a = as<double>(params["var_a"]);
	double Vb = as<double>(params["Vb"]);
	double dens_dep = 0.0;
	double dens_dep_secondary = 0.0;
	double R0_growth = R0;

	double (*R_function)(double, double, double, double, double);
	switch (hash_str(growth_fun)) {
		case c_default:
			R_function = &R_bar_wrap;
			break;
		case c_ceiling:
			R_function = &R_bar_ceiling_wrap;
			dens_dep = as<double>(params["K0"]);
			break;
		case c_thetalogistic:
			R_function = &R_bar_thetalog;
			dens_dep = as<double>(params["K0"]);
			dens_dep_secondary = as<double>(params["thetalog"]);
			break;
		case c_gompertz:
			R_function = &R_bar_gompertz_wrap;
                        dens_dep = as<double>(params["K0"]);
			break;
		default:
			R_function = &R_bar_wrap;
			break;
	}


	double oz2 = pow(omegaz, 2);

	arma::mat XX(T + 1, len);
	XX.row(0) = X; // col names of XX = row names of X

	// tmp vars for env
	double eps_sel;
	double eps_dev;
	//Rfprintf("Making env .......................")
	arma::mat eps_all = make_env(T, env_args);

	for(int t=0; t <= T; t++) {
		// get environments of dev, sel for this gen
		eps_dev =  eps_all(0, t);
		eps_sel =  eps_all(1, t);

		// compute current-generation values for RN params, environmental optimum, and fitness
		abart = XX( t, 1); // column: reaction norm intercept
		bbart = XX( t, 2); // column: reaction norm slope

		//also save overall trait and optimum -- these depend on
		//environment, which we don't save
		zbart = abart + bbart * eps_dev;
		XX( t, 0) = zbart; //column: overall trait

		thetat = A + B * eps_sel;
		XX( t, 5) = thetat; //column: optimum

		// evolutionary change with and population growth
		wbart = W_bar(zbart, thetat, oz2, gamma, FALSE);
		XX( t, 3) = R0 * wbart; //column: absolute fitness

		if(t < T) {
			// compute next gen values for RN params, population size
			// \delta (a, b) = \beta (var_a, var_b) = gamma ( z - theta) (var_a, var_b)
			// given def betat, bbart line below defines Vb(env) = Vb * env ^2
			//
			// N = N * wbar
			betat = gamma * (zbart - thetat);
			XX( t + 1, 1) = abart - betat * var_a;
			XX( t + 1, 2) = bbart - betat * eps_dev * Vb ;
			XX( t + 1, 4) = R_function(R0, wbart, XX( t, 4), dens_dep, dens_dep_secondary); //column: population size
		}
	}
	return XX;
}
