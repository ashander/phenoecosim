// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Log_W_bar
double Log_W_bar(double zbar, double theta, double oz2, double gamma);
RcppExport SEXP phenoecosim_Log_W_bar(SEXP zbarSEXP, SEXP thetaSEXP, SEXP oz2SEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type zbar(zbarSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type oz2(oz2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    __result = Rcpp::wrap(Log_W_bar(zbar, theta, oz2, gamma));
    return __result;
END_RCPP
}
// W_bar
double W_bar(double zbar, double theta, double oz2, double gamma, bool LOG);
RcppExport SEXP phenoecosim_W_bar(SEXP zbarSEXP, SEXP thetaSEXP, SEXP oz2SEXP, SEXP gammaSEXP, SEXP LOGSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type zbar(zbarSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type oz2(oz2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type LOG(LOGSEXP);
    __result = Rcpp::wrap(W_bar(zbar, theta, oz2, gamma, LOG));
    return __result;
END_RCPP
}
// R_bar
double R_bar(double R0, double Wbar, double N);
RcppExport SEXP phenoecosim_R_bar(SEXP R0SEXP, SEXP WbarSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type Wbar(WbarSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    __result = Rcpp::wrap(R_bar(R0, Wbar, N));
    return __result;
END_RCPP
}
// R_bar_ceiling
double R_bar_ceiling(double R0, double Wbar, double N, double K);
RcppExport SEXP phenoecosim_R_bar_ceiling(SEXP R0SEXP, SEXP WbarSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type Wbar(WbarSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    __result = Rcpp::wrap(R_bar_ceiling(R0, Wbar, N, K));
    return __result;
END_RCPP
}
// R_bar_thetalog
double R_bar_thetalog(double R0, double Wbar, double N, double K0, double thetaL);
RcppExport SEXP phenoecosim_R_bar_thetalog(SEXP R0SEXP, SEXP WbarSEXP, SEXP NSEXP, SEXP K0SEXP, SEXP thetaLSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type Wbar(WbarSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< double >::type thetaL(thetaLSEXP);
    __result = Rcpp::wrap(R_bar_thetalog(R0, Wbar, N, K0, thetaL));
    return __result;
END_RCPP
}
// R_bar_gompertz
double R_bar_gompertz(double R0, double Wbar, double N, double K0);
RcppExport SEXP phenoecosim_R_bar_gompertz(SEXP R0SEXP, SEXP WbarSEXP, SEXP NSEXP, SEXP K0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type Wbar(WbarSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type K0(K0SEXP);
    __result = Rcpp::wrap(R_bar_gompertz(R0, Wbar, N, K0));
    return __result;
END_RCPP
}
// Va
double Va(arma::vec env, arma::mat GG);
RcppExport SEXP phenoecosim_Va(SEXP envSEXP, SEXP GGSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type env(envSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type GG(GGSEXP);
    __result = Rcpp::wrap(Va(env, GG));
    return __result;
END_RCPP
}
// N_e
double N_e(double R0, double N);
RcppExport SEXP phenoecosim_N_e(SEXP R0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    __result = Rcpp::wrap(N_e(R0, N));
    return __result;
END_RCPP
}
// SHC
double SHC(double sigma_g2, double omegaz2, double N, double R0);
RcppExport SEXP phenoecosim_SHC(SEXP sigma_g2SEXP, SEXP omegaz2SEXP, SEXP NSEXP, SEXP R0SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type sigma_g2(sigma_g2SEXP);
    Rcpp::traits::input_parameter< double >::type omegaz2(omegaz2SEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type R0(R0SEXP);
    __result = Rcpp::wrap(SHC(sigma_g2, omegaz2, N, R0));
    return __result;
END_RCPP
}
// Env_shift_cpp
arma::vec Env_shift_cpp(int t, List env_args);
RcppExport SEXP phenoecosim_Env_shift_cpp(SEXP tSEXP, SEXP env_argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type env_args(env_argsSEXP);
    __result = Rcpp::wrap(Env_shift_cpp(t, env_args));
    return __result;
END_RCPP
}
// make_env
arma::mat make_env(int T, List env_args);
RcppExport SEXP phenoecosim_make_env(SEXP TSEXP, SEXP env_argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< List >::type env_args(env_argsSEXP);
    __result = Rcpp::wrap(make_env(T, env_args));
    return __result;
END_RCPP
}
// simulate_pheno_ts
arma::mat simulate_pheno_ts(int T, arma::rowvec X, List params, List env_args, std::string growth_fun, bool poisson, bool varying_g);
RcppExport SEXP phenoecosim_simulate_pheno_ts(SEXP TSEXP, SEXP XSEXP, SEXP paramsSEXP, SEXP env_argsSEXP, SEXP growth_funSEXP, SEXP poissonSEXP, SEXP varying_gSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type env_args(env_argsSEXP);
    Rcpp::traits::input_parameter< std::string >::type growth_fun(growth_funSEXP);
    Rcpp::traits::input_parameter< bool >::type poisson(poissonSEXP);
    Rcpp::traits::input_parameter< bool >::type varying_g(varying_gSEXP);
    __result = Rcpp::wrap(simulate_pheno_ts(T, X, params, env_args, growth_fun, poisson, varying_g));
    return __result;
END_RCPP
}
