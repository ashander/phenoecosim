#' run a simulation on a parameter grid and summarize it
#' @param summary_function a function that returns a data frame based on
#' simulation output
#' @param d a data.frame that defines the grid of parameters to run the
#' simulation on
#' @param Tlim the length of time to run
#' @param ... should include params for simulate_timeseries
#' @details
#' d - data frame with columns of rho, delta, omegaz, alpha, Vb
#' value a data.frame with possibly many reps per parameter or some summary
#' statistic. actual value depends on return value of summary_function.
#' built-in summary functions include the identity (no summary) `summary_id`
#' @export
run_sims_on_grid <- function(summary_function,  d,  Tlim, ...) {
  if (!is.function(summary_function))
    stop("must be function")
  if (is.na(Tlim))
    stop("NA in time domain")

  as.data.frame(runner(d, summary_function, Tlim=Tlim, ...))
}

#' @importFrom dplyr do_ rowwise
runner <- function(d, summarizer, ...){
  do_(rowwise(d), ~ data.frame(. , do.call(summarizer, c(., ...))))
}

#' unsummarized raw data
#' @param ... parameters used by simulate timeseries
#' @export
summary_id <- function( ...) {
  ## can do something like function(Npop0, ...)
  ## do stuff with Npop0 or other args
  ## then must call
  ## simulate_timeseries(Npop0=Npop0, ...) # argslist includes all args from above
  simulate_timeseries(...)
}

#' simulate a timeseries
#' @param nrep reps
#' @param Tlim max time
#' @param delta shift
#' @param rho predictability of the environmental noise
#' @param alpha intial mean plasticity
#' @param omegaz selection
#' @param Vb variance in RN slope
#' @param sigma_xi variance in the optimum
#' @param Npop0 initial population
#' @param var_a variance in trait
#' @param Ve microenvronemntal variance
#' @param fractgen developmental delay
#' @param omega_Wmax max productivity
#' @param A reference environment intercept
#' @param B environmental sensitivity of selection
#' @param K0 carrying capacity (default=10000)
#' @param stationarity (logical) generate stationary IC if TRUE
#' @param density density dependence ('independent', 'gompertz', 'thetalogistic', 'ceiling')
#' @param thetalog theta parameter for theta-logistic
#' @param poisson (logical) reproduction occurs as N(t+1) = Poisson(f(N(t)))
#' @param ... avoid throwing errors if we pass too many? (bad idea)
#' @details  uses `simulate_pheno_ts` under the hood. Poisson only implemented
#' for density='independent'
#' @export
simulate_timeseries <- function(nrep, Tlim, delta, rho, alpha, omegaz, Vb,
                                sigma_xi, Npop0, var_a,  Ve=Ve,
                                fractgen=fractgen, omega_Wmax=omega_Wmax, A=A,
                                B=B, stationarity = TRUE, K0=10000, density="independent",
                                thetalog = NA, poisson=FALSE, ...) {
  ## assumes 0 is the time at which the environment shifts
  ## this because the Vz is set at outset based on delta
  Vz <- Vz_(var_a, Vb, delta, Ve);
  gamma_sh <-  1 / (omegaz ^ 2 + Vz) # gamma after the environmental shift

  ## gamma in env 0 (or neglecting the effect of variance in plasticity)
  ## would be useful if T_change =/= 0
  ## gamma_0 <-   1 / (omegaz ^ 2 + var_a + Ve);
  ## r discounted by standing load after env shift
  rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))

  ## Initial conditions

  abar0 <- A
  bbar0 <- alpha * B

  ## cpp sim init
  if (density == "independent")
    param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                      A=A, B=B, R0=omega_Wmax, var_a=var_a, Vb=Vb, Ve=Ve)
  else if (density != "thetalogistic")
    param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                      A=A, B=B, R0=omega_Wmax, var_a=var_a, Vb=Vb, Ve=Ve, K0=K0)
  else if (density == "thetalogistic")
    param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                      A=A, B=B, R0=omega_Wmax, var_a=var_a, Vb=Vb, Ve=Ve, K0=K0, 
                      thetalog=thetalog)

  env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho, fractgen=fractgen)
  X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
  ic_generator <- initial_state_generator(X0, stationarity, nrep, param.list, env.list, alpha)
  all.out <- lapply(1:nrep, function(r)
                    {
                      cbind(rep(r, Tlim + 1), 0:Tlim,
                            simulate_pheno_ts(Tlim, ic_generator(r), param.list, env.list,
                                              density, poisson))
                    })
  all.out <- do.call(rbind, all.out)
  colnames(all.out) <- c("rep", "time", names(X0))
  return (all.out)
}

## wrapper, to generate states.
##
## provides a closure that
## (1) either just always returns the state X0
## (2) or, if stationarity, generates a set of stationary states and provides
##     them when asked
initial_state_generator <- function(X0, stationarity, nrep, param.list, env.list, alpha) {
  if (!stationarity){
    function(r) return(X0)
  } else {
    X0_list <- stationary_draws(nrep, X0, param.list, env.list, alpha)
    function(r) X0_list[[r]]
  }
}

## run the system from initial state X0
## amend the trait states (the absolute reaction norm elevation and relative
## reaction norm slope)
## with values from the end of the time sieres
## assuming Tmin is long enough to reach stationarity
stationary_draws <- function(nrep, X0, param.list, env.list, alpha, Tmin=15000, Tmax=30000) {
  if (nrep > Tmax - Tmin)
    stop("provide arguments Tmin and Tmax, where Tmax - Tmin > ", nrep)
  env.list[['delta']] <- 0.01
  env.list[['rho_tau']] <- alpha ## assume population initially matches predictability
  tmp <- simulate_pheno_ts(Tmax, X0, param.list, env.list)
  tmp <- tmp[Tmin:Tmax, ]
  colnames(tmp) <- names(X0)
  bbars <- rnorm(n=nrep, mean = mean(tmp[,'bbar']), sd = sd(tmp[,'bbar']))
  abars <- rnorm(n=nrep, mean = mean(tmp[,'abar']), sd = sd(tmp[,'abar']))
  lapply(1:nrep, function(i) { out <- X0
                             out[['abar']] <- abars[i]
                             out[['bbar']] <- bbars[i]
                             out
                           })
}

#' variance load
#' @param siga2 squared genetic variance in reaction norm intercept
#' @param sigb2 squared genetic variance in reaction norm slope
#' @param delta shift in the mean environment
#' @param sige2 squared microenvronmental variance
#' @aliases variance_load
#' @export
Vz_ <- function(siga2, sigb2, delta, sige2)
    siga2 + sigb2 * delta ^ 2 + sige2

#' G matrix
#' @param Gaa additive genetic variance in RN intercept
#' @param Gbb additive genetic variance in RN slope
#' @param Gab additive genetic covariance between RN slope and height
#' @aliases G_matrix
#' @export
G <- function(Gaa, Gbb, Gab)
    matrix( c(Gaa, Gab, Gab, Gbb), nrow=2) #eqn 3a
