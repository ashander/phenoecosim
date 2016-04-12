
library(phenoecosim)
context("Simulation testing")

test_envr <- list(nrep=1,
                  delta=2,
                  rho=0.5,
                  alpha=0.3,
                  omegaz=sqrt(25),
                  omega_Wmax=1.05,
                  A=0,
                  B=2,
                  K0=100000000,
                  Tchange=0,
                  sigma_xi=1,
                  Npop0=1e4,
                  var_a=0.1,
                  Vb=0.005,
                  Ve=0.01,
                  fractgen=5,
                  Tstart=0,
                  Nc=50,
                  Tlim=10000)


test_that("Sim ts produce right object", {
    invisible(list2env(test_envr, envir=environment()))
    set.seed(111)
    tmp <- simulate_timeseries(nrep=nrep, Tlim=Tlim, delta=delta, rho=rho, alpha=alpha, omegaz=omegaz,
                  sigma_xi=sigma_xi, Npop=Npop0, var_a=var_a, Vb=Vb,
                  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                  A=A, B=B, K0=K0, stationarity = FALSE)
    expect_equal_to_reference(tmp, 'sim-ts.rds', scale=1, tolerance=0.1)
})

# ITS A DATA FRAME
test_that("run_sims_on_grid traj produce right object", {
    invisible(list2env(test_envr, envir=environment()))
    set.seed(111)

    alpha <- 0.5
    omegaz <- sqrt(20)
    traj.data <- data.frame(
        rho = c(0.3, 0.7, 0.7),
        delta=c(2.5, 2.5, 5)
        )
    traj.data$obs <- 1:nrow(traj.data)
    traj.data$omegaz <- omegaz
    traj.data$Vb  <- Vb
    traj.data$alpha <- alpha
    Tlim <-1500
    example_traj_reps <- 10
    time_calc <- system.time(traj.data <- run_sims_on_grid(summary_id, traj.data,  nrep=example_traj_reps, var_a=var_a,
                           Tlim=Tlim, seed=111,  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                           Npop0=Npop0, sigma_xi=sigma_xi, K0=K0, B=B, A=A, stationarity = FALSE))
    expect_equal_to_reference(traj.data, "calc_many_traj.rds", scale=1, tolerance=0.1)
    cat("\nTime traj calc: ", time_calc['elapsed'], "\n")
    expect_less_than(time_calc['elapsed'], 0.6)
})

test_that("run_sims_on_grid traj speed", {
    invisible(list2env(test_envr, envir=environment()))
    set.seed(111)

    alpha <- 0.5
    omegaz <- sqrt(20)
    traj.data <- data.frame(
        rho = c(0.3, 0.7, 0.7),
        delta=c(2.5, 2.5, 5)
        )
    traj.data$obs <- 1:nrow(traj.data)
    traj.data$omegaz <- omegaz
    traj.data$Vb  <- Vb
    traj.data$alpha <- alpha
    time_calc <- system.time(run_sims_on_grid(summary_id, traj.data,  nrep=20, var_a=var_a,
                           Tlim=200, seed=111,  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                           Npop0=Npop0, sigma_xi=sigma_xi, K0=K0, B=B, A=A, stationarity = FALSE))
    cat("\nTime traj calc: ", time_calc['elapsed'], "\n")
    expect_less_than(time_calc['elapsed'], 0.6)

})


################################################################################
test_that("Stationarity checks", {
    invisible(list2env(test_envr, envir=environment()))

  set.seed(111)
  Vz <- Vz_(var_a, Vb, delta, Ve);
  gamma_sh <-  1 / (omegaz ^ 2 + Vz) # gamma after the environmental shift
  rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
  abar0 <- A
  bbar0 <- alpha * B
  param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                     A=A, B=B, R0=omega_Wmax, var_a=var_a, Vb=Vb, Ve=Ve)
  env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho,
                   fractgen=fractgen)

  nrep <- 100
  X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
  X0_list <- stationary_draws(nrep, X0, param.list, env.list, alpha)
  tmp <- do.call(rbind, X0_list)
  tm <- colMeans(tmp)
  expect_equal(unname(tm['abar']), 0.0, tolerance = 0.1, scale =1)
  expect_equal(unname(tm['bbar'])/B, alpha, tolerance = 0.1, scale =1)
  expect_equal_to_reference(X0_list, 'stationarity.rds', scale=1, tolerance=0.1)
})


test_that("Sim ts produce right object with stationarity", {
    invisible(list2env(test_envr, envir=environment()))

    set.seed(111)
    tmp <- simulate_timeseries(nrep=nrep, Tlim=Tlim, delta=delta, rho=rho, alpha=alpha, omegaz=omegaz,
                  sigma_xi=sigma_xi, Npop=Npop0, var_a=var_a, Vb=Vb,
                  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                  A=A, B=B, K0=K0, stationarity = TRUE)
    expect_equal_to_reference(tmp, 'sim-ts-stat.rds', scale=1, tolerance=0.1)
})


test_that("run_sims_on_grid traj produce right object with stationarity", {
    invisible(list2env(test_envr, envir=environment()))
    set.seed(111)

    alpha <- 0.5
    omegaz <- sqrt(20)
    traj.data <- data.frame(
        rho = c(0.3, 0.7, 0.7),
        delta=c(2.5, 2.5, 5)
        )
    traj.data$obs <- 1:nrow(traj.data)
    traj.data$omegaz <- omegaz
    traj.data$Vb  <- Vb
    traj.data$alpha <- alpha
    Tlim <-1500
    example_traj_reps <- 10
    time_calc <- system.time(traj.data <- run_sims_on_grid(summary_id, traj.data,  nrep=example_traj_reps, var_a=var_a,
                           Tlim=Tlim, seed=111,  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                           Npop0=Npop0, sigma_xi=sigma_xi, K0=K0, B=B, A=A, stationarity = TRUE))
    expect_equal_to_reference(traj.data, "calc_many_traj-stat.rds", scale=1, tolerance=0.1)
    cat("\nTime persist calc w/ stationarity: ", time_calc['elapsed'], "\n")
    expect_less_than(time_calc['elapsed'], 0.6)
})

test_that("run_sims_on_grid traj speed with stationarity", {
    invisible(list2env(test_envr, envir=environment()))
    set.seed(111)

    alpha <- 0.5
    omegaz <- sqrt(20)
    traj.data <- data.frame(
        rho = c(0.3, 0.7, 0.7),
        delta=c(2.5, 2.5, 5)
        )
    traj.data$obs <- 1:nrow(traj.data)
    traj.data$omegaz <- omegaz
    traj.data$Vb  <- Vb
    traj.data$alpha <- alpha
    time_calc <- system.time(run_sims_on_grid(summary_id, traj.data,  nrep=20, var_a=var_a,
                           Tlim=200, seed=111,  Ve=Ve, fractgen=fractgen, omega_Wmax=omega_Wmax,
                           Npop0=Npop0, sigma_xi=sigma_xi, K0=K0, B=B, A=A, stationarity = TRUE))
    cat("\nTime persist calc w/ stationarity: ", time_calc['elapsed'], "\n")
    expect_less_than(time_calc['elapsed'], 0.6)
})
