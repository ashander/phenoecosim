library(phenoecosim)

## variance load
Vz_ <- function(siga2, sigb2, delta, sige2)
    siga2 + sigb2 * delta ^ 2 + sige2

## unvarying parameters
A <- 0
B <- 2 #environmental cline in the optimum phenotype
fractgen <- 5
tau <- 1 / fractgen
Ve <- 0.5
Tchange <- 0 #time for a change in env regime
Npop0 <- 1e4
rho1 <- 0.5
sigma_xi <- 1
omega_z <- sqrt(20)
omega_Wmax <- 1.10 #maximum fitness across environments

context("Phenotypic timeseries regression tests")
for( delta in c(1.5, 5, 8))
  for( Tlim in c(250, 500, 1000))
    for( Va in c(0.05, 0.1))
      for( Vb in c(0.025, 0.05))
        for( omegaz in sqrt(c(10, 25)))
          for( rho0 in c(0.3, 0.7))
            for( alpha in c(0.3, 0.7)) {
              test_that(paste("right object across parameters constant g"), {
                          Vz <- Vz_(Va, Vb, delta, Ve)
                          gamma_sh <-  1 / (omegaz ^ 2 + Vz)
                          rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
                          abar0 <- A;
                          bbar0 <- alpha * B

                          param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                             A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve)

                          env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)
                          X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)

                          set.seed(111)
                          tmp <- simulate_pheno_ts(Tlim, X0, param.list, env.list)
                          target_parms <- paste(Tlim, Va, Vb, omegaz, omega_Wmax, rho0, rho1,
                                                alpha, delta, sigma_xi, sep='-')
                          target_filename <- paste0('simulate', target_parms, '.rds')
                          expect_equal_to_reference(tmp, target_filename, scale=1, tolerance=0.1,
                                                    label=paste('Pars:',
                                                                paste(paste(names(env.list), env.list, sep='='),
                                                                      paste(names(param.list), param.list, sep='='),
                                                                      sep=',', collapse=' ')
                                                                ))

              })
              test_that(paste("right object across parameters constant g, ceiling"), {
                          Vz <- Vz_(Va, Vb, delta, Ve)
                          gamma_sh <-  1 / (omegaz ^ 2 + Vz)
                          rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
                          abar0 <- A;
                          bbar0 <- alpha * B

                          env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)
                          X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
                          param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                             A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve,
                                             K0=10000)
                          set.seed(111)
                          tmp <- simulate_pheno_ts(Tlim, X0, param.list, env.list, "ceiling")
                          target_parms <- paste(Tlim, Va, Vb, omegaz, omega_Wmax, rho0, rho1,
                                                alpha, delta, sigma_xi, sep='-')
                          target_filename <- paste0('simulate-ceiling', target_parms, '.rds')
                          expect_equal_to_reference(tmp, target_filename, scale=1, tolerance=0.1,
                                                    label=paste('Pars:',
                                                                paste(paste(names(env.list), env.list, sep='='),
                                                                      paste(names(param.list), param.list, sep='='),
                                                                      sep=',', collapse=' ')
                                                                ))
              })
              test_that(paste("right object across parameters constant g, gompertz"), {
                          Vz <- Vz_(Va, Vb, delta, Ve)
                          gamma_sh <-  1 / (omegaz ^ 2 + Vz)
                          rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
                          abar0 <- A;
                          bbar0 <- alpha * B

                          env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)
                          X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
                          param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                             A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve,
                                             K0=10000)
                          set.seed(111)
                          tmp <- simulate_pheno_ts(Tlim, X0, param.list, env.list, "gompertz")
                          target_parms <- paste(Tlim, Va, Vb, omegaz, omega_Wmax, rho0, rho1,
                                                alpha, delta, sigma_xi, sep='-')
                          target_filename <- paste0('simulate-gompertz', target_parms, '.rds')
                          expect_equal_to_reference(tmp, target_filename, scale=1, tolerance=0.1,
                                                    label=paste('Pars:',
                                                                paste(paste(names(env.list), env.list, sep='='),
                                                                      paste(names(param.list), param.list, sep='='),
                                                                      sep=',', collapse=' ')
                                                                ))
                                                    })

              test_that(paste("right object across parameters constant g, thetalogistic"), {
                          Vz <- Vz_(Va, Vb, delta, Ve)
                          gamma_sh <-  1 / (omegaz ^ 2 + Vz)
                          rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
                          abar0 <- A;
                          bbar0 <- alpha * B

                          env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)
                          X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
                          param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                             A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve,
                                             K0=10000, thetalog=1.0)

                          set.seed(111)
                          tmp <- simulate_pheno_ts(Tlim, X0, param.list, env.list, "thetalogistic")
                          target_parms <- paste(Tlim, Va, Vb, omegaz, omega_Wmax, rho0, rho1,
                                                alpha, delta, sigma_xi, sep='-')
                          target_filename <- paste0('simulate-theta-logistic', target_parms, '.rds')
                          expect_equal_to_reference(tmp, target_filename, scale=1, tolerance=0.1,
                                                    label=paste('Pars:',
                                                                paste(paste(names(env.list), env.list, sep='='),
                                                                      paste(names(param.list), param.list, sep='='),
                                                                      sep=',', collapse=' ')
                                                                ))
                                                    })

              test_that(paste("right object across parameters varying g"), {
                          Vz <- Vz_(Va, Vb, delta, Ve)
                          gamma_sh <-  1 / (omegaz ^ 2 + Vz)
                          rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
                          abar0 <- A;
                          bbar0 <- alpha * B

                          param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                             A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve)

                          env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)
                          X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)

                          set.seed(111)
                          tmp <- simulate_pheno_ts(Tlim, X0, param.list, env.list, varying_g=TRUE)
                          target_parms <- paste(Tlim, Va, Vb, omegaz, omega_Wmax, rho0, rho1,
                                                alpha, delta, sigma_xi, sep='-')
                          target_filename <- paste0('simulate-varying-g', target_parms, '.rds')
                          expect_equal_to_reference(tmp, target_filename, scale=1, tolerance=0.1,
                                                    label=paste('Pars:',
                                                                paste(paste(names(env.list), env.list, sep='='),
                                                                      paste(names(param.list), param.list, sep='='),
                                                                      sep=',', collapse=' ')
                                                                ))

              })
            }


context("Poisson mean/variance without plasticity")
delta <- 0
Vb <- 0
alpha <- 0

## remove env variance!
B <- 0
sigma_xi <- 0.0

for( Va in c(0.05, 0.1))
  for( omegaz in sqrt(c(10, 25)))
    for( rho0 in c(0.3, 0.7))
      for( R0 in c(1.1, 1.3, 2)) 
        for( Npop0 in c(5, 50, 300)) {
          ## for some reason Npop0 = 20 cause trouble...
            Vz <- Vz_(Va, Vb, delta, Ve)
            gamma_sh <-  1 / (omegaz ^ 2 + Vz)
            rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
            abar0 <- A;
            bbar0 <- alpha * B

            param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                               A=A, B=B, R0=R0, var_a=Va, Vb=Vb, Ve=Ve, thetalog=5, K0=1000)
            env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)

            Tlim <- 1
            X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
            lambda <- R0 * W_bar(abar0 + bbar0 * delta,
                                 A + B * delta,
                                 omegaz ^ 2, gamma_sh, LOG=FALSE) * Npop0
            Npop_idx <- 5
            reps <- 100

            set.seed(111)
            our_sims <- replicate(reps,
                                  simulate_pheno_ts(Tlim, X0,
                                                    param.list, env.list,
                                                    poisson=TRUE)[2 , Npop_idx])
            our_sims_dd <- replicate(reps,

                                  simulate_pheno_ts(Tlim, X0,
                                                    param.list, env.list, growth_fun="thetalogistic",
                                                    poisson=TRUE)
                                  [2 , Npop_idx])
            set.seed(111)
            rpois_sims <- replicate(reps, rpois(1, lambda))
            test_that(paste(paste("correct Poisson output for lambda~", R0 * Npop0, ":"),
                            paste(
                            paste(names(X0), X0, sep='=', collapse=','),
                            paste(names(param.list), param.list, sep='=', collapse=','),
                            paste(names(env.list), env.list, sep='=', collapse=','),
                            sep='\n', collapse=' ')),
            {
              # for debug --  expect_equal(our_sims, rpois_sims)

              # use se of mean for poisson -- but the actual error is large ...
              se_pois <- sqrt(mean(our_sims) / reps)
              expect_equal(mean(our_sims), mean(rpois_sims), scale=1, tolerance = 4 * se_pois)
              expect_equal(mean(our_sims_dd), mean(rpois_sims), scale=1, tolerance = 4 * se_pois)

              # use variance of sample variance for normal to set expected
              # precision, the se of sample var is about
              se_var_pois <- sqrt(2 / (reps - 1) * var(our_sims) ^ 2)
              # also shift expected variance down for dd (because it's capped?)
              expect_equal(var(our_sims),  var(rpois_sims), scale=1, tolerance = 1.5 * se_var_pois)
              expect_equal(var(our_sims_dd),  0.95 * var(rpois_sims), scale=1, tolerance = 1.5 * se_var_pois)
            })
}
