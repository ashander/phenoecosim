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
  test_that(paste("right object across parameters for delta=",delta), {
              ## varying params
              for( Tlim in c(250, 500, 1000))
                for( Va in c(0.05, 0.1))
                  for( Vb in c(0.025, 0.05))
                    for( omegaz in sqrt(c(10, 25)))
                      for( rho0 in c(0.3, 0.7))
                        for( alpha in c(0.3, 0.7))
                        {
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

                          ## ceiling
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
                          ## gompertz
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

                          ## thetalogistic
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
                        }
})
