library(phenoecosim)

dd_forms <- c("ceiling",  "thetalogistic", "gompertz")
relative_error_from_Npop_eq <- c(ceiling=0.0025,
                                 thetalogistic=0.05, ## 
                                 gompertz=0.05 ## 
                                 )
debug_eq <- FALSE
debug_similar <- FALSE

## unvarying parameters
A <- 0
B <- 2 #environmental cline in the optimum phenotype
fractgen <- 5
tau <- 1 / fractgen
Ve <- 0.5
Tchange <- 0 #time for a change in env regime
Npop0 <- 1e4
rho1 <- 0.0
sigma_xi <- 0.1
omega_z <- sqrt(20)
omega_Wmax <- 1.7 #maximum fitness across environments

for (form in dd_forms) {
  context(paste(form, "near K for long time"))
  delta = 0
  alpha = 0.0
  ## varying params
  for( Tlim in c(1000, 3000))
    for( Va in c(0.05, 0.1))
      for( Vb in c(0.025, 0.05))
        for( omegaz in sqrt(c(10, 25)))
          for( K0 in c(1000, 10000, 20000))
            for( rho0 in c(0.4, 0.7)){
              Vz <- Vz_(Va, Vb, delta, Ve)
              gamma_sh <-  1 / (omegaz ^ 2 + Vz)
              rmax <- log(omega_Wmax) -  1 / 2 *  log(1 + Vz / (omegaz ^ 2))
              abar0 <- A;
              bbar0 <- alpha * B

              env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho1, fractgen=fractgen)

              ## r discounted by standing load in ref \environment
              rmaxnosh <- log(omega_Wmax) -  1 / 2 *  log(1 + (Va + Ve ) / (omegaz ^ 2))

              if (form == 'gompertz') 
                Npop_eq <- K0 ^ (rmaxnosh / log(omega_Wmax))
              if (form == 'ceiling')
                Npop_eq <- K0
              if (form == 'thetalogistic')
                Npop_eq <- K0

              ## Initial conditions
              Npop0 <- Npop_eq

              X0 <- c(zbar=NA, abar=abar0, bbar=bbar0, Wbar=NA, Npop=Npop0, theta=NA)
              param.list <- list(Vz=Vz, gamma_sh=gamma_sh, rmax=rmax, omegaz=omegaz,
                                 A=A, B=B, R0=omega_Wmax, var_a=Va, Vb=Vb, Ve=Ve,
                                 K0=K0, thetalog=2)
              set.seed(111)
              traj_dd <- simulate_pheno_ts(Tlim, X0, param.list, env.list, form)

              test_that(paste("End of traj at carrying cap for K=", K0, ":", form),
                        {
                          Npop_idx = 5
                          tail_len = 200
                          tail_idxs = (Tlim - tail_len + 1):Tlim
                          expect_equal(traj_dd[tail_idxs, Npop_idx], rep(Npop_eq, tail_len),
                                       tolerance = relative_error_from_Npop_eq[form], scale = Npop_eq)

                          if (debug_eq) {
                            browser()
                            par(mfrow=c(1, 2))
                            plot(traj_dd[tail_idxs , Npop_idx], ylog=TRUE)
                            form = "thetalogistic"
                            traj_dd <- simulate_pheno_ts(Tlim, X0, param.list, env.list, form)
                            plot(traj_dd[tail_idxs , Npop_idx], ylog=TRUE)
                          }
                        })
              if (form == 'gompertz')
                test_that(paste("Trajectories similar for small thetalog and gompertz. K=", K0),
                          { 
                            Npop_idx = 5
                            set.seed(111)
                            traj_gom <- simulate_pheno_ts(Tlim, X0, param.list, env.list, form)
                            param.list['thetalog'] <- 0.1
                            X0['Npop'] <- K0
                            set.seed(111)
                            traj_tl <- simulate_pheno_ts(Tlim, X0, param.list, env.list,
                                                         'thetalogistic')
                            expect_equal(traj_gom[ , Npop_idx], traj_tl[ , Npop_idx],
                                         tolerance = .2, scale=Npop_eq)
                            if (debug_similar) {
                              browser()
                              par(mfrow=c(1, 2))
                              plot(traj_gom[ , 5], ylog=TRUE)
                              plot(traj_tl[ , 5], ylog=TRUE)
                            }
                          })
              if (form == 'ceiling')
                test_that(paste("Trajectories similar for large thetalog and ceiling K=", K0),
                          { 
                            Npop_idx = 5
                            set.seed(111)
                            traj_ceil <- simulate_pheno_ts(Tlim, X0, param.list, env.list, form)
                            param.list['thetalog'] <- 5
                            set.seed(111)
                            traj_tl <- simulate_pheno_ts(Tlim, X0, param.list, env.list,
                                                         'thetalogistic')
                            expect_equal(traj_ceil[ , Npop_idx], traj_tl[ , Npop_idx],
                                         tolerance = 0.2, scale = Npop_eq)
                          })

            }
}
