library(phenoecosim)

ERROR_TOL <- 0.15
for(rho in 1:8 * 0.1) {
  context(paste("Environmental simulation testing for rho =", rho))
  max_T = 1000
  start_T = 1
  for(delta in 1:8)
    for(fractgen in 5:10)
      for(sigma_xi in c(0.5, 0.75, 1, 3)) {
        env.list <- list(delta=delta, sigma_xi=sigma_xi, rho_tau=rho, fractgen=fractgen)
        env_dev_sel = make_env(T=max_T, env_args=env.list)
        par_names <- paste0(paste(names(env.list), env.list, sep='='), collapse=',')
        test_that(paste("Environmental simulation produces expected correlations, mean, sd for", par_names), {
          envcor <- cor(env_dev_sel[1, start_T:max_T], env_dev_sel[2,  start_T:max_T])
          envmeans <- rowMeans(env_dev_sel)
          envsd <- apply(env_dev_sel, 1, sd)
          expect_equal(envcor, rho, tolerance=ERROR_TOL, scale=sigma_xi)
          expect_equal(envmeans, rep(delta, 2), tolerance=ERROR_TOL, scale=sigma_xi)
          expect_equal(envsd, rep(sigma_xi, 2), tolerance=ERROR_TOL, scale=sigma_xi)
        })
        test_that(paste("Environmental simulation produces mean", delta, "for", par_names), {
        })
      }
}
