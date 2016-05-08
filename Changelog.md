0.2.0 / 2016-05-07
==================

  * Poisson sampling for density-independent
    - returns growth is Poisson(f(N(t)) where it was previously
      N(t+1) = f(N(t))
    - pass parameter `poisson = TRUE` (a logical)
    - tests of mean/variance agreement with `rpois` over one generation

  * Enable density-dependence
    - Gompertz, theta-logistic, ceiling growth
    - specify through argument `density` with values "ceiling", "gompertz",
      "thetalogistic", or "independent" (density-independent)
    - all require argument `K0`, which defaults to 10,000, but theta-logistic
      also requires parameter `thetalog`

  * Add simulation runners, with functions:
    - run_sims_on_grid: runs things, and is given a summarizer
    - summary_id: an 'identity' summarizer; does nothing
    - simulate_timeseries: a wrapper for simulate_pheno_ts that computes a few
      things and passes them simulate_pheno_ts. notably, it decides whether to use
      stationary initial conditions. if these are needed, they are generated
      separately, in R, and used to initialize the (c++) simulate_pheno_ts
