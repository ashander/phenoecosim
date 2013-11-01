#' Compute population growth rate under stabilizing selection
#' @param N number of individuals in this generation
#' @param R0 basic reproductive number
#' @param Wbar average fitness
#' @param K carrying capacity
#' @param thetaL theta-logistic parameter for density dependence
#' @details Assumes theta-logistic population regulation
#' would be good to have separate DD function specified after
#' chevin and lande 2010
R_bar <- function(N, R0, Wbar, K, thetaL)
    R0^(1 - (N/K)^thetaL) * Wbar

#' Compute log mean fitness under stabilizing selection
#' @param zbar mean trait prior to selection
#' @param theta environmental optimum
#' @param Oz2 strength of stabilizing selection
#' @param gamma 1/(Oz2 + Vz2), with Vz2 phenotypic variance 
#' @param log (=TRUE) whether to return the log of fitness  
#' @details no details
W_bar <- function(zbar, theta, Oz2, gamma, log=TRUE){
    if(log)
        return(0.5 * log(gamma * Oz2) -gamma /2 * (zbar - theta)^2)## compare with eqn 2c lande 2009
    else
        return( sqrt(gamma* Oz2) * exp(- gamma /2 * (zbar - theta)^2)) 
    ## old code
    ## Log.W.bar <- function(z.bar, theta, e=0)  
    ##     0.5 * log(gamma * omega2) - gamma / 2 * (z.bar - theta)^2 # eqn 2c
}


#' Mean fitness -- with constraint, need to make as function of environment...
#' @param x  the environment
#' @param LL lethal limit (symmetric)
#' @details fitness surf infinitesimal optimum  need to revisit
stab_curve<- function(x, LL) -(x-LL)*(x+LL) * 1/LL^2

#' Mean fitness -- with constraint, need to make as function of environment...
#' @param env the environment
#' @param LL lethal limit 
#' @details need to revisit this
Log_W_bar_lethal <- function(env, LL=4)   {
  if(env > -LL & env < LL)
    return(stab_curve(env, LL))
  else
    return(0)
}


#' G matrix
#' @param Gaa additive genetic variance in RN intercept
#' @param Gbb additive genetic variance in RN slope
#' @param Gab additive genetic covariance between RN slope and height
#' @export
G <- function(Gaa, Gbb, Gab)
    matrix( c(Gaa, Gab, Gab, Gbb), nrow=2) #eqn 3a

#' Selection on plasticity as as function of environment assuming stabilizing selection
#' @param gamma = 1/(Oz2 + Vz2), with Vz2 phenotypic variance
#' @param A the environmental optimum RN int
#' @param B the environmental optimum RN slope
#' @param a current value of RN int
#' @param b current value of RN slope
#' @param e.t environment now
#' @param e.plast env that cues plasticity
#' @export
Beta <- function(gamma, A, B, a, b, e.t, e.plast){
  ## eqn 3b
  beta1 <- a - A + b * e.plast - B * e.t
  -gamma * c( beta1 , 
             beta1 * e.plast)
}

#' Va additive genetic variance in the phenotype as a function of the environment
#' @param e.plast environment in which plastiicty is cued
#' @param GG additive genetic variance-covariance matrix
#' @export
Va <-function(e.plast, GG){
    ##    first three terms of eqn 1c
    ##if(nrow(GG) != 2)
    ##    stop("g is not a matrix of correct size")
    ee <- c(1, e.plast)
    ee %*% GG %*% ee
  }

#' Compute trait change under stabilizing selection as function of environment
#' @param t timestep
#' @param X parameters (a, b, env)
#' @param p (rho, gamma, A, B)
#' @param G the (constant) G matrix
#' @param env.fn function to compute environment with signature function(environment, time, ...)
#' @param env.args extra args for env.fn
#' @details function signature for use with deSolve and 'iterate'
#' @import deSolve
#' @export
Pheno_lande <- function(t, X, p, G, env.fn, env.args) {
  ## state vars
  a <- X[1] # elevation
  b <- X[2]  # slope
  env <- X[3]
  ## params
  rho <- p[1]
  gamma <- p[2]
  Oz2 <- p[3]
  A <- p[4]
  B <- p[5]
  EE <- env.fn(env, t, rho, env.args)
  env <- EE[1]
  e.plast <- EE[2]
  z.bar <- a + b * e.plast
  log.W.bar <- W_bar(z.bar, A + B*env, Oz2, gamma)
  beta <- Beta(gamma, A, B, a, b, env, e.plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e.plast, G) 
  list(c(X, env, z.bar, va))
}

#' Compute trait + demographic change under stabilizing selection as function of environment
#' @param t timestep
#' @param X parameters (a, b, env, N)
#' @param p (rho, gamma, A, B, R0)
#' @param G the (constant) G matrix
#' @param env.fn function to compute environment with signature function(environment, time, ...)
#' @param env.args extra args for env.fn
#' @details function signature for use with deSolve and 'iterate', imposes fitness based on
#' assumed tolerance curve
#' @import deSolve
#' @export
Pheno_demo <- function(t, X, p, G, env.fn, env.args) {
    ## TODO fix error in demo, returns NA after Ni 10000 at time zero
    ## first issue is R0 missing from call... 
  ## state vars
  a <- X[1] # elevation
  b <- X[2]  # slope
  env <- X[3]
  N = X[4]

  ## params
  rho <- p[1]
  gamma <- p[2]
  Oz2 <- p[3]

  A <- p[4]
  B <- p[5]
  R0 <- p[6]
  
  EE <- env.fn(env, t, rho, env.args)
  env <- EE[1]
  e.plast <- EE[2]

  log.W.bar <- Log_W_bar_lethal(env)
  N <-   log.W.bar * R0 * N ## LOG demo update based on mean fitness
  
  beta <- Beta(gamma, A, B, a, b, env, e.plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e.plast, G) 
  list(c(X, env, N, va))
}


#' Compute trait + demographic change under stabilizing selection as function of environment via Lande Chevin
#' @param t timestep
#' @param X parameters (a, b, env, N)
#' @param p (rho, gamma, A, B, R0, K, thetaL)
#' @param G the (constant) G matrix
#' @param env.fn function to compute environment with signature function(environment, time, ...)
#' @param env.args extra args for env.fn
#' @details function signature for use with deSolve and 'iterate', imposes fitness based on Lande
#' and demography after CL 2010
#' @import deSolve
#' @export
Pheno_demo_lande <- function(t, X, p, G,  env.fn, env.args) {
  ## state vars
  a <- X[1] # elevation
  b <- X[2]  # slope
  env <- X[3]
  N = X[4]

  ## params
  rho <- p[1]
  gamma <- p[2]
  Oz2 <- p[3]

  A <- p[4]
  B <- p[5]
  R0 <- p[6]
  K <- p[7]
  thetaL <- p[8]
  
  EE <- env.fn(env, t, rho, env.args)
  env <- EE[1]
  e.plast <- EE[2]
  
  z.bar <- a + b * e.plast
  log.W.bar <- W_bar(z.bar, A + B*env, Oz2, gamma)
  N <- R_bar(N, R0, exp(log.W.bar), K, thetaL) * N ## demo update based on mean fitness
  
  beta <- Beta(gamma, A, B, a, b, env, e.plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e.plast, G) 
  list(c(X, env, N, va))
}

#' Compute trait + demographic change under stabilizing selection as function of environment via Lande Chevin
#' but with intervention
#' @param t timestep
#' @param X parameters (a, b, env, N)
#' @param p (rho, gamma, A, B, R0, K, thetaL)
#' @param G the (constant) G matrix
#' @param env.fn function to compute environment with signature function(environment, time, ...)
#' @param env.args extra args for env.fn
#' @param interv a function of the time, policy to modify the cue
#' @details function signature for use with deSolve and 'iterate', imposes fitness based on Lande
#' and demography after CL 2010
#' @import deSolve
#' @export
Pheno_demo_econ <- function(t, X, p, G,  env.fn, env.args, interv) {
  ## state vars
  a <- X[1] # elevation
  b <- X[2]  # slope
  env <- X[3]
  N = X[4]

  ## params
  rho <- p[1]
  gamma <- p[2]
  Oz2 <- p[3]

  A <- p[4]
  B <- p[5]
  R0 <- p[6]
  K <- p[7]
  thetaL <- p[8]
  
  EE <- env.fn(env, t, rho, env.args)
  env <- EE[1]
  e.plast <- EE[2] + interv(t)
  
  z.bar <- a + b * e.plast
  log.W.bar <- W_bar(z.bar, A + B*env, Oz2, gamma)
  N <- R_bar(N, R0, exp(log.W.bar), K, thetaL) * N ## demo update based on mean fitness
  
  beta <- Beta(gamma, A, B, a, b, env, e.plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e.plast, G) 
  list(c(X, env, N, va))
}



## for examples use #' examples \dontrun{}


## @knitr tolerance-functions
## fitsurf Lande (mean) -- not sure what i'm doing with this
## lande.fit <- function(x) exp( - (x * (E.B-B))^2/(2 * 50))



