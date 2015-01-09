#' Mean fitness -- with constraint, need to make as function of environment...
#' @param x  the environment
#' @param LL lethal limit (symmetric)
#' @details fitness surf infinitesimal optimum  need to revisit
#' @export
stab_curve<- function(x, LL) -(x-LL)*(x+LL) * 1/LL^2

#' Mean fitness -- with constraint, need to make as function of environment...
#' @param env the environment
#' @param LL lethal limit 
#' @details need to revisit this
#' @export
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
  N <- R_bar( R0, exp(log.W.bar)) # for _dd N , K, thetaL) * N ## demo update based on mean fitness
  
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
  N <- R_bar(R0, exp(log.W.bar)) # for _dd N, K, thetaL) * N ## demo update based on mean fitness
  
  beta <- Beta(gamma, A, B, a, b, env, e.plast)
  X <- c(a, b) + G %*% beta
  va <- Va(e.plast, G) 
  list(c(X, env, N, va))
}



## for examples use #' examples \dontrun{}


##  tolerance-functions
## fitsurf Lande (mean) -- not sure what i'm doing with this
## lande.fit <- function(x) exp( - (x * (E.B-B))^2/(2 * 50))




