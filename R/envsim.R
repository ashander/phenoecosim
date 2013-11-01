#' Compute a white noise environment with a shift
#' @param env previous year's environment NOT used in current formulation
#' @param t the time point
#' @param rho the correlation between environment of selection and the cue
#' @param env.args other args:
#' t.jump the location for the jump,
#' delta env change that takes place at t.jump
#' sd noise in env
#' sdc noise in cue
#' @details nada
#' @import MASS
#' @export
Env_shift <- function(env, t, rho, env.args){
    ## signature function(environment, time, ...)
    t.jump <- env.args[['t.jump']]
    k <- env.args[['k']]
    sd <- env.args[['sd']]
    sdc <- env.args[['sdc']] 
    ## for red noise do somehting like k * env  + sqrt(1 - k^2) * rnorm (n=1, mean=0, sd=sd)) but hard
    ## to figure this out with the cue / correlation thing
    out <- mvrnorm(1, mu = c(0,0) , Sigma = matrix(c(sd, sd*sdc*rho, sd*sdc*rho, sdc), nrow=2))
    if( t < t.jump)
      return(out)
    else if(t >= t.jump)
        return(env.args[['delta']] + out )
    else
        return("You should never see this")
}

#' Simulate one step of a linearly changing environment, possibly with noise
#' @param env environment in previous time step
#' @param t the time point
#' @param rho the correlation between environment of selection and the cue
#' @param env.args other args:
#' eta rate of environmental change,
#' sd noise in env
#' sdc noise in cue
#' @details other args in env args
#' @import MASS
#' @export
Env_smooth <- function(env, t, rho, env.args){
    ## signature function(environment, time, ...)
    eta <- env.args[['eta']]
    sdc <- env.args[['sdc']]
    sd <- env.args[['sd']]
    ## cue <- env.args[['cue']]     ## want to use  as current cue, but ... argument passing harder, maybe
    eta + mvrnorm(1, mu = c(env,env) , Sigma = matrix(c(sd, sd*sdc*rho, sd*sdc*rho, sdc), nrow=2))
}



#' Env with autocorrelation: need to add this and makesure I can re-obtain results from democonstriants_ms
