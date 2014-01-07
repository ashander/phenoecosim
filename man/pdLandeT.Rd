\name{pdLandeT}
\alias{pdLandeT}
\title{Compute trait + demographic change over T generations under stabilizing selection as
function of environment (after Lande Chevin)}
\usage{
  pdLandeT(T, X, params, G, env_args)
}
\arguments{
  \item{T}{end time, assuming start time of 1}

  \item{X}{parameters (a, b, env, logN)}

  \item{params}{a list with (gamma, Oz2, AB, R0, K, theta)}

  \item{G}{the (constant) G matrix}

  \item{env_args}{extra args for env.fn}
}
\description{
  Compute trait + demographic change over T generations
  under stabilizing selection as function of environment
  (after Lande Chevin)
}
\details{
  function signature as for use with deSolve imposes
  fitness based on Lande and demography after CL 2010
}

