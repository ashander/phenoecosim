\name{pdLande}
\alias{pdLande}
\title{Compute trait + demographic change over one generation under stabilizing selection as
function of environment (after Lande Chevin)}
\usage{
  pdLande(t, X, params, G, env_args)
}
\arguments{
  \item{t}{timestep}

  \item{X}{parameters (a, b, env, logN)}

  \item{params}{a list with (Oz2, AB, R0, K, theta)}

  \item{G}{the (constant) G matrix}

  \item{env_args}{extra args for env.fn}
}
\description{
  Compute trait + demographic change over one generation
  under stabilizing selection as function of environment
  (after Lande Chevin)
}
\details{
  function signature as for use with deSolve imposes
  fitness based on Lande and demography after CL 2010
}

