\name{Pheno_lande}
\alias{Pheno_lande}
\title{Compute trait change under stabilizing selection as function of environment}
\usage{
  Pheno_lande(t, X, p, G, env.fn, env.args)
}
\arguments{
  \item{t}{timestep}

  \item{X}{parameters (a, b, env)}

  \item{p}{(rho, gamma, A, B)}

  \item{G}{the (constant) G matrix}

  \item{env.fn}{function to compute environment with
  signature function(environment, time, ...)}

  \item{env.args}{extra args for env.fn}
}
\description{
  Compute trait change under stabilizing selection as
  function of environment
}
\details{
  function signature for use with deSolve and 'iterate'
}

