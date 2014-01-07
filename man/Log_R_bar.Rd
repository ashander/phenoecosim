\name{Log_R_bar}
\alias{Log_R_bar}
\title{Compute LOG population growth rate under stabilizing selection}
\usage{
  Log_R_bar(N, R0, logWbar, K, thetaL)
}
\arguments{
  \item{logWbar}{log average fitness}

  \item{N}{number of individuals in this generation}

  \item{R0}{basic reproductive number}

  \item{K}{carrying capacity}

  \item{thetaL}{theta-logistic parameter for density
  dependence}
}
\description{
  Compute LOG population growth rate under stabilizing
  selection
}
\details{
  Assumes theta-logistic population regulation would be
  good to have separate DD function specified after chevin
  and lande 2010
}

