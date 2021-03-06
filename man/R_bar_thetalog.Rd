% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{R_bar_thetalog}
\alias{R_bar_thetalog}
\title{Compute population growth rate under stabilizing selection and theta-logistic regulation}
\usage{
R_bar_thetalog(R0, Wbar, N, K0, thetaL)
}
\arguments{
\item{R0}{basic reproductive number}

\item{Wbar}{average fitness}

\item{N}{number of individuals in this generation}

\item{K0}{carrying capacity with optimal trait}

\item{thetaL}{theta-logistic parameter for density dependence}
}
\description{
Compute population growth rate under stabilizing selection and theta-logistic regulation
}
\details{
Assumes theta-logistic population regulation
would be good to have separate DD function specified after
chevin and lande 2010
}

