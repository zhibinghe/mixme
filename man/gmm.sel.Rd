% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{gmm.sel}
\alias{gmm.sel}
\title{Component selection in Gaussian mixture model}
\usage{
gmm.sel(
  X,
  M,
  alpha,
  lambda,
  muK,
  sigK,
  maxit = 5000,
  tol = 1e-04,
  verb = FALSE
)
}
\arguments{
\item{X}{Surrogate data or observations in main study}

\item{M}{initial number of components}

\item{alpha}{tuning parameter controlling degree of penalty}

\item{lambda}{A vector of size K containing initial value of mixing probabilities}

\item{muK}{A matrix of size K * p containing initial values of component mean, K-means center is specified if NULL}

\item{sigK}{An array with size p*p*K containing K p*p matrix, each of which is initial values of component variance,
K-means variance is specified if NULL}

\item{maxit}{The maximum number of iterations}

\item{tol}{The convergence criterion}

\item{verb}{If true, then various updates are printed during each iteration of the algorithm}
}
\value{
A list with the following elements:
\itemize{
\item{pi} The final mixing probabilities, only containing selected nonzero lambda
\item{mu} The final mean vectors
\item{sigK} The final variance matrix of each component
\item{posterior} Membership probability
\item{iter} Number of iterations
\item{loglik} log-likelihood of observed data
\item{bic_loglik} bic criteria 
}
}
\description{
Component selection in Gaussian mixture model
}
