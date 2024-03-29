% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{gendat}
\alias{gendat}
\title{Simulate Data for Gaussian Mixture Model Combined with Measurement Error and Linear Model}
\usage{
gendat(
  n1,
  n2,
  lambda,
  muK,
  sigK,
  alpha,
  A,
  sigE,
  beta,
  sigma2,
  is.diff = FALSE
)
}
\arguments{
\item{n1}{Number of units in the main study}

\item{n2}{Number of units in the validation study}

\item{lambda}{A vector of size K consisting of mixing probabilities}

\item{muK}{A matrix of size K*p containing Gaussian component mean}

\item{sigK}{An array of size p*p*K, each element is the Gaussian component variance}

\item{alpha}{A vector of size p consisting of intercept in measurement error model}

\item{A}{A matrix of size p*p consisting of slope in measurement error model}

\item{sigE}{A matrix of size p*P consisting of variance of noise in measurement error model}

\item{beta}{A vector of size K consisting of linear coefficients in linear model}

\item{sigma2}{A scalar of variance of random noise in linear model}

\item{is.diff}{If true, the distributions of data in main study and validation study are different}
}
\value{
A list with the following elements
\itemize{
\item{Zm} Surrogate data in the main study
\item{Xm} True observations in the main study
\item{Zv} Surrogate data in the validation study
\item{Xv} True observations in the validation study
\item{y} The outcome of linear model
\item{cluster.m} The true cluster labels of units in the main study
}
}
\description{
Simulate Data for Gaussian Mixture Model Combined with Measurement Error and Linear Model
}
