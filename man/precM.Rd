% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{precM}
\alias{precM}
\title{Precision Matrix Estimation Wrapper}
\usage{
precM(X, method = "glasso", standardize = TRUE)
}
\arguments{
\item{X}{Data matrix \eqn{n \times p}.}

\item{method}{One of \code{"glasso"}, \code{"NL"}, or \code{"QO"}.}

\item{standardize}{Whether to standardize \code{X}.}
}
\value{
A \eqn{p \times p} precision matrix estimate.
}
\description{
Precision Matrix Estimation Wrapper
}
