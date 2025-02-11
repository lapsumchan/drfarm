% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{DrFARM.whole}
\alias{DrFARM.whole}
\title{Fit all \code{DrFARM} models at once}
\usage{
DrFARM.whole(
  X,
  Y,
  Theta0,
  precM,
  k,
  lambda1.opt,
  lambda2.opt,
  K = NULL,
  C = NULL,
  standardize = TRUE,
  gamma = 1,
  thres = 1e-04,
  rotate = "none",
  scores = "regression",
  fm = "ml"
)
}
\arguments{
\item{X}{A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed}

\item{Y}{A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed}

\item{Theta0}{A size \eqn{q \times p} initial coefficient matrix}

\item{precM}{A size \eqn{p \times p} precision matrix}

\item{k}{A positive integer specifying the number of latent factors used in DrFARM}

\item{lambda1.opt}{The chosen lasso tuning parameter from \code{remMap.whole}}

\item{lambda2.opt}{The chosen group-lasso tuning parameter from \code{remMap.whole}}

\item{K}{An optional size \eqn{n \times n} kinship matrix. Default is \code{NULL}}

\item{C}{A size \eqn{q \times p} integer matrix specifying model inclusion and penalty rules:
\itemize{
  \item 0 = excluded from model
  \item 1 = penalized
  \item 2 = unpenalized
}
If \code{null}, all entries are penalized (1) by default}

\item{standardize}{A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}}

\item{gamma}{A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{1}}

\item{thres}{A numeric value for the convergence threshold. Default is \eqn{1 \times 10^{-4}}}

\item{rotate}{Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"none"})}

\item{scores}{Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"regression"})}

\item{fm}{Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"ml"})}
}
\value{
A list containing:
\item{Theta}{The EBIC-chosen \eqn{q \times p} coefficient matrix}
\item{B}{A size \eqn{q \times k} factor loading matrix}
\item{E.Z}{A size \eqn{n \times k} latent factor score matrix}
\item{diag.Psi}{A length-\eqn{q} vector of uniquenesses (diagonal of \eqn{\Psi}).}
\item{lambda1.opt}{Optimal \eqn{\lambda_1} value}
\item{lambda2.opt}{Optimal \eqn{\lambda_2} value}
}
\description{
Searches over a 2D grid of (\code{lambda1}, \code{lambda2}) pairs, computes \code{EBIC} for each, and
returns the \eqn{q x p} coefficient matrix that minimizes the EBIC, along with corresponding
factor loadings, factor scores, and uniquenesses
}
