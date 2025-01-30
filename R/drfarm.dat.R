#' Example data for DrFARM
#'
#' A list containing example design and response matrices, along with a true coefficient matrix,
#' intended for demonstration of the DrFARM model
#'
#' @format A list with the following components:
#' \describe{
#'   \item{X}{A \eqn{500 \times 10} numeric matrix of predictors}
#'   \item{Y}{A \eqn{500 \times 5} numeric matrix of responses (outcomes)}
#'   \item{Theta.t}{A \eqn{10 \times 5} numeric matrix representing the true coefficients used in data generation}
#' }
#'
#' @details
#' The data were simulated under a DrFARM model. \code{X} and \code{Y} were generated using a true
#' coefficient matrix \code{Theta.t}, plus a rank-2 factor structure and random noise. Specifically:
#' \enumerate{
#'   \item \eqn{X} is generated from a standard normal distribution, \eqn{\mathrm{N}(0, I_p)}
#'   \item A subset of entries in \eqn{\Theta.t} are nonzero, introduced randomly, to represent true signals.
#'   \item Residual factors (\eqn{Z}) and factor loadings (\eqn{B}) are introduced so that
#'     \eqn{Y = X \Theta^T + Z B^T + E}
#'   \item \eqn{E} is an independent noise matrix
#' }
#'
#' @examples
#' \dontrun{
#' data(drfarm.dat)
#' str(drfarm.dat)
#'
#' # Access design matrix and outcome matrix
#' X <- drfarm.dat$X
#' Y <- drfarm.dat$Y
#'
#' # True coefficient matrix
#' Theta.true <- drfarm.dat$Theta.t
#' }
#'
#' @source Simulated data
"drfarm.dat"
