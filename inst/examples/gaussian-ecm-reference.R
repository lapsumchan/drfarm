# Gaussian optimization baseline on an explicit small, asymmetric starting state.
# Run with the installed development candidate:
# Rscript --vanilla inst/examples/gaussian-ecm-reference.R
local({
  library(drfarm)

  X <- matrix(c(-1, -1, 1, 1), ncol = 1)
  Y <- rbind(c(-0.2, 0.8), c(-2.8, -1.8),
             c(2.2, 3.2), c(0.8, -2.2))
  Theta0 <- matrix(0, 2, 1)
  B0 <- matrix(c(1, 0.5), 2, 1)
  psi0 <- c(1, 2)

  fit <- gaussian.ecm.reference(
    X, Y, Theta0, B0, psi0, lambda1 = 1, lambda2 = 0,
    control = list(max.iter = 500L, objective.tol = 1e-8,
                   score.tol = 1e-6),
    coefficient.control = list(tol = 1e-10)
  )
  print(fit$Theta)
  print(fit$diag.Psi)
  print(fit$diagnostics[c("method", "converged", "termination", "iterations",
                          "objective", "stationarity", "inference.status")])
  print(tail(fit$diagnostics$history, 3))

  # Fresh marginal mean and implied covariance, on exactly the supplied scale.
  predictor.mean <- X %*% t(fit$Theta)
  Sigma <- tcrossprod(fit$B) + diag(fit$diag.Psi)
  residual <- Y - predictor.mean
  objective <- nrow(Y) * as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) / 2 +
    sum(residual * t(solve(Sigma, t(residual)))) / 2 + sum(abs(fit$Theta))
  cat("Independently evaluated observed penalized objective:", objective, "\n")
  stopifnot(
    identical(dim(fit$Theta), c(2L, 1L)),
    all(is.finite(c(fit$Theta, fit$B, fit$diag.Psi))),
    all(fit$diag.Psi > 0),
    abs(objective - fit$diagnostics$objective) < 1e-8,
    identical(fit$diagnostics$inference.status, "not_provided")
  )

  # This reference omits DrFARM's inner debiasing. Inspect termination and the
  # observed score; finite objective agreement is not inferential validation.
})
