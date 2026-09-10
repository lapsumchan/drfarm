# Exact unequal-variance coefficient fixture for the opt-in weighted solver.
# Run against the installed candidate:
# Rscript --vanilla inst/examples/weighted-update.R
local({
  library(drfarm)

  X <- matrix(c(1, -1) / sqrt(2), ncol = 1)
  Y <- X %*% matrix(c(8, 5.5), nrow = 1)
  sigma <- c(1, 0.25)  # Response variances, not standard deviations.
  lambda1 <- 2
  lambda2 <- 5

  weighted <- remMap.weighted(
    X, Y, lambda1 = lambda1, lambda2 = lambda2, sigma = sigma,
    control = list(tol = 1e-8, max.sweeps = 1000L,
                   root.tol = 1e-14, root.maxit = 200L)
  )
  print(weighted$Theta0)
  print(weighted$diagnostics)

  # Independent reference: h=(1,4), soft(s,2)=(6,20), t=1 gives b=(3,4).
  b <- drop(weighted$Theta0)
  coefficient.error <- max(abs(b - c(3, 4)))
  residual <- Y - X %*% t(weighted$Theta0)
  objective <- sum(colSums(residual^2) / (2 * sigma)) +
    lambda1 * sum(abs(b)) + lambda2 * sqrt(sum(b^2))
  gradient <- -drop(crossprod(X, residual)) / sigma +
    lambda1 * sign(b) + lambda2 * b / sqrt(sum(b^2))
  cat("Maximum coefficient error:", coefficient.error, "\n")
  cat("Objective (analytic target 56):", objective, "\n")
  cat("Independent KKT gradient infinity norm:", max(abs(gradient)), "\n")
  stopifnot(
    identical(dim(weighted$Theta0), c(2L, 1L)),
    coefficient.error < 1e-7,
    abs(objective - 56) < 1e-7,
    max(abs(gradient)) < 1e-7
  )

  # This is a coefficient-subproblem check. It does not certify the outer
  # DrFARM factor/debiasing iteration or its inference procedures.
})
