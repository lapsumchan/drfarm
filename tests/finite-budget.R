# A finite outer budget is a failure/status contract, not a convergence proof.
library(drfarm)
set.seed(20260909)
n <- 80L
p <- 3L
q <- 6L
X <- scale(matrix(rnorm(n * p), n, p))
Z <- matrix(rnorm(n), n, 1)
B <- matrix(c(0.6, -0.8, 0.5, 1, -0.7, 0.9), q, 1)
Theta <- matrix(c(0.3, 0, 0, 0.2, 0, 0, 0, -0.25, 0, 0, 0.15, 0,
                  0, 0, 0.2, 0, 0, -0.1), q, p)
Y <- X %*% t(Theta) + Z %*% t(B) + matrix(rnorm(n * q), n, q)
Theta0 <- remMap.one(X, Y, standardize = FALSE, lambda1 = 8, lambda2 = 8)
M <- solve(crossprod(X) / n)
fit <- DrFARM.one(X, Y, Theta0, M, k = 1,
                   lambda1 = 8, lambda2 = 8,
                   standardize = FALSE, max.iter = 1)
stopifnot(identical(dim(fit$Theta), c(q, p)),
          identical(dim(fit$B), c(q, 1L)),
          identical(dim(fit$E.Z), c(n, 1L)),
          all(is.finite(fit$Theta)), all(is.finite(fit$B)),
          all(is.finite(fit$E.Z)), all(is.finite(fit$diag.Psi)),
          all(fit$diag.Psi > 0),
          is.list(fit$diagnostics),
          isFALSE(fit$diagnostics$converged),
          is.numeric(fit$diagnostics$iterations),
          fit$diagnostics$iterations <= 1,
          is.character(fit$diagnostics$termination),
          length(fit$diagnostics$termination) == 1,
          nzchar(fit$diagnostics$termination),
          is.finite(fit$diagnostics$loss),
          fit$diagnostics$max.iter == 1)
cat("PASS DF-06 finite outer budget retains a finite candidate and reports",
    fit$diagnostics$termination, "after", fit$diagnostics$iterations,
    "iteration(s); no convergence claim.\n")
