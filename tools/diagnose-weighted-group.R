# Diagnostic of the historical unequal-variance group update.
# Kept outside package tests because the maintenance patch preserves the update.
# Usage: Rscript --vanilla weighted_group_counterexample.R /path/to/package-library
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the library containing the chosen drfarm installation.")
library(drfarm, lib.loc = args[1])
cat("Package:", as.character(packageVersion("drfarm", lib.loc = args[1])), "\n")
cat("Library:", normalizePath(args[1]), "\n")
X <- matrix(c(1, -1) / sqrt(2), 2, 1)
a <- c(3, 4)
Y <- X %*% matrix(a, 1, 2)
sigma <- c(1, 2)
lambda_group <- 1
fit <- drfarm:::remMap(X, Y, lamL1 = 0, lamL2 = lambda_group, sigma = sigma)
b <- drop(fit$phi)
gradient <- colSums(X^2) * (b - a) / sigma +
  lambda_group * b / sqrt(sum(b^2))
expected <- c(2.4, 2.4)
expected_gradient <- c(-0.6, -0.8) + 1 / sqrt(2)
cat("Historical coefficients:", format(b, digits = 16), "\n")
cat("Expected historical coefficients:", format(expected, digits = 16), "\n")
cat("Coefficient max absolute discrepancy:", max(abs(b - expected)), "\n")
cat("Natural weighted-objective gradient:", format(gradient, digits = 16), "\n")
cat("Gradient max absolute discrepancy:", max(abs(gradient - expected_gradient)), "\n")
cat("Gradient infinity norm:", max(abs(gradient)), "\n")
cat("Stationarity for RSS/(2 sigma) + ||b||_2:",
    if (max(abs(gradient)) < 1e-10) "SATISFIED" else "NOT SATISFIED", "\n")
if (!is.null(fit$diagnostics)) {
  cat("Native coefficient-change diagnostics (distinct from stationarity):\n")
  print(fit$diagnostics)
}
stopifnot(max(abs(b - expected)) < 1e-10,
          max(abs(gradient - expected_gradient)) < 1e-10,
          max(abs(gradient)) > 0.1)
cat("Counterexample reproduced; historical operation retained. No general-method claim.\n")
