# Deterministic mathematical fixtures. Run with R CMD check or:
# Rscript --vanilla tests/numerical-contracts.R (after installing drfarm).
# These tests establish the specified finite cases, not inferential validity.
library(drfarm)

assert_close <- function(actual, expected, label, tolerance = 1e-10) {
  if (!identical(dim(actual), dim(expected)) ||
      length(actual) != length(expected) ||
      any(!is.finite(actual)) || any(!is.finite(expected))) {
    stop(label, ": shape or finite-value contract failed")
  }
  discrepancy <- max(abs(actual - expected))
  if (discrepancy > tolerance) {
    stop(label, ": max absolute error ", format(discrepancy),
         " exceeds ", format(tolerance))
  }
  cat("PASS", label, "; max absolute error =", format(discrepancy), "\n")
}

assert_error <- function(expr, label) {
  result <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(result, "error")) stop(label, ": expected an error")
  cat("PASS", label, "; message =", conditionMessage(result), "\n")
}

# DF-01/02, P05/P06, N01/N02. Orthogonal, centered, unequal-norm X;
# n=8, p=2, q=3 makes transposed output impossible to hide. The noise
# is orthogonal to X, so the exact OLS coefficient is A (independent QR
# provides a second reference). The package objective has RSS/2, not RSS/(2n).
u <- rep(c(-1, 1), 4)
v <- rep(c(-1, -1, 1, 1), 2)
w <- c(rep(-1, 4), rep(1, 4))
z <- u * v
X <- cbind(u, 2 * v)
A <- rbind(c(1.5, -0.7, 2), c(-0.5, 1.2, 0.35))
noise <- outer(w, c(0.5, -0.25, 0.75))
Y <- X %*% A + noise
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
ols <- remMap.one(X, Y, standardize = FALSE, lambda1 = 0, lambda2 = 0)
assert_close(ols, t(A), "DF-01/02 orthogonal nonsquare OLS truth")
assert_close(ols, t(qr.coef(qr(X), Y)), "DF-02 independent QR reference")
assert_close(X %*% t(ols), Y - noise, "DF-01 prediction orientation")

# Independent solution of each orthogonal predictor block:
# min_b d_j ||b-a_j||^2/2 + lambda1 ||b_pen||_1 + lambda2 ||b_pen||_2,
# with C=0 constrained to zero and C=2 absent from both penalties.
# Only equal response variance gives this Euclidean sparse-group prox.
C <- matrix(c(1L, 2L, 0L, 1L, 1L, 0L), q, p, byrow = TRUE)
d <- colSums(X^2)
sparse_group_reference <- function(lambda1, lambda2, common_sigma = 1) {
  result <- matrix(0, p, q)
  for (j in seq_len(p)) {
    penalized <- C[, j] == 1L
    unpenalized <- C[, j] == 2L
    a <- A[j, penalized]
    soft <- sign(a) * pmax(abs(a) - lambda1 * common_sigma / d[j], 0)
    magnitude <- sqrt(sum(soft^2))
    if (magnitude > 0) {
      result[j, penalized] <- soft *
        max(0, 1 - lambda2 * common_sigma / (d[j] * magnitude))
    }
    result[j, unpenalized] <- A[j, unpenalized]
  }
  t(result)
}
lasso <- remMap.one(X, Y, FALSE, 3, 0, C = C)
assert_close(lasso, sparse_group_reference(3, 0),
             "DF-02 C=0/1/2 and zero group penalty")
sg <- remMap.one(X, Y, FALSE, 3, 4, C = C)
assert_close(sg, sparse_group_reference(3, 4),
             "DF-02 equal-sigma sparse-group analytic solution")
# Stationarity in the active penalized block is independent of the prox
# implementation; this fixture keeps all active entries away from zero.
for (j in seq_len(p)) {
  active <- C[, j] == 1L
  bj <- sg[active, j]
  stopifnot(all(bj != 0))
  kkt <- d[j] * (bj - A[j, active]) + 3 * sign(bj) +
    4 * bj / sqrt(sum(bj^2))
  assert_close(kkt, rep(0, length(kkt)),
               paste("DF-02 active-block KKT", j))
}
group_zero <- remMap.one(X, Y, FALSE, 3, 1000, C = C)
assert_close(group_zero, sparse_group_reference(3, 1000),
             "DF-02 inactive penalized groups retain C=2 entries")
# Internal variance parameter: common sigma rescales both penalties.
# The unequal-sigma group update has no claimed joint-prox certificate here.
scaled_sigma <- drfarm:::remMap(X, Y, 3, 4, C.m = t(C), sigma = rep(2, q))$phi
assert_close(t(scaled_sigma), sparse_group_reference(3, 4, 2),
             "DF-02 common residual-variance scaling")
sigma <- c(0.5, 1.5, 2)
weighted_lasso <- drfarm:::remMap(X, Y, 3, 0, C.m = t(C), sigma = sigma)$phi
weighted_expected <- t(sparse_group_reference(0, 0))
for (j in seq_len(p)) {
  active <- C[, j] == 1L
  weighted_expected[j, active] <- sign(A[j, active]) *
    pmax(abs(A[j, active]) - 3 * sigma[active] / d[j], 0)
}
assert_close(weighted_lasso, weighted_expected,
             "DF-02 unequal sigma valid when group penalty is zero")

# DF-01, N02. Separate raw affine units from package scale(X)/scale(Y).
# With zero penalty, Theta_standard[j,i] = A[i,j] sd(X_i)/sd(Y_j).
# Raw predictions additionally need the intercept; no intercept is returned.
Xraw <- sweep(sweep(X, 2, c(2, 0.3), "*"), 2, c(10, -3), "+")
Yraw <- sweep(sweep(Y, 2, c(5, 0.4, 2), "*"), 2, c(7, -11, 20), "+")
Xs <- scale(Xraw)
Ys <- scale(Yraw)
standard <- remMap.one(Xraw, Yraw, TRUE, 0, 0)
manual <- remMap.one(Xs, Ys, FALSE, 0, 0)
assert_close(standard, manual, "DF-01 explicit default standardization")
raw_coefficient <- sweep(sweep(t(standard), 1, attr(Xs, "scaled:scale"), "/"),
                         2, attr(Ys, "scaled:scale"), "*")
raw_intercept <- colMeans(Yraw) - drop(colMeans(Xraw) %*% raw_coefficient)
prediction_raw <- sweep(Xraw %*% raw_coefficient, 2, raw_intercept, "+")
prediction_scaled <- sweep(sweep(Xs %*% t(standard), 2,
                                attr(Ys, "scaled:scale"), "*"),
                          2, attr(Ys, "scaled:center"), "+")
assert_close(prediction_raw, prediction_scaled,
             "DF-01 coefficient/intercept back-transformation", 1e-9)

# DF-03/05, N02. Correlated predictors Cxx=[[1,.8],[.8,1]], q=3.
# Noise and supplied latent product are orthogonal to X. The exact inverse
# correction recovers the known generating coefficient from any Theta.
# Covariance is independently computed from each observation's linear weight.
Xc <- cbind(u, 0.8 * u + 0.6 * v)
Cxx <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
truth <- rbind(c(0, 0.7, -0.2), c(0.5, -0.1, 0.9))
initial <- rbind(c(0, 0.4, 0), c(0.2, 0, 0.6))
latent <- cbind(z, v * w)
loadings <- matrix(c(0.2, -0.1, 0.4, 0.3, 0.15, -0.2), q, 2)
noise_scale <- c(0.4, 0.7, 0.3)
Yc <- Xc %*% truth + latent %*% t(loadings) + outer(w, noise_scale)
delta <- truth - initial
variance <- vapply(seq_len(q), function(j) {
  n * (drop(t(delta[, j]) %*% Cxx %*% delta[, j]) + noise_scale[j]^2) /
    (n - sum(initial[, j] != 0))
}, numeric(1))
reference_entry <- function(M) {
  corrected <- initial + M %*% Cxx %*% delta
  weights <- M %*% t(Xc) / n
  se <- sqrt(outer(rowSums(weights^2), variance))
  t(2 * pnorm(-abs(corrected / se)))
}
Minv <- solve(Cxx)
assert_close(initial + Minv %*% Cxx %*% delta, truth,
             "DF-03 correlated inverse correction recovers generating truth")
pc <- entry.pvalue(Xc, Yc, t(initial), loadings, latent, Minv, FALSE)
assert_close(pc, reference_entry(Minv),
             "DF-03/05 debiasing and matching covariance")
# A nonsymmetric supplied precision estimate detects transposing M on the
# wrong side: M Cxx M' is the covariance, even when M itself is asymmetric.
M <- matrix(c(1.2, -0.3, 0.15, 0.9), 2, 2, byrow = TRUE)
pnonsym <- entry.pvalue(Xc, Yc, t(initial), loadings, latent, M, FALSE)
assert_close(pnonsym, reference_entry(M),
             "DF-03 nonsymmetric correction and covariance orientation")
angle <- 0.41
rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
assert_close(entry.pvalue(Xc, Yc, t(initial), loadings %*% rotation,
                         latent %*% rotation, M, FALSE), pnonsym,
             "DF-04/05 factor-product rotation invariance")

# DF-04. EBIC's documented K path expects scores in the K eigenbasis.
# This checks that evaluator on independently rotated inputs; it does NOT
# assert DrFARM.whole's untouched K selection path uses that basis correctly.
Psi <- c(0.8, 1.3, 0.6)
residual <- Yc - Xc %*% initial - latent %*% t(loadings)
expected_ebic <- sum(colSums(residual^2) / Psi) + n * sum(log(Psi)) +
  log(n) * sum(initial != 0) +
  2 * log(sum(choose(q, rowSums(initial != 0))))
assert_close(DrFARM.EBIC(Xc, Yc, t(initial), loadings, latent, Psi,
                         standardize = FALSE), expected_ebic,
             "DF-04 independent EBIC original basis")
for (K in list(diag(n), toeplitz(0.6^(0:(n - 1))) + diag(seq(0.1, 0.8, length.out = n)))) {
  U <- eigen(K)$vectors
  value <- DrFARM.EBIC(Xc, Yc, t(initial), loadings, t(U) %*% latent,
                       Psi, K = K, standardize = FALSE)
  assert_close(value, expected_ebic, "DF-04 EBIC consistently rotated participant basis")
}

# DF-05, P01/N02. Preserve and expose the historical TWO-SIDED Cauchy
# formula, not a standard one-sided ACAT claim. For q=1, analytic identity:
# 2*P(Cauchy <= -|cot(pi*p)|) = 2*min(p,1-p).
# Construct Gaussian z statistics via beta=z/sqrt(n-z^2), with RSS/n=beta^2+1.
n_tail <- 64
xt <- rep(c(-1, 1), n_tail / 2)
et <- rep(c(-1, -1, 1, 1), n_tail / 4)
Xt <- cbind(xt, rep(c(rep(-1, 4), rep(1, 4)), n_tail / 8))
Theta_t <- matrix(0, 1, 2)
Bt <- matrix(0, 1, 1)
Zt <- matrix(0, n_tail, 1)
for (entry_probability in c(1e-8, 0.25, 0.5, 0.75, 1 - 1e-8, 1)) {
  z_target <- -qnorm(entry_probability / 2)
  beta <- z_target / sqrt(n_tail - z_target^2)
  Yt <- matrix(xt * beta + et, n_tail, 1)
  entry <- entry.pvalue(Xt, Yt, Theta_t, Bt, Zt, diag(2), FALSE)[1, 1]
  combined <- pleio.pvalue(Xt, Yt, Theta_t, Bt, Zt, diag(2), FALSE)[1]
  assert_close(entry, entry_probability, "DF-05 constructed entry tail", 1e-10)
  assert_close(combined, 2 * min(entry_probability, 1 - entry_probability),
               "DF-05 historical single-outcome two-sided Cauchy identity", 1e-10)
}

# Invalid inputs should fail before native code or authoritative-looking
# inference. Positive degrees of freedom and variance are necessary guards,
# not sufficient conditions for inferential validity.
bad_X <- X
bad_X[1, 1] <- NA_real_
assert_error(remMap.one(bad_X, Y, FALSE, 0, 0), "DF-01 nonfinite input")
assert_error(remMap.one(X[-1, , drop = FALSE], Y, FALSE, 0, 0),
             "DF-01 mismatched participant counts")
assert_error(remMap.one(cbind(1, X[, 2]), Y, TRUE, 0, 0),
             "DF-01 constant standardized predictor")
assert_error(remMap.one(X, cbind(1, Y[, 2:3]), TRUE, 0, 0),
             "DF-01 constant standardized response")
assert_error(remMap.one(cbind(0, X[, 2]), Y, FALSE, 0, 0),
             "DF-01 zero-norm raw predictor")
assert_error(remMap.one(X, Y, FALSE, -1, 0), "DF-02 negative penalty")
assert_error(remMap.one(X, Y, FALSE, 0, 0, C = t(C)),
             "DF-01/02 transposed inclusion matrix")
assert_error(remMap.one(X, Y, FALSE, 0, 0, C = matrix(3L, q, p)),
             "DF-02 invalid inclusion code")
assert_error(entry.pvalue(Xc, Xc %*% truth, t(truth),
                          matrix(0, q, 1), matrix(0, n, 1), Minv, FALSE),
             "DF-05 zero residual variance")
assert_error(entry.pvalue(Xc, Yc, t(initial), loadings, latent,
                          matrix(0, p, p), FALSE),
             "DF-05 zero sandwich variance")
assert_error(entry.pvalue(diag(2), matrix(c(2, 3), 2, 1),
                          matrix(1, 1, 2), matrix(0, 1, 1),
                          matrix(0, 2, 1), diag(2), FALSE),
             "DF-05 nonpositive residual degrees of freedom")

# DF-06: an opt-in diagnostic must preserve the historical coefficient
# result. Native iterations count coordinate updates, not outer fit steps.
diagnostic <- remMap.one(X, Y, FALSE, 3, 4, C = C, diagnostics = TRUE)
assert_close(diagnostic$Theta0, sg, "DF-06 diagnostics preserve coefficients")
stopifnot(is.list(diagnostic$diagnostics),
          is.numeric(diagnostic$diagnostics$iterations),
          diagnostic$diagnostics$iterations >= 0,
          isTRUE(diagnostic$diagnostics$converged),
          is.character(diagnostic$diagnostics$termination),
          is.finite(diagnostic$diagnostics$final.delta),
          diagnostic$diagnostics$final.delta <= diagnostic$diagnostics$threshold)
cat("PASS DF-06 native diagnostic finite-case contract\n")

cat("All deterministic numerical contract fixtures completed.\n")
