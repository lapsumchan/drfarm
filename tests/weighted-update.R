# Fixed-variance convex coefficient update: independent references and KKT.
# Inputs, sigma (variances), penalties and Theta0 all use the supplied scale.
# There is no automatic centering/standardization or intercept in this API.
library(drfarm)

w_close <- function(actual, expected, label, tol = 1e-7) {
  stopifnot(identical(dim(actual), dim(expected)), length(actual) == length(expected),
            all(is.finite(actual)), all(is.finite(expected)))
  err <- max(abs(actual - expected))
  if (err > tol) stop(label, ": max absolute discrepancy ", err, " > ", tol)
  cat("PASS", label, "; max absolute discrepancy =", format(err), "\n")
}
w_error <- function(expr, label) {
  value <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(value, "error")) stop(label, ": expected an error")
  cat("PASS", label, "; message =", conditionMessage(value), "\n")
}
w_objective <- function(T, X, Y, sigma, C, lambda1, lambda2) {
  residual <- Y - X %*% t(T)
  groups <- vapply(seq_len(ncol(T)), function(j) {
    sqrt(sum(T[C[, j] == 1L, j]^2))
  }, numeric(1))
  sum(residual * sweep(residual, 2L, sigma, "/")) / 2 +
    lambda1 * sum(abs(T[C == 1L])) + lambda2 * sum(groups)
}

# Global KKT residual computed from the complete model residual, independently
# of block root/stopping code. Active groups use coordinate subgradient
# distances; zero groups use dist(-gradient, lambda1*[-1,1]^m)<=lambda2.
# C=2 is unpenalized and C=0 is constrained, including within mixed groups.
w_kkt <- function(T, X, Y, sigma, C, lambda1, lambda2) {
  G <- -t(crossprod(X, sweep(Y - X %*% t(T), 2L, sigma, "/")))
  violations <- c(0, abs(T[C == 0L]), abs(G[C == 2L]))
  for (j in seq_len(ncol(T))) {
    active <- C[, j] == 1L
    b <- T[active, j]
    g <- G[active, j]
    norm_b <- sqrt(sum(b^2))
    if (norm_b == 0) {
      violations <- c(violations,
                      max(sqrt(sum(pmax(abs(g) - lambda1, 0)^2)) - lambda2, 0))
    } else {
      r <- g + lambda2 * b / norm_b
      violations <- c(violations,
                      abs(r[b != 0] + lambda1 * sign(b[b != 0])),
                      pmax(abs(r[b == 0]) - lambda1, 0))
    }
  }
  max(violations)
}

# Independent full-matrix BFGS reference. Only used on fixtures whose free
# coordinates and penalized group norms are nonzero at the solution, so the
# objective is differentiable locally. Full-rank X gives strict convexity:
# a small independently checked global KKT residual identifies the unique fit.
w_reference <- function(X, Y, sigma, C, lambda1, lambda2, start) {
  q <- ncol(Y)
  p <- ncol(X)
  free <- C != 0L
  unpack <- function(par) { T <- matrix(0, q, p); T[free] <- par; T }
  fn <- function(par) w_objective(unpack(par), X, Y, sigma, C, lambda1, lambda2)
  gr <- function(par) {
    T <- unpack(par)
    G <- -t(crossprod(X, sweep(Y - X %*% t(T), 2L, sigma, "/")))
    G[C == 1L] <- G[C == 1L] + lambda1 * sign(T[C == 1L])
    for (j in seq_len(p)) {
      active <- C[, j] == 1L
      norm_b <- sqrt(sum(T[active, j]^2))
      if (norm_b > 0) G[active, j] <- G[active, j] + lambda2 * T[active, j] / norm_b
    }
    G[free]
  }
  fit <- optim(start[free], fn, gr, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 5000L))
  T <- unpack(fit$par)
  stopifnot(fit$convergence == 0, all(abs(T[free]) > 1e-4),
            w_kkt(T, X, Y, sigma, C, lambda1, lambda2) < 2e-6)
  T
}

ctl <- list(tol = 1e-11, max.sweeps = 1000L, root.tol = 1e-14, root.maxit = 200L)

# DF-02, P05/P06, N01/N02. Re-derived exact heterogeneous-variance case:
# d=1, h=(1,4), score=(8,22), soft score=(6,20), lambda2=5.
# beta=(3,4), norm(beta)=5 gives t=lambda2/norm(beta)=1 and zero KKT.
X1 <- matrix(c(1, -1) / sqrt(2), 2L, 1L)
Y1 <- X1 %*% matrix(c(8, 5.5), 1L, 2L)
sigma1 <- c(1, 0.25)
C1 <- matrix(1L, 2L, 1L)
exact <- remMap.weighted(X1, Y1, lambda1 = 2, lambda2 = 5,
                         sigma = sigma1, control = ctl)
w_close(exact$Theta0, matrix(c(3, 4), 2L, 1L),
        "DF-02 exact unequal-variance nonzero block", 1e-9)
stopifnot(isTRUE(exact$diagnostics$converged),
          w_kkt(exact$Theta0, X1, Y1, sigma1, C1, 2, 5) < 1e-8)
cat("PASS DF-02 exact block global KKT\n")

# A zero entry inside an active group obeys an L1 subgradient interval,
# rather than the nonzero-coordinate sign equation or whole-group zero test.
Yactive_zero <- X1 %*% matrix(c(8, 5.5, 2), 1L, 3L)
active_zero <- remMap.weighted(X1, Yactive_zero, 2, 5,
                               sigma = c(1, 0.25, 2), control = ctl)
w_close(active_zero$Theta0, matrix(c(3, 4, 0), 3L, 1L),
        "DF-02 zero entry inside active unequal-variance group", 1e-9)
stopifnot(w_kkt(active_zero$Theta0, X1, Yactive_zero, c(1, 0.25, 2),
                matrix(1L, 3L, 1L), 2, 5) < 1e-8)

# Previously executed historical counterexample: new fit must improve its
# chosen weighted objective and match a separate two-variable optimizer.
Yold <- X1 %*% matrix(c(3, 4), 1L, 2L)
sigma_old <- c(1, 2)
corrected <- remMap.weighted(X1, Yold, 0, 1, sigma = sigma_old, control = ctl)
reference_old <- w_reference(X1, Yold, sigma_old, C1, 0, 1,
                             matrix(c(3, 4), 2L, 1L))
w_close(corrected$Theta0, reference_old,
        "DF-02 historical counterexample versus independent BFGS", 1e-6)
historical <- matrix(c(2.4, 2.4), 2L, 1L)
stopifnot(w_objective(corrected$Theta0, X1, Yold, sigma_old, C1, 0, 1) <
            w_objective(historical, X1, Yold, sigma_old, C1, 0, 1) - 1e-5,
          w_kkt(corrected$Theta0, X1, Yold, sigma_old, C1, 0, 1) < 1e-8)
cat("PASS DF-02 corrected weighted objective improves historical counterexample\n")

# Exact zero-group boundary after L1: score=(4,5), soft=(3,4), norm=5.
# X=(1,0) keeps this boundary exactly representable in floating point.
Xzero <- matrix(c(1, 0), 2L, 1L)
Yzero <- matrix(c(4, 0, 10, 0), 2L, 2L)
at_zero <- remMap.weighted(Xzero, Yzero, 1, 5, sigma = c(1, 2), control = ctl)
w_close(at_zero$Theta0, matrix(0, 2L, 1L), "DF-02 zero-group equality boundary", 0)
above_zero <- remMap.weighted(Xzero, Yzero, 1, 5.01, sigma = c(1, 2), control = ctl)
w_close(above_zero$Theta0, matrix(0, 2L, 1L), "DF-02 inactive group beyond boundary", 0)
below_zero <- remMap.weighted(Xzero, Yzero, 1, 4.99, sigma = c(1, 2), control = ctl)
stopifnot(sqrt(sum(below_zero$Theta0^2)) > 0,
          w_kkt(below_zero$Theta0, Xzero, Yzero, c(1, 2), C1, 1, 4.99) < 1e-8)
cat("PASS DF-02 active group below boundary\n")

# DF-01/02, P05/P06. Orthogonal nonsquare design, all mask values,
# heterogeneous sigma, and zero group penalty: independent scalar lasso.
u <- rep(c(-1, 1), 4)
v <- rep(c(-1, -1, 1, 1), 2)
w <- c(rep(-1, 4), rep(1, 4))
X <- cbind(u, 2 * v)
A <- rbind(c(1.5, -0.7, 2), c(-0.5, 1.2, 0.35))
Y <- X %*% A + outer(w, c(0.5, -0.25, 0.75))
C <- matrix(c(1L, 2L, 0L, 1L, 1L, 0L), 3L, 2L, byrow = TRUE)
sigma <- c(0.5, 1.4, 2.3)
expected <- t(A)
expected[C == 0L] <- 0
for (j in seq_len(ncol(X))) {
  selected <- C[, j] == 1L
  expected[selected, j] <- sign(A[j, selected]) *
    pmax(abs(A[j, selected]) - 3 * sigma[selected] / sum(X[, j]^2), 0)
}
lasso <- remMap.weighted(X, Y, 3, 0, sigma, C, control = ctl)
w_close(lasso$Theta0, expected, "DF-01/02 masked heterogeneous-variance lasso", 1e-9)
ols <- remMap.weighted(X, Y, 0, 0, sigma = sigma, control = ctl)
w_close(ols$Theta0, t(qr.coef(qr(X), Y)), "DF-01/02 zero-penalty QR coefficient orientation", 1e-9)
equal <- remMap.weighted(X, Y, 3, 4, sigma = rep(2, 3), C = C, control = ctl)
legacy_equal <- drfarm:::remMap(X, Y, 3, 4, C.m = t(C), sigma = rep(2, 3))$phi
w_close(equal$Theta0, t(legacy_equal), "DF-02 equal-variance orthogonal historical agreement", 1e-9)

# DF-02/06; P05/P09/N01. Strictly convex correlated, heterogeneous global
# fit compared with full-matrix BFGS, then checked with full-residual KKT.
Xc <- cbind(u, 0.55 * u + sqrt(1 - 0.55^2) * v)
Yc <- Xc %*% A + outer(w, c(0.5, -0.25, 0.75))
Call <- matrix(1L, 3L, 2L)
global <- remMap.weighted(Xc, Yc, 0.3, 0.7, sigma, control = ctl)
global_ref <- w_reference(Xc, Yc, sigma, Call, 0.3, 0.7, t(A))
w_close(global$Theta0, global_ref, "DF-02 correlated global BFGS reference", 1e-6)
stopifnot(isTRUE(global$diagnostics$converged),
          w_kkt(global$Theta0, Xc, Yc, sigma, Call, 0.3, 0.7) < 1e-7)
cat("PASS DF-02 correlated global KKT\n")
masked <- remMap.weighted(Xc, Yc, 0.3, 0.7, sigma, C, control = ctl)
masked_ref <- w_reference(Xc, Yc, sigma, C, 0.3, 0.7, t(A))
w_close(masked$Theta0, masked_ref, "DF-02 mixed-mask correlated BFGS reference", 1e-6)
stopifnot(w_kkt(masked$Theta0, Xc, Yc, sigma, C, 0.3, 0.7) < 1e-7)
cat("PASS DF-02 mixed-mask global KKT\n")
w_close(global$diagnostics$objective,
        w_objective(global$Theta0, Xc, Yc, sigma, Call, 0.3, 0.7),
        "DF-06 independently recomputed full objective", 1e-9)
stopifnot(all(diff(global$diagnostics$objective.trace) <=
                1e-10 * max(1, abs(global$diagnostics$initial.objective))))
cat("PASS DF-06 objective trace is nonincreasing to recorded tolerance\n")

# DF-01/N02/N03. Only common-unit changes preserve this isotropic group
# penalty with scalar lambdas. No arbitrary outcome-wise scaling equivalence.
response_scale <- 4
rescaled_y <- remMap.weighted(Xc, Yc * response_scale,
                              0.3 / response_scale, 0.7 / response_scale,
                              sigma * response_scale^2, control = ctl)
w_close(rescaled_y$Theta0, global$Theta0 * response_scale,
        "DF-01 response/variance/penalty unit transformation", 1e-6)
w_close(rescaled_y$diagnostics$objective, global$diagnostics$objective,
        "DF-01 response-unit objective identity", 1e-9)
predictor_scale <- 3
rescaled_x <- remMap.weighted(Xc * predictor_scale, Yc,
                              0.3 * predictor_scale, 0.7 * predictor_scale,
                              sigma, control = ctl)
w_close(rescaled_x$Theta0, global$Theta0 / predictor_scale,
        "DF-01 predictor/penalty unit transformation", 1e-7)
order <- rev(seq_len(nrow(Xc)))
reordered <- remMap.weighted(Xc[order, ], Yc[order, ], 0.3, 0.7, sigma, control = ctl)
w_close(reordered$Theta0, global$Theta0, "N03 matched participant reordering", 1e-9)

# Different feasible warm starts must converge to the same unique optimum.
warm <- remMap.weighted(Xc, Yc, 0.3, 0.7, sigma,
                        Theta0 = global$Theta0 + matrix(c(0.1, -0.2, 0.3), 3L, 2L),
                        control = ctl)
w_close(warm$Theta0, global$Theta0, "DF-06 warm-start unique-solution agreement", 1e-7)
stopifnot(isTRUE(warm$diagnostics$converged),
          w_kkt(warm$Theta0, Xc, Yc, sigma, Call, 0.3, 0.7) < 1e-7)
cat("PASS DF-06 warm-start global KKT\n")

# Preserve bounded failures as finite candidates with explicit diagnostics.
short <- remMap.weighted(Xc, Yc, 0.3, 0.7, sigma,
                         control = list(tol = 1e-13, max.sweeps = 1L))
stopifnot(isFALSE(short$diagnostics$converged),
          identical(short$diagnostics$termination, "max_sweeps"),
          short$diagnostics$sweeps == 1L, all(is.finite(short$Theta0)),
          is.finite(short$diagnostics$objective),
          short$diagnostics$objective <= short$diagnostics$initial.objective)
cat("PASS DF-06 finite sweep budget retains explicit nonconvergence\n")
root_limited <- remMap.weighted(X1, Y1, 2, 5, sigma1,
                                control = list(root.maxit = 1L, root.tol = 1e-14))
stopifnot(isFALSE(root_limited$diagnostics$converged),
          identical(root_limited$diagnostics$termination, "root_failure"),
          all(is.finite(root_limited$Theta0)),
          nzchar(root_limited$diagnostics$failure))
cat("PASS DF-06 finite root budget retains explicit failure\n")

w_error(remMap.weighted(X, Y, 1, 1, sigma = c(1, 0, 1)), "N02 nonpositive variance rejected")
w_error(remMap.weighted(X, Y, 1, 1, sigma = c(1, Inf, 1)), "N02 nonfinite variance rejected")
w_error(remMap.weighted(X, Y, 1, 1, sigma = c(1, 2)), "DF-01 variance length rejected")
w_error(remMap.weighted(X, Y, 1, 1, Theta0 = matrix(0, 2, 3)), "DF-01 transposed warm start rejected")
w_error(remMap.weighted(X, Y, 1, 1, control = list(max.sweeps = 0)), "DF-06 invalid sweep budget rejected")

# DF-06/N03. Bounded outer integration on the existing bundled data. Use one
# outer step, an explicit precision matrix and a fixed initial coefficient;
# no tuning grid, precision search or scientific inference is needed here.
data("drfarm.dat", package = "drfarm", envir = environment())
Xb <- scale(drfarm.dat$X)
Yb <- scale(drfarm.dat$Y)
Tb <- matrix(0, ncol(Yb), ncol(Xb))
Mb <- solve(crossprod(Xb) / nrow(Xb))
set.seed(20260909)
outer_default <- DrFARM.one(Xb, Yb, Tb, Mb, k = 2, lambda1 = 1, lambda2 = 2,
                            standardize = FALSE, max.iter = 1)
set.seed(20260909)
outer_historical <- DrFARM.one(Xb, Yb, Tb, Mb, k = 2, lambda1 = 1, lambda2 = 2,
                               standardize = FALSE, max.iter = 1,
                               coefficient.update = "historical")
for (field in c("Theta", "B", "E.Z", "diag.Psi")) {
  stopifnot(identical(outer_default[[field]], outer_historical[[field]]))
}
cat("PASS DF-06 default and explicit historical routing: four arrays exactly equal\n")
set.seed(20260909)
outer_weighted <- DrFARM.one(Xb, Yb, Tb, Mb, k = 2, lambda1 = 1, lambda2 = 2,
                             standardize = FALSE, max.iter = 1,
                             coefficient.update = "weighted",
                             weighted.control = list(tol = 1e-9))
inner_history <- outer_weighted$diagnostics$coefficient.history
stopifnot(identical(outer_weighted$diagnostics$coefficient.update, "weighted"),
          identical(outer_weighted$diagnostics$inference.status, "unvalidated_for_weighted_update"),
          isFALSE(outer_weighted$diagnostics$converged),
          length(inner_history) == 1L,
          isTRUE(inner_history[[1]]$converged),
          inner_history[[1]]$kkt.scaled <= inner_history[[1]]$threshold,
          inner_history[[1]]$objective <= inner_history[[1]]$initial.objective,
          all(is.finite(outer_weighted$Theta)))
cat("PASS DF-06 weighted outer routing preserves coefficient history and inference status\n")
set.seed(20260909)
outer_failed <- tryCatch(DrFARM.one(Xb, Yb, Tb, Mb, k = 2, lambda1 = 1, lambda2 = 2,
                                    standardize = FALSE, max.iter = 1,
                                    coefficient.update = "weighted",
                                    weighted.control = list(tol = 1e-13, max.sweeps = 1L)),
                          error = identity)
stopifnot(inherits(outer_failed, "error"),
          grepl("Weighted coefficient subproblem failed", conditionMessage(outer_failed), fixed = TRUE))
cat("PASS DF-06 failed weighted inner solve stops the outer fit; message =",
    conditionMessage(outer_failed), "\n")
cat("Weighted fixed-variance coefficient fixtures completed.\n")
