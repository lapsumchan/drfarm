# Separately named Gaussian optimization baseline. These deterministic tests
# concern a supplied-scale, independent-participant Gaussian objective, not
# DrFARM inference, default initialization or global nonconvex optimization.
library(drfarm)

ecm_checks <- 0L
ecm_errors <- numeric()
ecm_check <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label)
  ecm_checks <<- ecm_checks + 1L
  cat("PASS", label, "\n")
}
ecm_close <- function(actual, expected, label, tolerance = 1e-9) {
  ecm_check(identical(dim(actual), dim(expected)) &&
              length(actual) == length(expected) &&
              all(is.finite(c(actual, expected))) &&
              max(abs(actual - expected)) <= tolerance, label)
  ecm_errors <<- c(ecm_errors, max(abs(actual - expected)))
}
ecm_error <- function(expr, label) {
  answer <- tryCatch({ force(expr); NULL }, error = identity)
  ecm_check(inherits(answer, "error"), label)
}

# Independent reference: condition using the q-dimensional response covariance
# rather than the latent precision. Determinant/solve and augmented QR differ
# from the production Cholesky objective and loading normal-equation solve.
ecm_penalty <- function(Theta, lambda1, lambda2, C) {
  if (any(Theta[C == 0L] != 0)) return(Inf)
  penalized <- Theta
  penalized[C != 1L] <- 0
  lambda1 * sum(abs(penalized)) +
    lambda2 * sum(sqrt(colSums(penalized^2)))
}
ecm_observed <- function(X, Y, Theta, B, psi, lambda1, lambda2, C) {
  R <- Y - X %*% t(Theta)
  Sigma <- tcrossprod(B) + diag(psi, ncol(Y))
  nrow(Y) * as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) / 2 +
    sum(R * t(solve(Sigma, t(R)))) / 2 +
    ecm_penalty(Theta, lambda1, lambda2, C)
}
ecm_posterior <- function(X, Y, Theta, B, psi) {
  Sigma <- tcrossprod(B) + diag(psi, ncol(Y))
  gain <- solve(Sigma, B)
  list(M = (Y - X %*% t(Theta)) %*% gain,
       V = diag(ncol(B)) - crossprod(B, gain))
}
ecm_full_rss <- function(X, Y, Theta, B, posterior) {
  residual <- Y - X %*% t(Theta) - posterior$M %*% t(B)
  colSums(residual^2) +
    nrow(Y) * diag(B %*% posterior$V %*% t(B))
}
ecm_Q <- function(X, Y, Theta, B, psi, posterior, lambda1, lambda2, C) {
  nrow(Y) * sum(log(psi)) / 2 +
    sum(ecm_full_rss(X, Y, Theta, B, posterior) / psi) / 2 +
    ecm_penalty(Theta, lambda1, lambda2, C)
}
ecm_factor_QR <- function(X, Y, Theta, posterior, psi.min = 0) {
  augmented_X <- rbind(posterior$M, sqrt(nrow(Y)) * chol(posterior$V))
  augmented_Y <- rbind(Y - X %*% t(Theta), matrix(0, ncol(posterior$M), ncol(Y)))
  B <- t(qr.coef(qr(augmented_X), augmented_Y))
  psi <- pmax(colSums((augmented_Y - augmented_X %*% t(B))^2) / nrow(Y), psi.min)
  list(B = B, psi = psi)
}
ecm_orthogonal_coefficient <- function(X, Y, B, psi, posterior, lambda1, C) {
  # This exact scalar lasso reference is used only for lambda2=0 and X'X
  # diagonal. It explicitly distinguishes C=0, C=1 and C=2.
  gram <- crossprod(X)
  stopifnot(max(abs(gram - diag(diag(gram), ncol(X)))) < 1e-12)
  T <- t(sweep(crossprod(X, Y - posterior$M %*% t(B)), 1L, diag(gram), "/"))
  for (j in seq_len(ncol(X))) {
    penalized <- C[, j] == 1L
    T[penalized, j] <- sign(T[penalized, j]) *
      pmax(abs(T[penalized, j]) - lambda1 * psi[penalized] / gram[j, j], 0)
  }
  T[C == 0L] <- 0
  T
}

# Adapted from the previously executed independent Gaussian fixtures. No data
# generation, seed search or psych factor initializer is used here.
ecm_asymmetric <- local({
  X <- matrix(c(-1, -1, 1, 1), 4L, 1L)
  z <- c(1, -1, 1, -1)
  u <- c(1, -1, -1, 1)
  list(X = X, Y = X %*% matrix(c(1.5, 0.5), 1L, 2L) +
         outer(z, c(1, 2)) + outer(u, c(0.3, -0.7)),
       Theta0 = matrix(0, 2L, 1L), B0 = matrix(c(1, 0.5), 2L, 1L),
       psi0 = c(1, 2), lambda1 = 1, lambda2 = 0)
})
ecm_stationary <- local({
  X <- matrix(c(-1, -1, 1, 1), 4L, 1L)
  z <- c(1, -1, 1, -1)
  u <- c(1, -1, -1, 1)
  list(X = X, Y = X %*% matrix(c(7/4, 11/4), 1L, 2L) +
         outer(z, rep(sqrt(15)/4, 2)) + outer(u, c(1, -1)/sqrt(2)),
       Theta0 = matrix(c(1, 2), 2L, 1L), B0 = matrix(1, 2L, 1L),
       psi0 = c(1, 1), lambda1 = 1, lambda2 = 0)
})

ecm_one_control <- list(max.iter = 1L, objective.tol = 1e-12, score.tol = 1e-12)
ecm_coefficient_control <- list(tol = 1e-12, root.tol = 1e-14)
ecm_run <- function(a, control = ecm_one_control, coefficient.control = ecm_coefficient_control) {
  suppressWarnings(do.call(gaussian.ecm.reference,
                           c(a, list(control = control, coefficient.control = coefficient.control))))
}

# DF-03/04/06: the earlier stationary state must stay stationary when all
# operations use its single coefficient tuple. The independent Gaussian
# score is zero at this fixed state; a historical debiasing cycle moves it.
a <- ecm_stationary
C <- matrix(1L, 2L, 1L)
Sigma <- tcrossprod(a$B0) + diag(a$psi0)
R <- a$Y - a$X %*% t(a$Theta0)
ecm_close(crossprod(R) / nrow(R), Sigma,
          "DF-04 stationary fixture residual covariance identity")
ecm_close(-crossprod(a$X, R) %*% solve(Sigma) + matrix(1, 1L, 2L),
          matrix(0, 1L, 2L), "P05 stationary fixture coefficient score")
stationary <- ecm_run(a)
for (key in c("Theta", "B", "diag.Psi")) {
  expected <- switch(key, Theta = a$Theta0, B = a$B0, diag.Psi = a$psi0)
  ecm_close(stationary[[key]], expected, paste("DF-06 initial stationary", key))
}
ecm_close(stationary$diagnostics$objective, 7 + 2 * log(3),
          "DF-06 stationary observed objective exact value")
ecm_check(isTRUE(stationary$diagnostics$converged) &&
            identical(stationary$diagnostics$termination, "initial_stationary") &&
            stationary$diagnostics$iterations == 0L &&
            stationary$diagnostics$attempted == 0L,
          "DF-06 initial stationarity explicitly returns zero cycles")

# P05/DF-02: positive group penalty at a nonzero observed-likelihood
# stationary group. Set delta=Sigma*(lambda1+lambda2*Theta/||Theta||)/n.
# Orthogonal residual contrasts then complete residual covariance to Sigma.
# Hence -X'R Sigma^{-1} exactly cancels both penalty gradients, while the
# loading/variance observed gradients vanish. This tests the new outer score,
# independently of the weighted solver's frozen-Q KKT implementation.
group_input <- local({
  X <- matrix(rep(c(-1, 1), 4L), 8L, 1L)
  Z <- cbind(rep(c(-1, -1, 1, 1), 2L), c(rep(-1, 4L), rep(1, 4L)))
  Theta0 <- matrix(c(1, 2), 2L, 1L)
  B0 <- matrix(c(1, 0.5), 2L, 1L)
  psi0 <- c(1, 2)
  Sigma <- tcrossprod(B0) + diag(psi0)
  penalty_gradient <- 0.4 + 0.7 * Theta0 / sqrt(sum(Theta0^2))
  delta <- Sigma %*% penalty_gradient / nrow(X)
  R <- X %*% t(delta) + Z %*% chol(Sigma - tcrossprod(delta))
  list(X = X, Y = X %*% t(Theta0) + R, Theta0 = Theta0,
       B0 = B0, psi0 = psi0, lambda1 = 0.4, lambda2 = 0.7)
})
group_stationary <- ecm_run(group_input)
ecm_close(group_stationary$diagnostics$gradients$coefficient,
          -0.4 - 0.7 * group_input$Theta0 / sqrt(sum(group_input$Theta0^2)),
          "P05 positive-group observed gradient balances both penalties")
ecm_check(isTRUE(group_stationary$diagnostics$converged) &&
            identical(group_stationary$diagnostics$termination, "initial_stationary") &&
            group_stationary$diagnostics$stationarity$maximum <= 1e-12,
          "P05 positive-group observed subgradient stationarity recognized")

# DF-02/03/06: one asymmetric cycle, using an exact orthogonal lasso
# coefficient followed by independent augmented QR and full posterior RSS.
a <- ecm_asymmetric
C <- matrix(1L, 2L, 1L)
posterior <- ecm_posterior(a$X, a$Y, a$Theta0, a$B0, a$psi0)
expected_Theta <- ecm_orthogonal_coefficient(a$X, a$Y, a$B0, a$psi0,
                                            posterior, a$lambda1, C)
ecm_close(expected_Theta, matrix(c(33/68, 0), 2L, 1L),
          "DF-02 rederived asymmetric scalar lasso reference")
expected_factor <- ecm_factor_QR(a$X, a$Y, expected_Theta, posterior)
asymmetric <- ecm_run(a)
ecm_close(asymmetric$Theta, expected_Theta, "DF-02 one-cycle coefficient CM")
ecm_close(asymmetric$B, expected_factor$B, "DF-04 one-cycle loading CM versus augmented QR")
ecm_close(asymmetric$diag.Psi, expected_factor$psi,
          "DF-03 one-cycle variance CM versus augmented residual sums")
ecm_check(any(abs(asymmetric$Theta - matrix(c(25/34, 2/17), 2L, 1L)) > 0.1),
          "DF-03 coefficient remains sparse rather than receiving inner debiasing")
expected_Q <- c(
  initial = ecm_Q(a$X, a$Y, a$Theta0, a$B0, a$psi0, posterior, 1, 0, C),
  coefficient = ecm_Q(a$X, a$Y, expected_Theta, a$B0, a$psi0, posterior, 1, 0, C),
  loadings = ecm_Q(a$X, a$Y, expected_Theta, expected_factor$B, a$psi0, posterior, 1, 0, C),
  variance = ecm_Q(a$X, a$Y, expected_Theta, expected_factor$B,
                    expected_factor$psi, posterior, 1, 0, C))
ecm_close(unname(asymmetric$diagnostics$trials[[1L]]$q), unname(expected_Q),
          "DF-06 four frozen-Q values independently recomputed")
ecm_check(all(diff(expected_Q) < -1e-4), "DF-06 each asymmetric CM decreases frozen Q")
expected_L <- ecm_observed(a$X, a$Y, asymmetric$Theta, asymmetric$B,
                           asymmetric$diag.Psi, 1, 0, C)
ecm_close(asymmetric$diagnostics$objective, expected_L,
          "DF-06 fresh observed objective equals independent dense likelihood")
expected_initial_L <- ecm_observed(a$X, a$Y, a$Theta0, a$B0, a$psi0, 1, 0, C)
ecm_close(asymmetric$diagnostics$trials[[1L]]$observed.before, expected_initial_L,
          "DF-06 trial begins at independently evaluated accepted objective")
ecm_close(asymmetric$diagnostics$trials[[1L]]$observed.after, expected_L,
          "DF-06 trial acceptance uses independently evaluated candidate objective")
ecm_close(asymmetric$diagnostics$history$objective, c(expected_initial_L, expected_L),
          "DF-06 history contains initial and accepted observed likelihoods")
fresh <- ecm_posterior(a$X, a$Y, asymmetric$Theta, asymmetric$B, asymmetric$diag.Psi)
ecm_close(asymmetric$E.Z, fresh$M, "DF-03 returned posterior mean uses final accepted tuple")
ecm_close(asymmetric$posterior.covariance, fresh$V,
          "DF-03 returned posterior covariance uses final accepted tuple")
ecm_check(max(abs(asymmetric$E.Z - posterior$M)) > 0.1,
          "DF-03 asymmetric posterior distinguishes fresh from frozen moments")
ecm_check(!asymmetric$diagnostics$converged &&
            identical(asymmetric$diagnostics$termination, "max_iter") &&
            asymmetric$diagnostics$iterations == 1L &&
            asymmetric$diagnostics$attempted == 1L,
          "DF-06 finite outer budget retains explicit nonconvergence")

# DF-01/02: nonsquare matrices and C=0/1/2 within a shared predictor group.
# The factor/variance QR reference continues to use the same masked coefficient.
u <- rep(c(-1, 1), 4L)
v <- rep(c(-1, -1, 1, 1), 2L)
w <- c(rep(-1, 4L), rep(1, 4L))
mixed <- list(X = cbind(u, 2 * v),
              Y = cbind(u, 2 * v) %*% rbind(c(1.5, -0.7, 2), c(-0.5, 1.2, 0.35)) +
                outer(w, c(0.5, -0.25, 0.75)),
              Theta0 = matrix(0, 3L, 2L), B0 = matrix(c(1, 0.5, -0.25), 3L, 1L),
              psi0 = c(0.7, 1.3, 2), lambda1 = 0.8, lambda2 = 0,
              C = matrix(c(0L, 1L, 2L, 1L, 2L, 1L), 3L, 2L))
a <- mixed
post_mixed <- ecm_posterior(a$X, a$Y, a$Theta0, a$B0, a$psi0)
theta_mixed <- ecm_orthogonal_coefficient(a$X, a$Y, a$B0, a$psi0,
                                         post_mixed, a$lambda1, a$C)
factor_mixed <- ecm_factor_QR(a$X, a$Y, theta_mixed, post_mixed)
masked <- ecm_run(a)
ecm_close(masked$Theta, theta_mixed, "DF-01/02 nonsquare C0/C1/C2 coefficient CM")
ecm_close(masked$B, factor_mixed$B, "DF-03 masked loading update uses identical coefficient")
ecm_close(masked$diag.Psi, factor_mixed$psi,
          "DF-03 masked variance update uses identical coefficient")
ecm_check(all(masked$Theta[a$C == 0L] == 0) &&
            all(abs(masked$Theta[a$C == 2L]) > 0.1),
          "DF-02 exclusion feasibility and nonzero unpenalized coefficients")
ecm_close(masked$diagnostics$objective,
          ecm_observed(a$X, a$Y, masked$Theta, masked$B, masked$diag.Psi,
                         a$lambda1, a$lambda2, a$C),
          "DF-02 observed monitor applies the mask-aware penalty")

# P04/N02: returned arrays are gradients of the smooth observed likelihood.
# Compare to central finite differences of the independent dense likelihood,
# allowing all coefficient coordinates to vary and setting penalties to zero.
# Predetermined step=1e-6 and absolute tolerance=2e-7 account for subtraction
# roundoff on this small well-scaled fixture; these are not optimizer bounds.
derivative_step <- 1e-6
derivative_tolerance <- 2e-7
gradient_point <- list(Theta = masked$Theta, B = masked$B, psi = masked$diag.Psi)
likelihood_at <- function(point) {
  ecm_observed(mixed$X, mixed$Y, point$Theta, point$B, point$psi,
                 0, 0, matrix(1L, 3L, 2L))
}
for (component in c("Theta", "B", "psi")) {
  difference <- gradient_point[[component]]
  for (j in seq_along(difference)) {
    up <- down <- gradient_point
    up[[component]][j] <- up[[component]][j] + derivative_step
    down[[component]][j] <- down[[component]][j] - derivative_step
    difference[j] <- (likelihood_at(up) - likelihood_at(down)) / (2 * derivative_step)
  }
  field <- switch(component, Theta = "coefficient", B = "loadings", psi = "variance")
  ecm_close(masked$diagnostics$gradients[[field]], difference,
            paste("P04 independent central derivative", field), derivative_tolerance)
}

# DF-04/N03: factor coordinates are unidentified. Compare fitted latent
# products and implied response covariance under a nontrivial k=2 rotation.
rotated_input <- mixed
rotated_input$B0 <- rbind(c(0.8, 0.2), c(0.3, -0.5), c(-0.2, 0.7))
rotated_input$lambda2 <- 0.3
rotation <- matrix(c(0.6, 0.8, -0.8, 0.6), 2L, 2L)
rotation_control <- list(max.iter = 3L, objective.tol = 1e-12, score.tol = 1e-12)
rotation_base <- ecm_run(rotated_input, rotation_control)
rotated_input$B0 <- rotated_input$B0 %*% rotation
rotation_fit <- ecm_run(rotated_input, rotation_control)
ecm_close(rotation_fit$Theta, rotation_base$Theta, "DF-04 factor rotation preserves coefficient")
ecm_close(rotation_fit$diag.Psi, rotation_base$diag.Psi, "DF-04 factor rotation preserves variance")
ecm_close(tcrossprod(rotation_fit$B), tcrossprod(rotation_base$B),
          "DF-04 factor rotation preserves loading covariance product")
ecm_close(rotation_fit$E.Z %*% t(rotation_fit$B),
          rotation_base$E.Z %*% t(rotation_base$B),
          "DF-04 factor rotation preserves posterior latent fitted product")
ecm_close(rotation_fit$diagnostics$history$objective,
          rotation_base$diagnostics$history$objective,
          "DF-04 factor rotation preserves accepted objective trajectory")
ecm_check(min(eigen(tcrossprod(rotation_fit$B) + diag(rotation_fit$diag.Psi),
                    symmetric = TRUE, only.values = TRUE)$values) > 0 &&
            min(eigen(rotation_fit$posterior.covariance,
                      symmetric = TRUE, only.values = TRUE)$values) > 0,
          "N02 implied response and posterior covariances are positive definite")
ecm_check(all(diff(rotation_fit$diagnostics$history$objective) <= 1e-10),
          "DF-06 accepted observed-likelihood history is nonincreasing")

# DF-06/N01: an unfinished coefficient solve cannot be promoted to an outer
# iterate. The initial accepted tuple and its freshly computed posterior survive.
inner_short <- mixed
inner_short$X <- cbind(u, 0.55 * u + sqrt(1 - 0.55^2) * v)
inner_short$lambda2 <- 0.3
limited <- ecm_run(inner_short,
                   coefficient.control = list(tol = 1e-13, max.sweeps = 1L))
ecm_check(!limited$diagnostics$converged &&
            identical(limited$diagnostics$termination, "coefficient_failure") &&
            limited$diagnostics$iterations == 0L && limited$diagnostics$attempted == 1L,
          "DF-06 finite inner budget rejects the unfinished outer trial")
ecm_close(limited$Theta, inner_short$Theta0, "DF-06 inner failure rolls back coefficient")
ecm_close(limited$B, inner_short$B0, "DF-06 inner failure rolls back loadings")
ecm_close(limited$diag.Psi, inner_short$psi0, "DF-06 inner failure rolls back variances")
ecm_check(!limited$diagnostics$trials[[1L]]$accepted &&
            !limited$diagnostics$trials[[1L]]$inner$converged,
          "DF-06 rejected trial retains unfinished coefficient diagnostic")

# Controller fault injection only: a local function copy receives a bogus
# finite coefficient despite a reported inner success. This does not mutate
# the installed namespace or supply evidence of natural failure frequency.
bad_coefficient_reference <- gaussian.ecm.reference
fault_environment <- new.env(parent = environment(bad_coefficient_reference))
original_coefficient_solver <- get("remMap.weighted", envir = fault_environment)
fault_environment$remMap.weighted <- function(...) {
  answer <- original_coefficient_solver(...)
  answer$Theta0 <- answer$Theta0 + 100
  answer
}
environment(bad_coefficient_reference) <- fault_environment
rejected <- suppressWarnings(do.call(bad_coefficient_reference,
                                      c(ecm_asymmetric, list(control = ecm_one_control))))
ecm_check(!rejected$diagnostics$converged &&
            identical(rejected$diagnostics$termination, "q_increase") &&
            rejected$diagnostics$iterations == 0L &&
            !rejected$diagnostics$trials[[1L]]$accepted,
          "DF-06 injected uphill coefficient is rejected by frozen-Q controller")
ecm_close(rejected$Theta, ecm_asymmetric$Theta0,
          "DF-06 injected trial rejection retains accepted coefficient")
ecm_close(rejected$B, ecm_asymmetric$B0,
          "DF-06 injected trial rejection retains accepted loading")
ecm_close(rejected$diag.Psi, ecm_asymmetric$psi0,
          "DF-06 injected trial rejection retains accepted variance")

# N02/DF-06: the unrestricted positive-variance domain has no finite Gaussian
# minimizer when the same-tuple fitted residual vanishes exactly. A user-selected
# positive floor is an explicit constrained objective and has a valid CM.
boundary_input <- list(X = matrix(c(-1, 1), 2L, 1L),
                        Y = matrix(c(-1, 1), 2L, 1L) %*% matrix(c(1, 2), 1L, 2L),
                        Theta0 = matrix(0, 2L, 1L), B0 = matrix(0, 2L, 1L),
                        psi0 = c(1, 1), lambda1 = 0, lambda2 = 0)
boundary <- ecm_run(boundary_input)
ecm_check(!boundary$diagnostics$converged &&
            identical(boundary$diagnostics$termination, "variance_boundary") &&
            boundary$diagnostics$iterations == 0L,
          "N02 zero residual variance returns explicit open-domain boundary failure")
ecm_close(boundary$Theta, boundary_input$Theta0, "DF-06 boundary failure rolls back complete tuple")
floor_fit <- ecm_run(boundary_input,
                     list(max.iter = 3L, objective.tol = 1e-12,
                          score.tol = 1e-12, psi.min = c(0.25, 0.5)))
ecm_close(floor_fit$Theta, matrix(c(1, 2), 2L, 1L), "DF-02 constrained exact-fit coefficients")
ecm_close(floor_fit$diag.Psi, c(0.25, 0.5), "N02 explicit vector variance floor attained")
ecm_close(floor_fit$diagnostics$objective, log(0.25) + log(0.5),
          "DF-06 constrained exact-fit observed objective")
ecm_check(isTRUE(floor_fit$diagnostics$converged) &&
            floor_fit$diagnostics$stationarity$maximum <= 1e-12,
          "P05 active variance-floor KKT has the correct one-sided sign")

# Public scope and input contracts. No precision matrix, debiased coefficient,
# or inferential result is offered by this separately named optimization API.
ecm_check(identical(asymmetric$diagnostics$inference.status, "not_provided") &&
            identical(asymmetric$diagnostics$method, "Gaussian ECM optimization reference") &&
            !any(c("Theta.db", "p.values", "p.value", "precM") %in% names(asymmetric)),
          "DF-03 explicit optimization-only method and inference status")
invalid <- ecm_asymmetric
invalid$psi0 <- c(1, 0)
ecm_error(ecm_run(invalid), "N02 nonpositive starting variance rejected")
invalid <- ecm_asymmetric
invalid$Theta0 <- t(invalid$Theta0)
ecm_error(ecm_run(invalid), "DF-01 transposed coefficient start rejected")
invalid <- mixed
invalid$Theta0[invalid$C == 0L] <- 1
ecm_error(ecm_run(invalid), "DF-02 infeasible C0 start rejected without projection")
ecm_error(ecm_run(ecm_asymmetric, list(max.iter = 0L)),
          "DF-06 zero iteration budget rejected")
ecm_error(ecm_run(ecm_asymmetric, list(psi.min = c(2, 2))),
          "N02 initial tuple below explicit variance floor rejected")
ecm_error(do.call(gaussian.ecm.reference, c(ecm_asymmetric, list(K = diag(4)))),
          "DF-04 optional kinship is outside the reference API")
cat("Gaussian ECM reference fixtures completed:", ecm_checks, "checks.\n")
cat("Largest independent numeric discrepancy:", format(max(ecm_errors), digits = 16),
    "; algebra/trajectory tolerance: 1e-9; central derivative tolerance: 2e-7 at step 1e-6.\n")
