# This separately named optimization reference does not call inner debiasing,
# precision estimation, or the DrFARM inference functions.
.ecm_control <- function(control, q) {
  defaults <- list(max.iter = 500L, objective.tol = 1e-8,
                   score.tol = 1e-6, psi.min = 0)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% names(defaults)) ||
       anyDuplicated(names(control)))))
    stop("control must be a named list of max.iter, objective.tol, score.tol, psi.min.", call. = FALSE)
  defaults[names(control)] <- control
  for (key in c("max.iter", "objective.tol", "score.tol")) {
    x <- defaults[[key]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
      stop(key, " must be a finite numeric scalar.", call. = FALSE)
  }
  if (defaults$max.iter <= 0 || defaults$max.iter != floor(defaults$max.iter) ||
      defaults$max.iter > .Machine$integer.max)
    stop("max.iter must be a positive integer within the R integer range.", call. = FALSE)
  if (defaults$objective.tol < 0 || defaults$objective.tol >= 1)
    stop("objective.tol must be nonnegative and smaller than one.", call. = FALSE)
  if (defaults$score.tol <= 0)
    stop("score.tol must be positive.", call. = FALSE)
  x <- defaults$psi.min
  if (!is.numeric(x) || !length(x) %in% c(1L, q) ||
      any(!is.finite(x)) || any(x < 0))
    stop("psi.min must be nonnegative and finite, of length one or q.", call. = FALSE)
  defaults$psi.min <- rep(x, length.out = q)
  defaults
}

.ecm_penalty <- function(Theta, C, lambda1, lambda2) {
  value <- lambda1 * sum(abs(Theta[C == 1L]))
  for (j in seq_len(ncol(Theta)))
    value <- value + lambda2 * .weighted_norm(Theta[C[, j] == 1L, j])
  if (!is.finite(value)) stop("Nonfinite ECM penalty.", call. = FALSE)
  value
}

# Minimum subgradient distance for each predictor group, with the gradient
# of the *observed* full-covariance likelihood, not the frozen-Q gradient.
.ecm_coefficient_score <- function(gradient, Theta, C, lambda1, lambda2) {
  residuals <- numeric(ncol(Theta))
  for (j in seq_len(ncol(Theta))) {
    pen <- C[, j] == 1L
    b <- Theta[pen, j]
    g <- gradient[j, pen]
    radius <- .weighted_norm(b)
    if (radius == 0) {
      distance <- max(0, .weighted_norm(sign(g) * pmax(abs(g) - lambda1, 0)) - lambda2)
    } else {
      v <- g + lambda2 * b / radius
      v <- ifelse(b != 0, v + lambda1 * sign(b),
                  sign(v) * pmax(abs(v) - lambda1, 0))
      distance <- .weighted_norm(v)
    }
    residuals[j] <- .weighted_norm(c(distance, gradient[j, C[, j] == 2L]))
  }
  max(residuals)
}

.ecm_expected_rss <- function(R, B, posterior) {
  # Mean residual plus posterior uncertainty. The sum-of-squares form avoids
  # subtracting nearly equal quadratic terms in the expanded RSS identity.
  mean.residual <- R - tcrossprod(posterior$M, B)
  covariance.root <- B %*% t(chol(posterior$V))
  value <- colSums(mean.residual^2) + nrow(R) * rowSums(covariance.root^2)
  if (any(!is.finite(value)) || any(value < 0))
    stop("Invalid full posterior expected residual sums of squares.", call. = FALSE)
  value
}

.ecm_Q <- function(X, Y, Theta, B, psi, posterior, C, lambda1, lambda2) {
  value <- nrow(Y) * sum(log(psi)) / 2 +
    sum(.ecm_expected_rss(Y - X %*% t(Theta), B, posterior) / psi) / 2 +
    .ecm_penalty(Theta, C, lambda1, lambda2)
  if (!is.finite(value)) stop("Nonfinite ECM Q objective.", call. = FALSE)
  value
}

.ecm_evaluate <- function(X, Y, Theta, B, psi, C, lambda1, lambda2, psi.min) {
  n <- nrow(Y)
  k <- ncol(B)
  R <- Y - X %*% t(Theta)
  Sigma <- tcrossprod(B) + diag(psi, ncol(Y))
  root <- chol(Sigma)
  whitened <- forwardsolve(t(root), t(R))
  value <- n * sum(log(diag(root))) + sum(whitened^2) / 2 +
    .ecm_penalty(Theta, C, lambda1, lambda2)
  W <- chol2inv(root)
  D <- R %*% W
  H <- n * W - crossprod(D)
  coefficient.gradient <- -crossprod(X, D)
  loading.gradient <- H %*% B
  variance.gradient <- diag(H) / 2
  if (any(!is.finite(c(value, coefficient.gradient, loading.gradient, variance.gradient))))
    stop("Nonfinite ECM observed objective or gradient.", call. = FALSE)
  active <- psi.min > 0 & psi == psi.min
  projected <- variance.gradient
  projected[active] <- pmin(projected[active], 0)
  raw <- list(coefficient = .ecm_coefficient_score(coefficient.gradient, Theta, C, lambda1, lambda2),
              loadings = .weighted_norm(loading.gradient),
              variance = max(abs(projected)))
  stationarity <- lapply(raw, function(x) x / n)
  stationarity$maximum <- max(unlist(stationarity))
  stationarity$raw <- raw
  weighted.B <- sweep(B, 1L, psi, "/")
  V <- chol2inv(chol(diag(k) + crossprod(B, weighted.B)))
  M <- R %*% weighted.B %*% V
  S <- crossprod(M) + n * V
  if (any(!is.finite(c(M, V, S))))
    stop("Nonfinite ECM posterior moments.", call. = FALSE)
  list(objective = value, stationarity = stationarity, active = active,
       gradients = list(coefficient = t(coefficient.gradient),
                        loadings = loading.gradient, variance = variance.gradient,
                        variance.projected = projected),
       posterior = list(M = M, V = V, S = S))
}

# A declared floating comparison allowance, not a rigorous arithmetic error
# bound. The same rule applies to each frozen-Q stage and fresh observed L.
.ecm_increase <- function(before, after) {
  after - before > 128 * .Machine$double.eps * max(1, abs(before), abs(after))
}

#' Fit a separately named Gaussian ECM optimization reference
#'
#' A supplied-scale Gaussian optimization baseline with independent participants,
#' a standard normal latent factor, and diagonal positive residual variances.
#' This changes the estimator by omitting DrFARM inner debiasing. It is not a
#' replacement for DrFARM or a validated inferential procedure.
#'
#' @param X A finite n-by-p predictor matrix with positive finite squared column
#'   norms and at least two rows. No centering or standardization is performed.
#' @param Y A finite n-by-q response matrix, with q at least two.
#' @param Theta0 A finite q-by-p starting coefficient matrix. Entries with C=0
#'   must already be exactly zero; infeasible starts are rejected.
#' @param B0 A finite q-by-k starting loading matrix, with 1 <= k < q.
#'   Factor coordinates are unidentified up to orthogonal rotation.
#' @param psi0 Positive finite starting residual variances of length q, at or
#'   above any explicit positive psi.min bound.
#' @param lambda1,lambda2 Nonnegative entry and predictor-group penalties for
#'   the summed loss (no division by n).
#' @param C A q-by-p mask: 0 excludes, 1 includes an entry in both penalties,
#'   2 leaves it unpenalized and outside the group norm. NULL penalizes all.
#' @param control A named list: max.iter (500), objective.tol (1e-8),
#'   score.tol (1e-6), psi.min (0). psi.min is a nonnegative scalar or length-q
#'   vector. Zero means the open positive-variance domain, with no automatic
#'   floor. A positive bound explicitly changes the constrained target.
#' @param coefficient.control Named controls passed to \code{remMap.weighted}.
#'
#' @details With beta=t(Theta), R=Y-X beta and Sigma=B B'+diag(psi), the target
#'   is the observed negative Gaussian log likelihood plus the mask-aware
#'   sparse-group penalty:
#'   \deqn{L=\frac n2\log|\Sigma|+\frac12\mathrm{tr}(R\Sigma^{-1}R')+
#'   \lambda_1\sum_{r,j:C_{rj}=1}|\Theta_{rj}|+
#'   \lambda_2\sum_j\sqrt{\sum_{r:C_{rj}=1}\Theta_{rj}^2}.}
#'   The constant nq log(2 pi)/2 is omitted. Latent factors have mean zero,
#'   identity covariance and are independent of X and independent residuals;
#'   the Gaussian conditional and marginal mean coefficient is Theta. Only
#'   K=NULL is supported. No intercept, initialization, precision estimation,
#'   standardization, tuning selection or inference is performed internally.
#'
#'   Each cycle freezes M=E[Z|Y] and V=Var(Z_i|Y_i) at the accepted tuple and
#'   sets S=M'M+nV. It updates Theta with \code{remMap.weighted}, then computes
#'   B=(Y-X Theta')' M S inverse, then updates each variance from the full
#'   posterior expected residual sum of squares divided by n (subject to an
#'   explicit bound). The same updated Theta is used in both latter steps.
#'   No inner-debiased coefficient enters this cycle. All four Q stage values
#'   are evaluated under the one frozen posterior. The observed likelihood,
#'   posterior and gradients are then recomputed at the candidate tuple.
#'
#'   Accepted cycles must not increase any Q stage or observed L beyond
#'   128 times machine epsilon times max(1, absolute before, absolute after).
#'   This is a floating comparison allowance, not a proven error bound.
#'   Coefficient failure, numerical failure, open-domain variance boundary,
#'   Q increase or observed increase rejects the entire trial and returns the
#'   last accepted tuple with a warning. Invalid inputs or unrepresentable
#'   initial calculations raise an error. Computed zero expected RSS at a zero
#'   bound reports variance_boundary; it is not silently floored. A positive
#'   computed RSS whose division by n underflows to zero is a numerical failure.
#'   Floating underflow can also affect the RSS calculation itself.
#'
#'   Stationarity uses fresh observed gradients: maximum predictor-block
#'   sparse-group subgradient distance for Theta, Frobenius loading gradient
#'   norm, and maximum absolute variance gradient (projected at active positive
#'   bounds). Each is divided by n; these retain working-unit gradient scales,
#'   not dimensionless or parameter-error bounds. The maximum must be at most
#'   score.tol. An initially stationary tuple returns immediately. Otherwise
#'   convergence also requires relative observed-objective change at most
#'   objective.tol, using max(1, absolute previous L) as denominator.
#'   First-order stationarity is not global optimality; a zero loading start
#'   can be an absorbing stationary saddle. This nonconvex procedure has no
#'   inherited DrFARM estimation or inference validation claims.
#'
#' @return A list of class gaussian_ecm_reference containing Theta (q-by-p),
#'   B (q-by-k), diag.Psi, E.Z (fresh n-by-k posterior mean),
#'   posterior.covariance (k-by-k), and diagnostics. Diagnostics include method,
#'   inference.status (not_provided), converged, termination, iterations
#'   (accepted cycles), attempted, initial.objective, objective, objective.change
#'   (last accepted relative change; zero before any cycle), stationarity with
#'   raw residual magnitudes, gradients (coefficient q-by-p, loadings q-by-k,
#'   variance and projected variance), variance.bound.active, history (initial
#'   and accepted states), trials, failure and controls. Each trial records
#'   acceptance, reason, named Q stages, observed
#'   before/after values, and coefficient diagnostics; unavailable trial values
#'   are NA. No debiased estimate, standard error or p-value is returned.
#'
#' @examples
#' X <- matrix(c(-1, -1, 1, 1), 4, 1)
#' z <- c(1, -1, 1, -1)
#' u <- c(1, -1, -1, 1)
#' Y <- X %*% matrix(c(7/4, 11/4), 1, 2) +
#'   outer(z, rep(sqrt(15)/4, 2)) + outer(u, c(1, -1)/sqrt(2))
#' fit <- gaussian.ecm.reference(X, Y, matrix(c(1, 2), 2, 1),
#'   matrix(1, 2, 1), c(1, 1), lambda1=1, lambda2=0)
#' fit$diagnostics$termination # initial_stationary on this analytic fixture
#' @export
gaussian.ecm.reference <- function(X, Y, Theta0, B0, psi0, lambda1, lambda2,
                                   C = NULL, control = list(),
                                   coefficient.control = list()) {
  .validate_xy(X, Y, FALSE)
  q <- ncol(Y)
  if (q < 2L) stop("Y must have at least two response columns.", call. = FALSE)
  .validate_matrix(Theta0, "Theta0", q, ncol(X))
  .validate_matrix(B0, "B0", q)
  if (ncol(B0) < 1L || ncol(B0) >= q)
    stop("B0 must have at least one and fewer than q factor columns.", call. = FALSE)
  .validate_penalty(lambda1, lambda2, C, q, ncol(X))
  control <- .ecm_control(control, q)
  coefficient.control <- .weighted_control(coefficient.control)
  if (is.null(C)) C <- matrix(1L, q, ncol(X))
  if (any(Theta0[C == 0L] != 0))
    stop("Theta0 must be zero wherever C=0.", call. = FALSE)
  if (!is.numeric(psi0) || length(psi0) != q || any(!is.finite(psi0)) ||
      any(psi0 <= 0) || any(psi0 < control$psi.min))
    stop("psi0 must be positive finite variances at or above psi.min.", call. = FALSE)
  Theta <- Theta0
  B <- B0
  psi <- psi0
  current <- .ecm_evaluate(X, Y, Theta, B, psi, C, lambda1, lambda2, control$psi.min)
  initial <- current$objective
  change <- 0
  history <- data.frame(iteration = 0L, objective = initial, change = 0,
                        stationarity = current$stationarity$maximum)
  trials <- list()
  iterations <- 0L
  attempted <- 0L
  failure <- NULL
  termination <- if (current$stationarity$maximum <= control$score.tol)
    "initial_stationary" else "max_iter"
  while (termination == "max_iter" && attempted < control$max.iter) {
    attempted <- attempted + 1L
    posterior <- current$posterior
    trial <- list(accepted = FALSE, reason = "numerical_failure",
                  q = c(initial = NA_real_, coefficient = NA_real_,
                        loadings = NA_real_, variance = NA_real_),
                  observed.before = current$objective, observed.after = NA_real_,
                  inner = NULL)
    result <- tryCatch({
      trial$q["initial"] <- .ecm_Q(X, Y, Theta, B, psi, posterior, C, lambda1, lambda2)
      inner <- tryCatch(suppressWarnings(remMap.weighted(X,
        Y - tcrossprod(posterior$M, B), lambda1, lambda2, sigma = psi,
        C = C, Theta0 = Theta, control = coefficient.control)), error = identity)
      if (inherits(inner, "error")) {
        trial$reason <- "coefficient_failure"
        stop(conditionMessage(inner), call. = FALSE)
      }
      trial$inner <- inner$diagnostics
      if (!isTRUE(inner$diagnostics$converged)) {
        trial$reason <- "coefficient_failure"
        stop("Weighted coefficient CM did not converge.", call. = FALSE)
      }
      next.Theta <- inner$Theta0
      trial$q["coefficient"] <- .ecm_Q(X, Y, next.Theta, B, psi, posterior, C, lambda1, lambda2)
      R <- Y - X %*% t(next.Theta)
      S.root <- chol(posterior$S)
      next.B <- t(backsolve(S.root, forwardsolve(t(S.root), crossprod(posterior$M, R))))
      trial$q["loadings"] <- .ecm_Q(X, Y, next.Theta, next.B, psi, posterior, C, lambda1, lambda2)
      rss <- .ecm_expected_rss(R, next.B, posterior)
      raw.psi <- rss / nrow(Y)
      if (any(rss > 0 & raw.psi == 0))
        stop("Positive expected RSS divided by n underflowed to zero.", call. = FALSE)
      if (any(raw.psi == 0 & control$psi.min == 0)) {
        trial$reason <- "variance_boundary"
        stop("Zero expected RSS has no positive-variance CM minimizer; trial rejected.", call. = FALSE)
      }
      next.psi <- pmax(raw.psi, control$psi.min)
      trial$q["variance"] <- .ecm_Q(X, Y, next.Theta, next.B, next.psi, posterior, C, lambda1, lambda2)
      next.state <- .ecm_evaluate(X, Y, next.Theta, next.B, next.psi,
                                   C, lambda1, lambda2, control$psi.min)
      trial$observed.after <- next.state$objective
      trial$stationarity <- next.state$stationarity
      trial$variance.bound.active <- next.state$active
      if (any(vapply(seq_len(3L), function(j)
          .ecm_increase(trial$q[j], trial$q[j + 1L]), logical(1)))) {
        trial$reason <- "q_increase"
        stop("A frozen-Q CM stage increased; trial rejected.", call. = FALSE)
      }
      if (.ecm_increase(current$objective, next.state$objective)) {
        trial$reason <- "observed_increase"
        stop("Fresh observed objective increased; trial rejected.", call. = FALSE)
      }
      list(Theta = next.Theta, B = next.B, psi = next.psi, state = next.state)
    }, error = identity)
    if (inherits(result, "error")) {
      termination <- trial$reason
      failure <- conditionMessage(result)
      trial$failure <- failure
      trials[[attempted]] <- trial
      break
    }
    change <- abs(result$state$objective - current$objective) / max(1, abs(current$objective))
    Theta <- result$Theta
    B <- result$B
    psi <- result$psi
    current <- result$state
    iterations <- iterations + 1L
    trial$accepted <- TRUE
    trial$reason <- "accepted"
    trials[[attempted]] <- trial
    history <- rbind(history, data.frame(iteration = iterations,
      objective = current$objective, change = change, stationarity = current$stationarity$maximum))
    if (change <= control$objective.tol && current$stationarity$maximum <= control$score.tol)
      termination <- "stationarity_and_objective"
  }
  converged <- termination %in% c("initial_stationary", "stationarity_and_objective")
  if (!converged)
    warning("Gaussian ECM reference stopped with ", termination,
            "; returning the last accepted tuple.", call. = FALSE)
  structure(list(Theta = Theta, B = B, diag.Psi = psi,
    E.Z = current$posterior$M, posterior.covariance = current$posterior$V,
    diagnostics = list(method = "Gaussian ECM optimization reference",
      inference.status = "not_provided", converged = converged, termination = termination,
      iterations = iterations, attempted = attempted, initial.objective = initial,
      objective = current$objective, objective.change = change,
      stationarity = current$stationarity, gradients = current$gradients,
      variance.bound.active = current$active, history = history, trials = trials,
      failure = failure, control = control, coefficient.control = coefficient.control)),
    class = "gaussian_ecm_reference")
}
