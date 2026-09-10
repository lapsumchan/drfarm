# Stable Euclidean norm for finite vectors; finite-computation failures are
# reported rather than interpreted as satisfied stopping criteria.
.weighted_norm <- function(x) {
  if (!length(x)) return(0)
  scale <- max(abs(x))
  if (!is.finite(scale)) stop("Nonfinite weighted-solver norm.", call. = FALSE)
  if (scale == 0) return(0)
  value <- scale * sqrt(sum((x / scale)^2))
  if (!is.finite(value)) stop("Weighted-solver norm overflow.", call. = FALSE)
  value
}

.weighted_control <- function(control) {
  defaults <- list(tol = 1e-8, max.sweeps = 1000L, root.tol = 1e-14,
                   root.maxit = 200L)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% names(defaults)) ||
       anyDuplicated(names(control)))))
    stop("control must be a named list of tol, max.sweeps, root.tol, root.maxit.", call. = FALSE)
  defaults[names(control)] <- control
  for (key in names(defaults)) {
    value <- defaults[[key]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(key, " must be a positive finite scalar.", call. = FALSE)
  }
  for (key in c("max.sweeps", "root.maxit")) {
    if (defaults[[key]] != floor(defaults[[key]]) || defaults[[key]] > .Machine$integer.max)
      stop(key, " must be a positive integer within the R integer range.", call. = FALSE)
  }
  if (defaults$tol >= 1 || defaults$root.tol >= 1)
    stop("tol and root.tol must be smaller than one.", call. = FALSE)
  defaults
}

# One penalized predictor block. For h>0 and u=soft(s,lambda1), a nonzero
# optimum is u/(h+t); t>0 solves ||t*u/(h+t)||=lambda2. Normalize the scalar
# root onto [0,1] via alpha=t/(H+t), H=max(h), to avoid an unbounded bracket.
.weighted_block <- function(h, score, lambda1, lambda2, control) {
  u <- sign(score) * pmax(abs(score) - lambda1, 0)
  magnitude <- .weighted_norm(u)
  if (magnitude <= lambda2)
    return(list(beta = rep(0, length(h)), iterations = 0L, residual = 0))
  if (lambda2 == 0)
    return(list(beta = u / h, iterations = 0L, residual = 0))
  H <- max(h)
  hn <- h / H
  if (any(hn <= 0) || any(!is.finite(hn)))
    stop("Curvature ratios exceed the supported floating-point range.", call. = FALSE)
  target <- lambda2 / magnitude
  f <- function(alpha) {
    fraction <- alpha / (hn * (1 - alpha) + alpha)
    .weighted_norm((u / magnitude) * fraction) - target
  }
  root <- stats::uniroot(f, c(0, 1), f.lower = -target, f.upper = 1 - target,
                         tol = control$root.tol, maxiter = control$root.maxit,
                         check.conv = TRUE)
  alpha <- root$root
  beta <- (u / H) * (1 - alpha) / (hn * (1 - alpha) + alpha)
  list(beta = beta, iterations = root$iter, residual = abs(f(alpha)))
}

.weighted_objective <- function(E, beta, sigma, mask, lambda1, lambda2) {
  value <- sum(sweep(E, 2L, sqrt(sigma), "/")^2) / 2 +
    lambda1 * sum(abs(beta[mask == 1]))
  for (j in seq_len(nrow(beta)))
    value <- value + lambda2 * .weighted_norm(beta[j, mask[j, ] == 1])
  if (!is.finite(value)) stop("Nonfinite weighted objective.", call. = FALSE)
  value
}

# Exact subgradient-distance formula for this finite convex objective.
# The raw residual has objective-gradient units. Scaling does not bound error
# in coefficients or inference, particularly for ill-conditioned designs.
.weighted_kkt <- function(X, E, beta, sigma, mask, lambda1, lambda2, scales) {
  gradient <- -sweep(crossprod(X, E), 2L, sigma, "/")
  if (any(!is.finite(gradient))) stop("Nonfinite weighted gradient.", call. = FALSE)
  residuals <- numeric(nrow(beta))
  for (j in seq_len(nrow(beta))) {
    penalized <- mask[j, ] == 1
    b <- beta[j, penalized]
    g <- gradient[j, penalized]
    radius <- .weighted_norm(b)
    if (radius == 0) {
      pen <- max(0, .weighted_norm(sign(g) * pmax(abs(g) - lambda1, 0)) - lambda2)
    } else {
      v <- g + lambda2 * (b / radius)
      v <- ifelse(b != 0, v + lambda1 * sign(b),
                  sign(v) * pmax(abs(v) - lambda1, 0))
      pen <- .weighted_norm(v)
    }
    residuals[j] <- .weighted_norm(c(pen, gradient[j, mask[j, ] == 2]))
  }
  list(raw = max(residuals), scaled = max(residuals / scales), by.block = residuals)
}

#' Fit a fixed-variance weighted sparse-group regression
#'
#' Solves a convex coefficient subproblem on the supplied working scale.
#' Response variances remain fixed throughout this fit. This is a separately
#' named alternative to the historical remMap coefficient update.
#'
#' @param X A finite n-by-p predictor matrix, with positive finite squared
#'   column norms and at least two rows. Rank deficiency is allowed, but
#'   coefficients need not then be uniquely identified by the objective.
#' @param Y A finite n-by-q response or factor-adjusted response matrix.
#' @param lambda1,lambda2 Nonnegative entry and predictor-group penalties,
#'   respectively. The residual loss is summed, not divided by n.
#' @param sigma Positive finite response variances, of length q. These are
#'   variances in the supplied Y units, not standard deviations.
#' @param C A q-by-p mask: 0 fixes a coefficient to zero; 1 includes it in
#'   both penalties; 2 leaves it unpenalized and outside the group norm.
#'   NULL penalizes all entries.
#' @param Theta0 An optional finite q-by-p starting coefficient matrix on the
#'   supplied working scale. Entries with C=0 are set to zero before fitting.
#' @param control A named list: tol (default 1e-8), max.sweeps (1000),
#'   root.tol (1e-14) and root.maxit (200). tol bounds the scaled KKT residual;
#'   root.tol is an absolute tolerance on a scalar transformed root in [0,1].
#'
#' @details No centering, standardization, intercept fitting, variance estimation
#'   or inference is performed. With beta=t(Theta), the objective is
#'   \deqn{\frac12\sum_{i,r}(Y_{ir}-(X\beta)_{ir})^2/\sigma_r
#'   +\lambda_1\sum_{j,r:C_{rj}=1}|\beta_{jr}|
#'   +\lambda_2\sum_j\sqrt{\sum_{r:C_{rj}=1}\beta_{jr}^2}.}
#'   Cyclic predictor-block descent uses a one-dimensional root for unequal
#'   response curvature. Zero-group and unpenalized solutions are explicit.
#'   The residual matrix is recomputed from X, Y and beta each sweep before
#'   checking the full KKT equations. The maximum per-block subgradient-distance
#'   residual is scaled by the maximum of 1, the norm of that block's weighted
#'   X'Y score, lambda1 times the square root of its penalized-entry count, and
#'   lambda2 (penalty scales are omitted for entirely unpenalized blocks).
#'   A small residual is a finite stopping check, not a coefficient-error bound.
#'   Exact zero subgradient is sufficient for global minimization of this convex
#'   subproblem; it establishes no global optimum for the outer DrFARM model.
#'   On an exhausted budget or numerical solve failure the last finite iterate
#'   is returned with a warning and converged=FALSE. Inspect termination and
#'   failure. If a trial cannot be evaluated, the previous fully evaluated sweep
#'   is returned. Invalid inputs or unrepresentable initialization calculations raise
#'   an error before iteration. Objective history contains evaluated iterates, not a proof of
#'   exact arithmetic or a bound on statistical error.
#'
#' @return A list with Theta0 (q-by-p) and diagnostics. Diagnostics include
#'   converged, termination, sweeps (attempted, including a partial or rolled-back
#'   final sweep), objective, initial.objective, kkt.residual
#'   (raw), kkt.scaled, threshold, max.sweeps, root.iterations,
#'   max.block.root.residual (normalized scalar-equation residual),
#'   objective.trace, failure and control. No p-values are produced.
#'
#' @examples
#' X <- matrix(c(1, -1) / sqrt(2), 2, 1)
#' Y <- X %*% matrix(c(8, 5.5), 1, 2)
#' fit <- remMap.weighted(X, Y, lambda1=2, lambda2=5, sigma=c(1, 0.25))
#' fit$Theta0 # analytic optimum: c(3, 4)
#' fit$diagnostics$kkt.residual
#' @export
remMap.weighted <- function(X, Y, lambda1, lambda2,
                            sigma = rep(1, ncol(Y)), C = NULL, Theta0 = NULL,
                            control = list()) {
  .validate_xy(X, Y, FALSE)
  .validate_penalty(lambda1, lambda2, C, ncol(Y), ncol(X))
  control <- .weighted_control(control)
  if (!is.numeric(sigma) || length(sigma) != ncol(Y) ||
      any(!is.finite(sigma)) || any(sigma <= 0))
    stop("sigma must contain one positive finite variance per response.", call. = FALSE)
  p <- ncol(X)
  q <- ncol(Y)
  mask <- if (is.null(C)) matrix(1L, p, q) else t(C)
  beta <- matrix(0, p, q)
  if (!is.null(Theta0)) {
    .validate_matrix(Theta0, "Theta0", q, p)
    beta <- t(Theta0)
  }
  beta[mask == 0] <- 0
  dimnames(beta) <- list(colnames(X), colnames(Y))
  d <- colSums(X^2)
  curvature <- outer(d, sigma, "/")
  if (any(!is.finite(curvature)) || any(curvature <= 0))
    stop("Weighted curvature is outside the supported floating-point range.", call. = FALSE)
  zero.score <- sweep(crossprod(X, Y), 2L, sigma, "/")
  scales <- vapply(seq_len(p), function(j) {
    npen <- sum(mask[j, ] == 1)
    max(1, .weighted_norm(zero.score[j, mask[j, ] != 0]),
        if (npen) lambda1 * sqrt(npen) else 0, if (npen) lambda2 else 0)
  }, numeric(1))
  if (any(!is.finite(scales))) stop("Nonfinite KKT scaling.", call. = FALSE)
  E <- Y - X %*% beta
  value <- .weighted_objective(E, beta, sigma, mask, lambda1, lambda2)
  initial <- value
  history <- value
  kkt <- .weighted_kkt(X, E, beta, sigma, mask, lambda1, lambda2, scales)
  sweeps <- 0L
  root.iterations <- 0
  root.residual <- 0
  termination <- "max_sweeps"
  failure <- NULL
  if (kkt$scaled <= control$tol) termination <- "kkt"
  while (termination == "max_sweeps" && sweeps < control$max.sweeps) {
    sweeps <- sweeps + 1L
    before <- value
    last.valid <- list(beta = beta, E = E, value = value, kkt = kkt)
    for (j in seq_len(p)) {
      old <- beta[j, ]
      partial <- E + tcrossprod(X[, j], old)
      score <- drop(crossprod(X[, j], partial)) / sigma
      if (any(!is.finite(score))) {
        termination <- "numerical_failure"
        failure <- "Nonfinite partial-residual score."
        break
      }
      next.beta <- rep(0, q)
      free <- mask[j, ] == 2
      next.beta[free] <- score[free] / curvature[j, free]
      pen <- mask[j, ] == 1
      if (any(pen)) {
        block <- tryCatch(.weighted_block(curvature[j, pen], score[pen],
                           lambda1, lambda2, control), error = identity)
        if (inherits(block, "error")) {
          termination <- "root_failure"
          failure <- conditionMessage(block)
          break
        }
        next.beta[pen] <- block$beta
        root.iterations <- root.iterations + block$iterations
        root.residual <- max(root.residual, block$residual)
      }
      if (any(!is.finite(next.beta))) {
        termination <- "numerical_failure"
        failure <- "Nonfinite coefficient update."
        break
      }
      beta[j, ] <- next.beta
      E <- partial - tcrossprod(X[, j], next.beta)
    }
    # Refresh residuals for the reported objective and *global* KKT check.
    E <- Y - X %*% beta
    evaluated <- tryCatch(list(
      objective = .weighted_objective(E, beta, sigma, mask, lambda1, lambda2),
      kkt = .weighted_kkt(X, E, beta, sigma, mask, lambda1, lambda2, scales)),
      error = identity)
    if (inherits(evaluated, "error")) {
      beta <- last.valid$beta
      E <- last.valid$E
      value <- last.valid$value
      kkt <- last.valid$kkt
      termination <- "numerical_failure"
      failure <- paste("Could not evaluate trial iterate:", conditionMessage(evaluated))
      break
    }
    value <- evaluated$objective
    kkt <- evaluated$kkt
    history <- c(history, value)
    if (termination != "max_sweeps") break
    if (value > before + 128 * .Machine$double.eps * max(1, abs(before))) {
      termination <- "objective_increase"
      failure <- "Weighted objective increased beyond the roundoff allowance."
      break
    }
    if (kkt$scaled <= control$tol) termination <- "kkt"
  }
  converged <- identical(termination, "kkt")
  diagnostics <- list(converged = converged, termination = termination,
    sweeps = sweeps, objective = value, initial.objective = initial,
    kkt.residual = kkt$raw, kkt.scaled = kkt$scaled,
    threshold = control$tol, max.sweeps = control$max.sweeps,
    root.iterations = root.iterations, max.block.root.residual = root.residual,
    objective.trace = history, failure = failure, control = control,
    criterion = "scaled subgradient distance for the fixed-variance coefficient objective")
  if (!converged) warning("Weighted coefficient solver stopped: ", termination,
                         if (!is.null(failure)) paste0("; ", failure), call. = FALSE)
  list(Theta0 = t(beta), diagnostics = diagnostics)
}
