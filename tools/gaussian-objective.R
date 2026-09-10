# Diagnostic-only Gaussian references; not a production DrFARM fitting API.
# Work on supplied (possibly K-rotated) data, with latent variances d_i fixed.
# All NLLs omit nq/2*log(2*pi); Q also omits the parameter-constant latent prior.
gaussian_penalty <- function(Theta, lambda1 = 0, lambda2 = 0, C = NULL) {
  if (is.null(C)) C <- matrix(1L, nrow(Theta), ncol(Theta))
  stopifnot(identical(dim(C), dim(Theta)), all(C %in% 0:2),
            lambda1 >= 0, lambda2 >= 0)
  if (any(Theta[C == 0] != 0)) return(Inf)
  lambda1 * sum(abs(Theta[C == 1])) + lambda2 * sum(vapply(
    seq_len(ncol(Theta)), function(j) sqrt(sum(Theta[C[, j] == 1, j]^2)), numeric(1)))
}
gaussian_validate <- function(X, Y, Theta, B, psi, d) {
  stopifnot(is.matrix(X), is.matrix(Y), is.matrix(Theta), is.matrix(B),
    nrow(X) == nrow(Y), nrow(Theta) == ncol(Y), ncol(Theta) == ncol(X),
    nrow(B) == ncol(Y), length(psi) == ncol(Y), length(d) == nrow(Y),
    all(is.finite(c(X, Y, Theta, B, psi, d))), all(psi > 0), all(d > 0))
}
gaussian_logdet <- function(A) 2 * sum(log(diag(chol(A))))
gaussian_observed <- function(X, Y, Theta, B, psi, lambda1 = 0, lambda2 = 0,
                              C = NULL, d = rep(1, nrow(X))) {
  gaussian_validate(X, Y, Theta, B, psi, d)
  E <- Y - X %*% t(Theta)
  value <- 0
  for (i in seq_len(nrow(Y))) {
    R <- chol(diag(psi, ncol(Y)) + d[i] * tcrossprod(B))
    z <- forwardsolve(t(R), E[i, ])
    value <- value + sum(log(diag(R))) + sum(z^2) / 2
  }
  value + gaussian_penalty(Theta, lambda1, lambda2, C)
}
gaussian_posterior <- function(X, Y, Theta, B, psi, d = rep(1, nrow(X))) {
  gaussian_validate(X, Y, Theta, B, psi, d)
  E <- Y - X %*% t(Theta)
  k <- ncol(B)
  means <- matrix(0, nrow(Y), k)
  covariance <- vector("list", nrow(Y))
  for (i in seq_len(nrow(Y))) {
    V <- solve(diag(1 / d[i], k) + crossprod(B, sweep(B, 1L, psi, "/")))
    covariance[[i]] <- V
    means[i, ] <- V %*% crossprod(B, E[i, ] / psi)
  }
  list(mean = means, covariance = covariance,
       second = crossprod(means) + Reduce(`+`, covariance))
}
gaussian_expected_rss <- function(X, Y, Theta, B, posterior) {
  residual <- Y - X %*% t(Theta) - posterior$mean %*% t(B)
  value <- colSums(residual^2)
  for (V in posterior$covariance) value <- value + diag(B %*% V %*% t(B))
  value
}
gaussian_Q <- function(X, Y, Theta, B, psi, posterior,
                       lambda1 = 0, lambda2 = 0, C = NULL) {
  stopifnot(all(is.finite(psi)), all(psi > 0))
  nrow(Y) * sum(log(psi)) / 2 +
    sum(gaussian_expected_rss(X, Y, Theta, B, posterior) / psi) / 2 +
    gaussian_penalty(Theta, lambda1, lambda2, C)
}
# Gaussian KL(old posterior || new posterior), independently useful for
# L(new)-L(old) = Q(new|old)-Q(old|old) - KL(old || new).
gaussian_posterior_kl <- function(old, new) {
  k <- ncol(old$mean)
  sum(vapply(seq_len(nrow(old$mean)), function(i) {
    A <- old$covariance[[i]]; D <- new$covariance[[i]]
    delta <- new$mean[i, ] - old$mean[i, ]
    (sum(diag(solve(D, A))) + drop(crossprod(delta, solve(D, delta))) - k +
       gaussian_logdet(D) - gaussian_logdet(A)) / 2
  }, numeric(1)))
}
# Uses the full expected RSS (squared posterior-mean residual PLUS posterior
# covariance), rather than copying the compressed production variance formula.
gaussian_factor_cm <- function(X, Y, Theta, posterior) {
  residual <- Y - X %*% t(Theta)
  B <- t(solve(posterior$second, crossprod(posterior$mean, residual)))
  psi <- gaussian_expected_rss(X, Y, Theta, B, posterior) / nrow(Y)
  stopifnot(all(is.finite(psi)), all(psi > 0))
  list(B = B, psi = psi)
}
historical_monitor <- function(X, Y, Theta, psi, lambda1 = 0, lambda2 = 0) {
  # Precisely the current all-entry penalty and diagonal-only residual score.
  nrow(Y) * sum(log(psi)) / 2 +
    sum(sweep((Y - X %*% t(Theta))^2, 2L, psi, "/")) / 2 +
    lambda1 * sum(abs(Theta)) + lambda2 * sum(sqrt(colSums(Theta^2)))
}

# Instruments a local copy of the ACTUAL installed DrFARM.one closure. Replaces
# only its initial fa result with a declared fixed state; adds observers after
# uniquely identified assignments inside the existing while body. No production
# numerical expression or package namespace is changed. This is not a new fit.
trace_outer_fixed <- function(X, Y, Theta, B, psi, precision, lambda1, lambda2,
                              mode, C = NULL, max.iter = 5L, thres = 1e-10) {
  original <- drfarm::DrFARM.one
  capture <- new.env(parent = environment(original))
  capture$fa <- function(...) list(loadings = B, uniquenesses = psi)
  records <- list()
  capture$.outer_observe <- function(stage, state) {
    take <- function(name) if (exists(name, envir = state, inherits = FALSE))
      get(name, envir = state, inherits = FALSE) else NULL
    records[[length(records) + 1L]] <<- list(
      stage = stage, iteration = take("attempted"), X = take("X"), Y = take("Y"),
      d = take("d"), Theta = t(take("Theta.t")), B = as.matrix(take("B")),
      psi = take("diag.Psi"), Theta.db = if (!is.null(take("Theta.db.t")))
        t(take("Theta.db.t")) else NULL,
      mean = if (!is.null(take("E.zt"))) t(take("E.zt")) else NULL,
      second = take("E.zzt"), posterior.seconds = take("E.zzt.List"),
      monitor = if (stage == "monitor") take("loss") else NULL)
    invisible(NULL)
  }
  f <- original
  environment(f) <- capture
  pieces <- as.list(body(f))
  is_while <- vapply(pieces, function(z) is.call(z) && identical(z[[1]], as.name("while")), logical(1))
  stopifnot(sum(is_while) == 1L)
  wi <- which(is_while)
  loop <- pieces[[wi]]
  statements <- as.list(loop[[3]])
  markers <- c("E.zzt.inv" = "posterior", "Theta.t" = "coefficient",
               "Theta.db.t" = "debias", "diag.Psi" = "variance", "loss" = "monitor")
  seen <- setNames(integer(length(markers)), names(markers))
  changed <- list(statements[[1]])
  for (z in statements[-1]) {
    changed[[length(changed) + 1L]] <- z
    if (is.call(z) && identical(z[[1]], as.name("<-")) && is.symbol(z[[2]])) {
      name <- as.character(z[[2]])
      if (name %in% names(markers)) {
        # Theta.t also has a restore at loop bottom. Observe only native result.
        if (name == "Theta.t" && !identical(z[[3]], quote(inner.fit$phi))) next
        if (name == "diag.Psi" && identical(z[[3]], quote(prev.Psi))) next
        seen[name] <- seen[name] + 1L
        changed[[length(changed) + 1L]] <- substitute(
          .outer_observe(STAGE, environment()), list(STAGE = unname(markers[name])))
      }
    }
  }
  stopifnot(all(seen == 1L))
  loop[[3]] <- as.call(changed)
  pieces[[wi]] <- loop
  body(f) <- as.call(pieces)
  args <- list(X = X, Y = Y, Theta0 = Theta, precM = precision, k = ncol(B),
               lambda1 = lambda1, lambda2 = lambda2, C = C,
               standardize = FALSE, max.iter = max.iter, thres = thres,
               coefficient.update = mode)
  if (mode == "weighted") args$weighted.control <- list(tol = 1e-11)
  warnings <- character()
  fit <- withCallingHandlers(do.call(f, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  # Matched uninstrumented fixed-initializer closure: observation must not alter
  # the returned numerical results, diagnostics or warnings.
  plain <- original; environment(plain) <- capture
  plain_warnings <- character()
  reference <- withCallingHandlers(do.call(plain, args), warning = function(w) {
    plain_warnings <<- c(plain_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  stopifnot(identical(fit, reference), identical(warnings, plain_warnings))
  list(fit = fit, records = records, warnings = warnings,
       instrumentation_parity = TRUE, source_body = body(original),
       scope = "Fixed declared fa initializer in a local closure; original update expressions unchanged")
}
