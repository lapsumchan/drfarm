#!/usr/bin/env Rscript
# Reconcile preserved production paths with an explicit Gaussian target.
# Deterministic small diagnostic fixtures; no new production fitting mode.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(output = "outer-reconciliation", library = NULL)
  stopifnot(length(args) %% 2L == 0L)
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[i]); stopifnot(key %in% names(opt))
    opt[[key]] <- args[i + 1L]
  }
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library), .libPaths()))
  if (dir.exists(opt$output) && length(list.files(opt$output, all.files = TRUE, no.. = TRUE)))
    stop("Use a new output directory; prior receipts are preserved.")
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  source(file.path(dirname(normalizePath(script)), "gaussian-objective.R"), local = environment())
  library(drfarm)
  log <- file(file.path(opt$output, "console.log"), "wt")
  sink(log, split = TRUE); on.exit({sink(); base::close(log)}, add = TRUE)
  start <- proc.time()
  count <- 0L
  check <- function(ok, label) {
    if (!isTRUE(ok)) stop("FAIL: ", label)
    count <<- count + 1L; cat("PASS:", label, "\n")
  }
  close <- function(a, b, label, tol = 1e-10) {
    check(identical(dim(a), dim(b)) && length(a) == length(b) &&
            all(is.finite(c(a, b))) && max(abs(a - b)) <= tol, label)
  }
  fixtures <- list(
    zero_factor = list(X = matrix(c(1, -1), 2, 1),
      Y = rbind(c(3, 4), c(-1, -2)), Theta = matrix(0, 2, 1),
      B = matrix(0, 2, 1), psi = c(5, 10), precision = matrix(1),
      lambda1 = 1, lambda2 = 0, max.iter = 5L),
    stationary_nonzero = local({
      X <- matrix(c(-1, -1, 1, 1), 4, 1)
      z <- c(1, -1, 1, -1); u <- c(1, -1, -1, 1)
      list(X = X, Y = X %*% matrix(c(7/4, 11/4), 1, 2) +
        outer(z, rep(sqrt(15)/4, 2)) + outer(u, c(1/sqrt(2), -1/sqrt(2))),
        Theta = matrix(c(1, 2), 2, 1), B = matrix(1, 2, 1), psi = c(1, 1),
        precision = matrix(1), lambda1 = 1, lambda2 = 0, max.iter = 1L)
    }),
    monitor_disagreement = list(X = matrix(c(-1, -1, 1, 1), 4, 1),
      Y = rbind(c(-.2, .8), c(-2.8, -1.8), c(2.2, 3.2), c(.8, -2.2)),
      Theta = matrix(0, 2, 1), B = matrix(c(1, .5), 2, 1), psi = c(1, 2),
      precision = matrix(1), lambda1 = 1, lambda2 = 0, max.iter = 1L))
  saveRDS(fixtures, file.path(opt$output, "fixtures.rds"))
  results <- list(); rows <- list()
  for (name in names(fixtures)) for (mode in c("historical", "weighted")) {
    a <- fixtures[[name]]
    run <- paste(name, mode, sep = "-")
    traced <- do.call(trace_outer_fixed, c(a, list(mode = mode)))
    check(traced$instrumentation_parity, paste(run, "observer/uninstrumented exact parity"))
    posterior <- NULL; anchor <- NULL
    for (r in traced$records) {
      if (r$stage == "posterior") {
        posterior <- gaussian_posterior(r$X, r$Y, r$Theta, r$B, r$psi, r$d)
        close(r$mean, posterior$mean, paste(run, r$iteration, "posterior mean"))
        close(r$second, posterior$second, paste(run, r$iteration, "posterior second moment"))
        anchor <- r
      }
      T <- if (r$stage == "debias") r$Theta.db else r$Theta
      L <- gaussian_observed(r$X, r$Y, T, r$B, r$psi, a$lambda1, a$lambda2)
      Q <- gaussian_Q(r$X, r$Y, T, r$B, r$psi, posterior, a$lambda1, a$lambda2)
      H <- historical_monitor(r$X, r$Y, T, r$psi, a$lambda1, a$lambda2)
      rows[[length(rows) + 1L]] <- data.frame(fixture = name, mode = mode,
        iteration = r$iteration, stage = r$stage, observed_penalized = L,
        frozen_Q_penalized = Q, historical_monitor = H,
        theta1 = T[1, 1], theta2 = T[2, 1], B1 = r$B[1, 1], B2 = r$B[2, 1],
        psi1 = r$psi[1], psi2 = r$psi[2])
      if (r$stage == "monitor") close(r$monitor, H, paste(run, "actual monitor equality"))
      if (r$stage == "variance") {
        full_db <- gaussian_expected_rss(r$X, r$Y, r$Theta.db, r$B, posterior) / nrow(r$Y)
        close(r$psi, full_db, paste(run, "compressed variance equals full posterior RSS at beta_db"))
        oldL <- gaussian_observed(anchor$X, anchor$Y, anchor$Theta, anchor$B,
          anchor$psi, a$lambda1, a$lambda2)
        oldQ <- gaussian_Q(anchor$X, anchor$Y, anchor$Theta, anchor$B,
          anchor$psi, posterior, a$lambda1, a$lambda2)
        next_post <- gaussian_posterior(r$X, r$Y, r$Theta, r$B, r$psi, r$d)
        kl <- gaussian_posterior_kl(posterior, next_post)
        close(L - oldL, Q - oldQ - kl, paste(run, "Gaussian Q-minus-posterior-KL identity"))
      }
    }
    results[[run]] <- traced
  }
  tab <- do.call(rbind, rows)
  write.table(tab, file.path(opt$output, "stages.tsv"), sep = "\t", row.names = FALSE)
  saveRDS(results, file.path(opt$output, "traces.rds"))
  # The closed trajectory is derived independently before execution.
  zero <- results[["zero_factor-weighted"]]$fit
  close(zero$Theta, matrix(c(1.5, 2.5), 2, 1), "exact converged sparse beta")
  close(zero$diag.Psi, c(1, 1), "exact converged psi")
  check(zero$diagnostics$converged && zero$diagnostics$iterations == 3 &&
        zero$diagnostics$termination == "loss_tolerance", "declared outer convergence after three attempts")
  a <- fixtures$zero_factor
  rss <- colSums((a$Y - a$X %*% t(zero$Theta))^2)
  score_psi <- nrow(a$Y)/(2 * zero$diag.Psi) - rss/(2 * zero$diag.Psi^2)
  close(score_psi, c(-.25, -.25), "nonzero Gaussian variance score at declared convergence")
  L <- local({ a0 <- a; fit0 <- zero
    function(psi) gaussian_observed(a0$X, a0$Y, fit0$Theta, as.matrix(fit0$B), psi, 1, 0)
  })
  finite_difference <- vapply(1:2, function(j) {
    up <- down <- zero$diag.Psi; up[j] <- up[j] + 1e-6; down[j] <- down[j] - 1e-6
    (L(up) - L(down))/2e-6
  }, numeric(1))
  close(finite_difference, score_psi, "independent finite-difference variance gradient", 1e-8)
  close(L(c(1.25, 1.25)), 6 + 2*log(1.25), "lower objective from matching sparse-residual variance")
  check(L(c(1.25, 1.25)) < L(c(1, 1)), "declared fixed point is not Gaussian stationary")
  zrows <- subset(tab, fixture == "zero_factor" & mode == "weighted" & stage == "monitor")
  close(zrows$observed_penalized, c(15, 6.5, 6.5), "exact three-step objective trajectory")

  # Nonzero factor fixture begins at Gaussian penalized stationarity.
  a <- fixtures$stationary_nonzero
  residual <- a$Y - a$X %*% t(a$Theta)
  Sigma <- tcrossprod(a$B) + diag(a$psi)
  close(crossprod(residual)/4, Sigma, "nonzero-factor empirical covariance equals model covariance")
  score <- -crossprod(a$X, residual) %*% solve(Sigma) + matrix(1, 1, 2)
  close(score, matrix(0, 1, 2), "initial Gaussian penalized coefficient score is zero")
  post <- gaussian_posterior(a$X, a$Y, a$Theta, a$B, a$psi)
  cm <- gaussian_factor_cm(a$X, a$Y, a$Theta, post)
  close(cm$B, a$B, "coherent factor step stays at nonzero stationary state")
  close(cm$psi, a$psi, "coherent variance step stays at nonzero stationary state")
  hybrid <- results[["stationary_nonzero-weighted"]]$fit
  close(as.matrix(hybrid$B), matrix(7/8, 2, 1), "hybrid loadings change to exact 7/8")
  close(hybrid$diag.Psi, rep(59/64, 2), "hybrid variances change to exact 59/64")
  full_sparse <- gaussian_expected_rss(a$X, a$Y, hybrid$Theta, as.matrix(hybrid$B), post)/4
  close(full_sparse, rep(65/64, 2), "returned sparse beta needs posterior variance 65/64")
  before <- gaussian_observed(a$X, a$Y, a$Theta, a$B, a$psi, 1, 0)
  after <- gaussian_observed(a$X, a$Y, hybrid$Theta, as.matrix(hybrid$B), hybrid$diag.Psi, 1, 0)
  check(after > before + .05, "complete hybrid leaves Gaussian stationary state with larger objective")

  # A legitimate frozen-Q CM comparison is a reference calculation only.
  a <- fixtures$monitor_disagreement
  post <- gaussian_posterior(a$X, a$Y, a$Theta, a$B, a$psi)
  coeff <- remMap.weighted(a$X, a$Y - post$mean %*% t(a$B), 1, 0,
                           sigma = a$psi, control = list(tol = 1e-11))$Theta0
  cm <- gaussian_factor_cm(a$X, a$Y, coeff, post)
  comparisons <- rbind(
    before = c(Q = gaussian_Q(a$X, a$Y, coeff, a$B, a$psi, post, 1, 0),
      L = gaussian_observed(a$X, a$Y, coeff, a$B, a$psi, 1, 0),
      H = historical_monitor(a$X, a$Y, coeff, a$psi, 1, 0)),
    coherent_cm = c(Q = gaussian_Q(a$X, a$Y, coeff, cm$B, cm$psi, post, 1, 0),
      L = gaussian_observed(a$X, a$Y, coeff, cm$B, cm$psi, 1, 0),
      H = historical_monitor(a$X, a$Y, coeff, cm$psi, 1, 0)))
  check(comparisons[2,"Q"] < comparisons[1,"Q"] &&
        comparisons[2,"L"] < comparisons[1,"L"] && comparisons[2,"H"] > comparisons[1,"H"],
        "historical monitor rises during legitimate Gaussian/Q descent")
  write.table(comparisons, file.path(opt$output, "monitor-disagreement.tsv"), sep = "\t", col.names = NA)

  # Basis/orientation check of objective only; this does not repair K fitting.
  K <- .5 ^ abs(outer(1:4, 1:4, `-`))
  eigenK <- eigen(K); U <- eigenK$vectors
  E <- a$Y - a$X %*% t(a$Theta)
  dense <- kronecker(tcrossprod(a$B), K) + kronecker(diag(a$psi), diag(4))
  denseL <- gaussian_logdet(dense)/2 + drop(crossprod(c(E), solve(dense, c(E))))/2
  rotateL <- gaussian_observed(t(U) %*% a$X, t(U) %*% a$Y, a$Theta,
                               a$B, a$psi, d = eigenK$values)
  close(rotateL, denseL, "K-eigenbasis objective equals dense original-basis Gaussian", 1e-9)
  # C2 penalty mismatch and internal C0 drift are distinguished from feasibility.
  mask <- matrix(c(0, 2), 2, 1)
  masked <- remMap.weighted(a$X, a$Y, 1, 1, sigma = a$psi, C = mask)$Theta0
  db <- masked + t(a$precision %*% crossprod(a$X, a$Y - a$X %*% t(masked)))/4
  check(masked[1,1] == 0 && abs(db[1,1]) > .1, "internal debiasing can leave a C0 constraint")
  check(is.infinite(gaussian_penalty(db, 1, 1, mask)), "constrained Gaussian target rejects off-mask debiased point")
  close(gaussian_penalty(masked, 1, 1, mask), 0, "C2 excluded from Gaussian coefficient penalty")
  check(sum(abs(masked)) + sqrt(sum(masked^2)) > 0, "historical monitor penalizes nonzero C2 entry")

  receipt <- list(status = "PASS", checks = count, command = commandArgs(),
    package = as.character(packageVersion("drfarm")), library = find.package("drfarm"),
    session = sessionInfo(), threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")),
    dimensions = lapply(fixtures, function(a) c(n = nrow(a$X), p = ncol(a$X), q = ncol(a$Y), k = ncol(a$B))),
    random_generation = "none; deterministic explicit fixtures", elapsed = proc.time() - start,
    scope = "Fixed initializer injected in local actual-production closure. Numeric steps unchanged; observer parity checked. No new production algorithm or inference.",
    findings = list(zero_point_variance_score = score_psi,
      zero_point_objective = L(c(1,1)), zero_point_improved_variance_objective = L(c(1.25,1.25)),
      nonzero_start_objective = before, nonzero_hybrid_objective = after,
      coherent_monitor_comparison = comparisons),
    warnings = lapply(results, `[[`, "warnings"))
  saveRDS(receipt, file.path(opt$output, "receipt.rds"))
  writeLines(capture.output(str(receipt, max.level = 4)), file.path(opt$output, "receipt.txt"))
  writeLines(capture.output(sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
  print(comparisons); cat("PASS:", count, "checks. Full outer Gaussian failures retained as findings.\n")
}
main()
