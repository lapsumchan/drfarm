#!/usr/bin/env Rscript
# Start-specific optimization comparison, not validation of DrFARM inference.
# Reuses explicit fixtures and preparation from the prior two executed slices.
# Run from the repository root; the --fixtures/--inputs/--preparation arguments
# bind the exact saved assets. Each execution requires a fresh output directory.
main <- function() {
  opt <- list(output = "gaussian-ecm-comparison", library = NULL,
              fixtures = NULL, inputs = NULL, preparation = NULL,
              `max-iter` = "500", `psi-min` = "0")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2L) stop("Use --name value pairs")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  if (is.null(opt$fixtures)) stop("Supply --fixtures with the recovered fixtures.rds")
  if (xor(is.null(opt$inputs), is.null(opt$preparation)))
    stop("Supply both --inputs and --preparation for the optional bundled grid")
  if (dir.exists(opt$output) && length(list.files(opt$output, all.files = TRUE, no.. = TRUE)))
    stop("Use a new output directory to preserve earlier receipts")
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output)
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library), .libPaths()))
  library(drfarm)
  source("tools/gaussian-objective.R", local = TRUE)
  control <- list(max.iter = as.numeric(opt[["max-iter"]]), objective.tol = 1e-8,
                  score.tol = 1e-6, psi.min = as.numeric(opt[["psi-min"]]))
  coefficient.control <- list(tol = 1e-10, max.sweeps = 1000L,
                              root.tol = 1e-14, root.maxit = 200L)
  log <- file(file.path(opt$output, "console.log"), "wt")
  sink(log, split = TRUE)
  on.exit({ sink(); close(log) }, add = TRUE)
  start <- proc.time(); stages <- list(); outputs <- list(); summary <- list()
  persist <- function() {
    if (length(stages)) write.table(do.call(rbind, stages), file.path(opt$output, "stages.tsv"),
      sep = "\t", row.names = FALSE, quote = TRUE)
    if (length(summary)) write.table(do.call(rbind, summary), file.path(opt$output, "comparison.tsv"),
      sep = "\t", row.names = FALSE, quote = TRUE)
    saveRDS(outputs, file.path(opt$output, "results.rds"))
  }
  stage <- function(name, expr) {
    cat("START ", name, "\n", sep = ""); flush.console()
    t0 <- proc.time(); warnings <- character(); error <- NULL
    value <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    elapsed <- proc.time() - t0
    stages[[length(stages) + 1L]] <<- data.frame(stage = name,
      execution_status = if (is.null(error)) "PASS" else "FAIL",
      elapsed_seconds = unname(elapsed["elapsed"]), user_seconds = unname(elapsed["user.self"]),
      system_seconds = unname(elapsed["sys.self"]), warnings = paste(warnings, collapse = " | "),
      error = if (is.null(error)) "" else error,
      result_bytes = if (is.null(value)) NA_real_ else as.numeric(object.size(value)))
    persist()
    cat(if (is.null(error)) "PASS " else "FAIL ", name, " ",
        sprintf("%.3fs", elapsed["elapsed"]), "\n", sep = "")
    if (!is.null(error)) cat("  ", error, "\n", sep = "")
    if (length(warnings)) cat(paste0("  WARNING: ", warnings), sep = "\n")
    flush.console()
    list(value = value, warnings = warnings, error = error, elapsed = elapsed)
  }
  # Independent dense marginal Gaussian derivatives, not implementation helpers.
  # Raw norms have physical units; these are not the API's scaled stopping score.
  observed_score <- function(a, fit) {
    theta <- fit$Theta; b <- fit$B; psi <- as.numeric(fit$diag.Psi)
    mask <- a$C
    if (is.null(mask)) mask <- matrix(1L, nrow(theta), ncol(theta))
    residual <- a$Y - a$X %*% t(theta)
    W <- solve(tcrossprod(b) + diag(psi, ncol(a$Y)))
    G <- (nrow(a$Y) * W - W %*% crossprod(residual) %*% W) / 2
    gtheta <- -t(crossprod(a$X, residual) %*% W)
    coefficient <- vapply(seq_len(ncol(theta)), function(j) {
      penalized <- mask[, j] == 1L; free <- mask[, j] == 2L
      u <- theta[penalized, j]; g <- gtheta[penalized, j]
      norm_u <- sqrt(sum(u^2))
      if (norm_u == 0) {
        r <- max(0, sqrt(sum(pmax(abs(g) - a$lambda1, 0)^2)) - a$lambda2)
      } else {
        v <- g + a$lambda2 * u / norm_u
        v[u != 0] <- v[u != 0] + a$lambda1 * sign(u[u != 0])
        v[u == 0] <- sign(v[u == 0]) * pmax(abs(v[u == 0]) - a$lambda1, 0)
        r <- sqrt(sum(v^2))
      }
      sqrt(r^2 + sum(gtheta[free, j]^2))
    }, numeric(1))
    gpsi <- diag(G); projected <- gpsi
    if (control$psi.min > 0) projected[psi == control$psi.min] <-
      pmin(projected[psi == control$psi.min], 0)
    list(coefficient = max(coefficient), factor = sqrt(sum((2 * G %*% b)^2)),
         variance = max(abs(projected)), raw.variance = gpsi,
         feasible = all(theta[mask == 0] == 0) && all(psi > 0) && all(psi >= control$psi.min))
  }
  # One plain actual DrFARM closure with a fixed fa return. No namespace edits.
  fixed_drfarm <- function(a, mode) {
    f <- drfarm::DrFARM.one
    e <- new.env(parent = environment(f))
    e$fa <- function(...) list(loadings = a$B, uniquenesses = a$psi)
    environment(f) <- e
    do.call(f, list(X = a$X, Y = a$Y, Theta0 = a$Theta, precM = a$precision,
      k = ncol(a$B), lambda1 = a$lambda1, lambda2 = a$lambda2, C = a$C,
      standardize = FALSE, max.iter = control$max.iter, thres = 1e-8,
      coefficient.update = mode,
      weighted.control = if (mode == "weighted") coefficient.control else list()))
  }
  fit_ecm <- function(a) gaussian.ecm.reference(a$X, a$Y, a$Theta, a$B, a$psi,
    a$lambda1, a$lambda2, C = a$C, control = control,
    coefficient.control = coefficient.control)
  evaluate <- function(name, mode, a, call) {
    run <- stage(paste(name, mode, sep = "/"), call)
    fit <- run$value
    if (is.null(fit)) {
      row <- data.frame(case = name, method = mode, initial.objective = NA_real_,
        observed.objective = NA_real_, objective.change = NA_real_, termination = "execution_error",
        reported.converged = FALSE, accepted.iterations = NA_integer_, attempted.iterations = NA_integer_,
        api.scaled.score = NA_real_, independent.per.observation.score = NA_real_, raw.coefficient.kkt = NA_real_, raw.factor.score = NA_real_,
        raw.variance.score = NA_real_, feasible = NA, elapsed.seconds = unname(run$elapsed["elapsed"]))
    } else {
      initial <- gaussian_observed(a$X, a$Y, a$Theta, a$B, a$psi, a$lambda1, a$lambda2, a$C)
      final <- gaussian_observed(a$X, a$Y, fit$Theta, fit$B, fit$diag.Psi, a$lambda1, a$lambda2, a$C)
      score <- observed_score(a, fit)
      d <- fit$diagnostics
      row <- data.frame(case = name, method = mode, initial.objective = initial,
        observed.objective = final, objective.change = final - initial,
        termination = d$termination, reported.converged = d$converged,
        accepted.iterations = if (mode == "gaussian.ecm.reference") d$iterations else NA_integer_,
        attempted.iterations = if (mode == "gaussian.ecm.reference") d$attempted else d$iterations,
        api.scaled.score = if (mode == "gaussian.ecm.reference") d$stationarity$maximum else NA_real_,
        independent.per.observation.score = max(score$coefficient, score$factor, score$variance) / nrow(a$Y),
        raw.coefficient.kkt = score$coefficient, raw.factor.score = score$factor,
        raw.variance.score = score$variance, feasible = score$feasible,
        elapsed.seconds = unname(run$elapsed["elapsed"]))
      run$independent.objective <- final; run$independent.score <- score
      if (mode == "gaussian.ecm.reference") stopifnot(
        abs(final - d$objective) <= 1e-9 * (1 + abs(final)),
        abs(row$independent.per.observation.score - d$stationarity$maximum) <=
          1e-9 * (1 + row$independent.per.observation.score))
    }
    outputs[[paste(name, mode, sep = "/")]] <<- run
    summary[[length(summary) + 1L]] <<- row
    persist()
    invisible(run)
  }
  fixtures <- readRDS(opt$fixtures)
  saveRDS(fixtures, file.path(opt$output, "fixtures.rds"))
  for (name in names(fixtures)) {
    a <- fixtures[[name]]
    for (mode in c("historical", "weighted", "gaussian.ecm.reference"))
      evaluate(name, mode, a, if (mode == "gaussian.ecm.reference") fit_ecm(a) else fixed_drfarm(a, mode))
  }
  common <- NULL
  if (!is.null(opt$inputs)) {
    inputs <- readRDS(opt$inputs); prep <- readRDS(opt$preparation)
    common <- stage("bundled/preparation_and_one_factor_initialization", {
      X <- scale(inputs$X); Y <- scale(inputs$Y); Theta <- prep$initial$Theta0
      theta_db <- t(Theta) + prep$precision %*% crossprod(X, Y - X %*% t(Theta)) / nrow(X)
      set.seed(20260909L)
      fa <- psych::fa(Y - X %*% theta_db, nfactors = 2, rotate = "none",
                      scores = "regression", fm = "ml", covar = TRUE)
      list(X = X, Y = Y, Theta = Theta, B = unclass(as.matrix(fa$loadings)),
           psi = fa$uniquenesses, precision = prep$precision, grid = prep$grid,
           initialization = "Historical debiased initializer only; reused once for all methods and cells")
    })$value
    if (!is.null(common)) {
      saveRDS(common, file.path(opt$output, "bundled-common-start.rds"))
      # First/second passes are cold/warm only within this loaded process;
      # OS/backend cache state is uncontrolled. Each mode restarts at common.
      for (pass in c("cold", "warm")) for (mode in c("historical", "weighted", "gaussian.ecm.reference")) {
        profile_name <- paste(mode, pass, sep = "-")
        Rprof(file.path(opt$output, paste0(profile_name, ".Rprof")), interval = 0.01,
              memory.profiling = TRUE)
        for (i in c(1L, seq_len(nrow(common$grid)))) {
          a <- common; a$lambda1 <- common$grid$lambda1[i]; a$lambda2 <- common$grid$lambda2[i]
          kind <- if (i == 1L && !exists("did_one", inherits = FALSE)) "one" else paste0("grid", i)
          did_one <- TRUE
          evaluate(paste("bundled", pass, kind, sep = "/"), mode, a,
                   if (mode == "gaussian.ecm.reference") fit_ecm(a) else fixed_drfarm(a, mode))
        }
        rm(did_one)
        Rprof(NULL)
        profile <- tryCatch(summaryRprof(file.path(opt$output, paste0(profile_name, ".Rprof")),
                                          memory = "both"), error = function(e) list(error = conditionMessage(e)))
        saveRDS(profile, file.path(opt$output, paste0(profile_name, "-profile.rds")))
        if (!is.null(profile$by.total)) write.table(profile$by.total,
          file.path(opt$output, paste0(profile_name, "-inclusive-functions.tsv")), sep = "\t",
          row.names = TRUE, col.names = NA)
      }
    }
  }
  bindings <- c(fixtures = opt$fixtures, inputs = opt$inputs, preparation = opt$preparation)
  receipt <- list(command = commandArgs(), source_commit = Sys.getenv("DRFARM_SOURCE_COMMIT"),
    package_path = find.package("drfarm"), package_version = as.character(packageVersion("drfarm")),
    session = sessionInfo(), source_files_md5 = tools::md5sum(bindings),
    control = control, coefficient.control = coefficient.control, seed = 20260909L,
    input_generation = "None; recovered fixed fixtures and optional bundled/preparation assets",
    comparison_scope = "Common supplied starts. Preserved DrFARM controls differ in monitored target and stopping meaning. No cross-method selection or inference.",
    profile_scope = "Rprof inclusive function costs overlap and include comparison/reporting/serialization overhead. Fit timing excludes those costs. One fit plus four cells, two passes per method; setup separately. No speedup claim.",
    prior_fixture_max_iter = lapply(fixtures, `[[`, "max.iter"),
    current_fixture_max_iter = control$max.iter,
    inference = "NOT RUN; Gaussian ECM is a separate optimizer and inherits no DrFARM inference claims",
    threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")),
    elapsed = proc.time() - start, memory = grep("^Vm(HWM|RSS):", readLines("/proc/self/status"), value = TRUE),
    gc = gc(), execution_failures = sum(vapply(stages, function(x) x$execution_status == "FAIL", logical(1))))
  saveRDS(receipt, file.path(opt$output, "receipt.rds"))
  writeLines(capture.output(str(receipt, max.level = 3)), file.path(opt$output, "receipt.txt"))
  writeLines(capture.output(sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
  print(do.call(rbind, summary), row.names = FALSE)
  if (receipt$execution_failures == 0) 0L else 1L
}
quit(status = main(), save = "no")
