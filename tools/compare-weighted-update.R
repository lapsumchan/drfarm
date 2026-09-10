#!/usr/bin/env Rscript
# A bounded fixed-sigma solver/outer-fit comparison, not a simulation study.
# Reuses the bundled-data preparation, invariants and Rprof receipt conventions
# from tools/profile.R at maintenance commit cbb1e8c9bacbd0ba5ae3b4bd5b8a2d1d2f7481dd.
# The historical and weighted paths intentionally optimize different coefficient
# updates. Agreement between these paths is neither expected nor required.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(output = "weighted-comparison", library = NULL, `max-outer` = "1000")
  if (length(args) %% 2L) stop("Arguments must be --name value pairs")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  max_outer <- as.numeric(opt[["max-outer"]])
  stopifnot(is.finite(max_outer), max_outer > 0, max_outer == floor(max_outer))
  if (dir.exists(opt$output) && length(list.files(opt$output, all.files = TRUE, no.. = TRUE))) {
    stop("Use a new output directory to preserve previous receipts")
  }
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output)
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library), .libPaths()))
  log <- file(file.path(opt$output, "console.log"), "wt")
  sink(log, split = TRUE)
  on.exit({ sink(); close(log) }, add = TRUE)
  seed <- 20260909L
  control <- list(tol = 1e-8, max.sweeps = 1000L, root.tol = 1e-14, root.maxit = 200L)
  rows <- list(); results <- list(); failure_count <- 0L
  persist <- function() {
    if (length(rows)) write.table(do.call(rbind, rows), file.path(opt$output, "stages.tsv"),
                                  sep = "\t", quote = TRUE, row.names = FALSE)
  }
  stage <- function(run, name, expr) {
    cat("START ", run, "/", name, "\n", sep = "")
    flush.console()
    start <- proc.time(); warnings <- character(); error <- NULL
    value <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      cat("WARNING ", run, "/", name, ": ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    timing <- proc.time() - start
    rows[[length(rows) + 1L]] <<- data.frame(
      run = run, stage = name, execution_status = if (is.null(error)) "PASS" else "FAIL",
      elapsed_seconds = unname(timing["elapsed"]), user_seconds = unname(timing["user.self"]),
      system_seconds = unname(timing["sys.self"]), warning_count = length(warnings),
      warning = paste(warnings, collapse = " | "),
      error = if (is.null(error)) "" else error,
      result_object_bytes = if (is.null(error)) as.numeric(object.size(value)) else NA_real_)
    persist()
    cat(run, "/", name, ": ", if (is.null(error)) "PASS" else "FAIL", " (",
        sprintf("%.3fs", timing["elapsed"]), ")\n", sep = "")
    flush.console()
    if (!is.null(error)) {
      failure_count <<- failure_count + 1L
      stop(error, call. = FALSE)
    }
    value
  }
  stage("setup", "package_load", library("drfarm", character.only = TRUE))
  metadata <- list(
    command = commandArgs(), source_commit = Sys.getenv("DRFARM_SOURCE_COMMIT", unset = NA_character_),
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    session = sessionInfo(), package_path = find.package("drfarm"),
    package_version = as.character(packageVersion("drfarm")), seed = seed, RNGkind = RNGkind(),
    weighted_control = control, max_outer = max_outer, outer_threshold = 1e-4,
    threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                          "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")),
    environment = Sys.info(),
    profile_scope = "Common preparation measured once; one extra fit plus a four-cell grid per mode and pass. Rprof 0.01s inclusive costs overlap and include profiling overhead. No speedup claim.",
    cold_warm_scope = "First and second traversal per mode in one loaded process, historical then weighted each pass. Seed reset before every traversal; common inputs/initialization reused exactly, no fitted-output cache. OS/backend caches uncontrolled.",
    inference = "NOT RUN: this comparison validates a fixed-sigma coefficient subproblem and records changed outer behavior; it does not validate weighted DrFARM inference.",
    optimizer_scope = "Execution PASS is distinct from outer convergence. Retain loss_increase/max_iter and all weighted inner diagnostics. The unchanged outer loss is not a likelihood or KKT certificate.")
  writeLines(capture.output(sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
  inputs <- stage("setup", "bundled_data", {
    env <- new.env(parent = emptyenv())
    data("drfarm.dat", package = "drfarm", envir = env)
    env$drfarm.dat
  })
  X <- inputs$X; Y <- inputs$Y
  metadata$dimensions <- c(n = nrow(X), p = ncol(X), q = ncol(Y), k = 2L)
  saveRDS(inputs, file.path(opt$output, "inputs.rds"))
  metadata$input_md5 <- unname(tools::md5sum(file.path(opt$output, "inputs.rds")))
  set.seed(seed)
  initial <- stage("setup", "historical_remMap_initialization_2x2", remMap.whole(X, Y, n.lambda = 2))
  precision <- stage("setup", "precision_glasso", precM(X))
  full_grid <- stage("setup", "factor_initialization_grid", DrFARM.grid(
    X, Y, initial$Theta0, precision, 2, initial$lambda1.opt, initial$lambda2.opt))
  pick <- function(x) { x <- unique(x); x[unique(c(1L, length(x)))] }
  grid <- expand.grid(lambda1 = pick(full_grid$lambda1), lambda2 = pick(full_grid$lambda2))
  stopifnot(nrow(grid) == 4L)
  metadata$grid <- grid
  saveRDS(list(initial = initial, precision = precision, grid = grid),
          file.path(opt$output, "common-preparation.rds"))

  # Existing failure fixture, re-derived from its convex first-order equations.
  # X'X=1 and no lasso: b_j=z_j/(d_j+t), d_j=1/sigma_j,
  # z_j=a_j/sigma_j, t=lambda_group/||b||. Solve this scalar secular equation.
  counterexample <- tryCatch(stage("fixed_sigma", "unequal_variance_counterexample", {
    x <- matrix(c(1, -1) / sqrt(2), 2, 1)
    a <- c(3, 4); sigma <- c(1, 2); lambda <- 1
    y <- x %*% matrix(a, 1, 2)
    d <- 1 / sigma; z <- a / sigma
    root <- uniroot(function(t) sum((t * z / (d + t)) ^ 2) - lambda ^ 2,
                    interval = c(0, 10), tol = 1e-12)
    exact_reference <- z / (d + root$root)
    old <- drfarm:::remMap(x, y, lamL1 = 0, lamL2 = lambda, sigma = sigma)
    new <- remMap.weighted(x, y, lambda1 = 0, lambda2 = lambda,
                           sigma = sigma, control = control)
    gradient <- function(b) (b - a) / sigma + lambda * b / sqrt(sum(b ^ 2))
    objective <- function(b) sum((b - a) ^ 2 / sigma) / 2 + lambda * sqrt(sum(b ^ 2))
    old_b <- as.numeric(old$phi); new_b <- as.numeric(new$Theta0)
    tab <- data.frame(
      method = c("historical", "weighted", "independent_secular_root"),
      b1 = c(old_b[1], new_b[1], exact_reference[1]),
      b2 = c(old_b[2], new_b[2], exact_reference[2]),
      weighted_objective = vapply(list(old_b, new_b, exact_reference), objective, numeric(1)),
      gradient_infinity = vapply(list(old_b, new_b, exact_reference),
                                 function(b) max(abs(gradient(b))), numeric(1)))
    write.table(tab, file.path(opt$output, "counterexample.tsv"), sep = "\t",
                quote = TRUE, row.names = FALSE)
    result <- list(table = tab, historical_diagnostics = old$diagnostics,
                   weighted_diagnostics = new$diagnostics, root = root,
                   max_coefficient_difference = max(abs(new_b - exact_reference)),
                   contract = "X'X=1; sigma=(1,2); lambda1=0; lambda2=1; all entries penalized; unique convex minimizer.")
    saveRDS(result, file.path(opt$output, "counterexample.rds"))
    stopifnot(isTRUE(new$diagnostics$converged),
              max(abs(new_b - exact_reference)) <= 1e-7,
              max(abs(gradient(new_b))) <= 1e-7,
              max(abs(gradient(old_b))) > 0.1,
              objective(new_b) <= objective(old_b) + 1e-10)
    print(tab)
    result
  }), error = function(e) list(error = conditionMessage(e)))
  metadata$counterexample <- counterexample

  summarize_fit <- function(fit) {
    if (is.null(fit)) return(list(execution = "FAIL"))
    history <- fit$diagnostics$coefficient.history
    list(outer = fit$diagnostics,
         max_inner_kkt = if (length(history)) max(vapply(history, `[[`, numeric(1), "kkt.residual")) else NA_real_,
         max_inner_scaled_kkt = if (length(history)) max(vapply(history, `[[`, numeric(1), "kkt.scaled")) else NA_real_,
         inner_failures = if (length(history)) sum(!vapply(history, `[[`, logical(1), "converged")) else NA_integer_)
  }
  for (pass in c("cold", "warm")) for (mode in c("historical", "weighted")) {
    run <- paste(mode, pass, sep = "-")
    set.seed(seed)
    start <- proc.time()
    Rprof(file.path(opt$output, paste0(run, ".Rprof")), interval = 0.01, memory.profiling = TRUE)
    result <- tryCatch({
      call_fit <- function(i) {
        fit_args <- list(X = X, Y = Y, Theta0 = initial$Theta0, precM = precision,
          k = 2, lambda1 = grid[i, 1], lambda2 = grid[i, 2],
          standardize = TRUE, thres = 1e-4, max.iter = max_outer,
          coefficient.update = mode)
        if (mode == "weighted") fit_args$weighted.control <- control
        do.call(DrFARM.one, fit_args)
      }
      one <- stage(run, "one_fit_including_initialization_and_iterations", call_fit(1L))
      fits <- vector("list", nrow(grid)); errors <- character(nrow(grid))
      for (i in seq_len(nrow(grid))) {
        fits[i] <- list(tryCatch(stage(run, paste0("grid_fit_", i), call_fit(i)),
                                 error = function(e) { errors[i] <<- conditionMessage(e); NULL }))
      }
      saveRDS(list(one = one, fits = fits, errors = errors),
              file.path(opt$output, paste0(run, "-fits.rds")))
      if (any(nzchar(errors))) stop("Grid fit error: selection NOT RUN; failed candidates are not excluded.")
      selection <- stage(run, "historical_EBIC_selection", {
        ebic <- vapply(fits, function(f) DrFARM.EBIC(X, Y, f$Theta, f$B, f$E.Z, f$diag.Psi), numeric(1))
        stopifnot(all(is.finite(ebic)))
        list(EBIC = ebic, index = which.min(ebic), lambda = grid[which.min(ebic), ])
      })
      out <- list(one = one, fits = fits, selection = selection,
                   diagnostics = lapply(c(list(one), fits), summarize_fit))
      stage(run, "serialization", saveRDS(out, file.path(opt$output, paste0(run, "-results.rds"))))
      out
    }, error = function(e) list(error = conditionMessage(e)))
    Rprof(NULL)
    result$whole_path_seconds <- proc.time() - start
    result$linux_process_memory <- if (file.exists("/proc/self/status")) {
      grep("^Vm(HWM|RSS):", readLines("/proc/self/status"), value = TRUE)
    } else "NOT RUN: Linux /proc unavailable"
    results[[run]] <- result
    profile <- tryCatch(summaryRprof(file.path(opt$output, paste0(run, ".Rprof")), memory = "both"),
                         error = function(e) list(error = conditionMessage(e)))
    saveRDS(profile, file.path(opt$output, paste0(run, "-profile.rds")))
    if (!is.null(profile$by.total)) write.table(profile$by.total,
      file.path(opt$output, paste0(run, "-inclusive-functions.tsv")), sep = "\t",
      row.names = TRUE, col.names = NA, quote = TRUE)
    saveRDS(results, file.path(opt$output, "results.rds"))
  }
  invariant <- function(result) {
    list(Theta = result$one$Theta, factor_product = result$one$E.Z %*% t(result$one$B),
         response_covariance = tcrossprod(result$one$B) + diag(as.numeric(result$one$diag.Psi)),
         selection = result$selection)
  }
  metadata$cold_warm_agreement <- lapply(c("historical", "weighted"), function(mode) {
    cold <- results[[paste(mode, "cold", sep = "-")]]
    warm <- results[[paste(mode, "warm", sep = "-")]]
    if (!is.null(cold$error) || !is.null(warm$error)) return("NOT RUN: failed fit path")
    all.equal(invariant(cold), invariant(warm), tolerance = 1e-10)
  })
  names(metadata$cold_warm_agreement) <- c("historical", "weighted")
  metadata$run_summary <- lapply(results, function(result) {
    result[c("whole_path_seconds", "linux_process_memory", "error", "diagnostics", "selection")]
  })
  if (all(vapply(results, function(z) is.null(z$error), logical(1)))) {
    a <- results[["historical-cold"]]$one
    b <- results[["weighted-cold"]]$one
    metadata$between_method_changes <- list(
      one_Theta_max_abs = max(abs(a$Theta - b$Theta)),
      one_factor_product_max_abs = max(abs(a$E.Z %*% t(a$B) - b$E.Z %*% t(b$B))),
      historical_selected_row = results[["historical-cold"]]$selection$index,
      weighted_selected_row = results[["weighted-cold"]]$selection$index,
      scope = "Different coefficient algorithms; differences are reported, not equivalence failures. Historical EBIC and outer monitored loss retained.")
  }
  metadata$gc_high_water <- gc()
  metadata$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  metadata$stage_failures <- failure_count
  metadata$execution_status <- if (failure_count == 0L &&
    all(vapply(metadata$cold_warm_agreement, isTRUE, logical(1))) &&
    all(vapply(results, function(z) is.null(z$error), logical(1)))) "PASS" else "FAIL"
  saveRDS(metadata, file.path(opt$output, "receipt.rds"))
  writeLines(capture.output(str(metadata, max.level = 5)), file.path(opt$output, "receipt.txt"))
  cat("EXECUTION ", metadata$execution_status, "; outer statuses are reported separately.\n", sep = "")
  if (metadata$execution_status == "PASS") 0L else 1L
}
quit(status = main(), save = "no")
