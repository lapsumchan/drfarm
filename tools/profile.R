#!/usr/bin/env Rscript
# Matched cold/warm session runs. No result cache, optimization, or speed claim.
# Rprof inclusive function times overlap; they are not additive stage costs.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(output = "profile", library = NULL, kinship = "no")
  if (length(args) %% 2L) stop("Arguments must be --name value pairs")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  stopifnot(opt$kinship %in% c("yes", "no"))
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output)
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library), .libPaths()))
  log <- file(file.path(opt$output, "console.log"), "wt")
  sink(log, split = TRUE)
  on.exit({ sink(); close(log) }, add = TRUE)
  rows <- list(); results <- list(); counts <- list(); fit_status <- list()
  stage <- function(run, name, expr) {
    warnings <- character(); error <- NULL
    start <- proc.time()
    result <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    timing <- proc.time() - start
    rows[[length(rows) + 1L]] <<- data.frame(
      run = run, stage = name, status = if (is.null(error)) "PASS" else "FAIL",
      elapsed_seconds = unname(timing["elapsed"]), user_seconds = unname(timing["user.self"]),
      system_seconds = unname(timing["sys.self"]),
      result_object_bytes = if (is.null(error)) as.numeric(object.size(result)) else NA_real_,
      warning_count = length(warnings), warnings = paste(warnings, collapse = " | "),
      error = if (is.null(error)) "" else error)
    write.table(do.call(rbind, rows), file.path(opt$output, "stages.tsv"),
                sep = "\t", row.names = FALSE, quote = TRUE)
    cat(run, " ", name, " ", tail(rows, 1L)[[1]]$status, " ",
        sprintf("%.3fs", timing["elapsed"]), "\n", sep = "")
    flush.console()
    if (!is.null(error)) stop(error, call. = FALSE)
    result
  }
  stage("setup", "package_load", library("drfarm", character.only = TRUE))
  writeLines(capture.output(sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
  metadata <- list(seed = 20260909L, RNGkind = RNGkind(), session = sessionInfo(),
                   source_commit = Sys.getenv("DRFARM_SOURCE_COMMIT", unset = NA_character_),
                   package_path = find.package("drfarm"), sys = Sys.info(),
                   threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")),
                   options = list(standardize = TRUE, K = NULL, thres = 1e-4, k = 2,
                                  rotate = "none", scores = "regression", fm = "ml"),
                   cold_warm = "First versus second complete path in one loaded R process; seed reset, no model-result cache or input reuse; OS/backend caches uncontrolled.",
                   profile = "Rprof 0.01s sampling; overlapping inclusive costs. Native work is attributed to its R caller. Timing includes profiling overhead.",
                   memory = "Result object sizes are not process memory; GC high-water marks cover this R process; external time -v gives process peak RSS when available.",
                   interpretation = "One fit plus four-cell diagnostic grid; not the original 25-cell selection, a scaling study, or a speedup benchmark.")
  metadata$installed_data_file_md5 <- tools::md5sum(list.files(
    system.file("data", package = "drfarm"), full.names = TRUE))
  status <- 0L
  for (run in c("cold", "warm")) {
    set.seed(metadata$seed)
    profile_path <- file.path(opt$output, paste0(run, ".Rprof"))
    Rprof(profile_path, interval = 0.01, memory.profiling = TRUE)
    run_start <- proc.time()
    result <- tryCatch({
      dat <- stage(run, "data_preparation", {
        env <- new.env(parent = emptyenv())
        data("drfarm.dat", package = "drfarm", envir = env)
        env$drfarm.dat
      })
      X <- dat$X; Y <- dat$Y
      metadata$dimensions <- c(n = nrow(X), p = ncol(X), q = ncol(Y), k = 2L)
      init <- stage(run, "remMap_initialization_2x2", remMap.whole(X, Y, n.lambda = 2))
      precision <- stage(run, "precision_glasso", precM(X))
      full_grid <- stage(run, "factor_initialization_grid", DrFARM.grid(
        X, Y, init$Theta0, precision, 2, init$lambda1.opt, init$lambda2.opt))
      # Preserve first-occurrence ordering of each source-grid lambda axis.
      pick <- function(x) { x <- unique(x); x[unique(c(1L, length(x)))] }
      grid <- expand.grid(lambda1 = pick(full_grid$lambda1), lambda2 = pick(full_grid$lambda2))
      if (nrow(grid) != 4L) stop("Source grid does not contain two distinct points on both axes")
      counts[[run]] <- list(remMap_grid_cells = 4L, DrFARM_grid_cells = nrow(grid),
                            one_fit_extra = 1L, grid = grid)
      fit_call <- function(i) DrFARM.one(X, Y, init$Theta0, precision, 2,
                                        grid[i, 1], grid[i, 2])
      one <- stage(run, "one_fit_including_initialization_and_iterations", fit_call(1L))
      fits <- vector("list", nrow(grid))
      errors <- character(nrow(grid))
      for (i in seq_len(nrow(grid))) {
        fits[i] <- list(tryCatch(stage(run, paste0("grid_fit_", i), fit_call(i)),
                                 error = function(e) { errors[i] <<- conditionMessage(e); NULL }))
      }
      fit_status[[run]] <- lapply(c(list(one), fits), function(x) {
        if (is.null(x)) "FAILED"
        else if (!is.null(x$diagnostics)) x$diagnostics
        else "Returned fit; historical API does not expose convergence."
      })
      if (any(nzchar(errors))) stop("One or more grid fits failed; selection and inference NOT RUN (no failed-cell exclusion).")
      selection <- stage(run, "grid_selection", {
        ebic <- vapply(fits, function(f) DrFARM.EBIC(X, Y, f$Theta, f$B, f$E.Z, f$diag.Psi), numeric(1))
        if (any(!is.finite(ebic))) stop("Nonfinite EBIC; selection NOT RUN")
        list(EBIC = ebic, index = which.min(ebic), lambdas = grid[which.min(ebic), ])
      })
      fit <- fits[[selection$index]]
      inference <- stage(run, "entry_and_historical_predictor_inference", list(
        entry = entry.pvalue(X, Y, fit$Theta, fit$B, fit$E.Z, precision),
        pleio = pleio.pvalue(X, Y, fit$Theta, fit$B, fit$E.Z, precision)))
      value <- list(one = one, grid = grid, fits = fits, selection = selection,
                    inference = inference, precision = precision, initial = init,
                    input_file_md5 = metadata$installed_data_file_md5)
      stage(run, "serialization", saveRDS(value, file.path(opt$output, paste0(run, "-results.rds"))))
      value
    }, error = function(e) {
      status <<- 1L
      list(error = conditionMessage(e), skipped_stages = "Subsequent stages NOT RUN; consult ordered stages.tsv")
    })
    Rprof(NULL)
    timing <- proc.time() - run_start
    result$whole_path_seconds <- timing
    results[[run]] <- result
    profile <- tryCatch(summaryRprof(profile_path, memory = "both"), error = function(e) list(error = conditionMessage(e)))
    saveRDS(profile, file.path(opt$output, paste0(run, "-profile.rds")))
    if (!is.null(profile$by.total)) {
      write.table(profile$by.total, file.path(opt$output, paste0(run, "-inclusive-functions.tsv")),
                  sep = "\t", row.names = TRUE, col.names = NA, quote = TRUE)
      # fa includes initialization; solve/eigen/svd decompose/solve. DrFARM.one
      # includes initialization and iterations; its cost cannot be subtracted
      # from these overlapping samples to infer an exact iteration-only cost.
      roles <- c(factor_analysis = '"fa"', solves = '"solve"', eigen = '"eigen"',
                 svd = '"svd"', iteration_kernel = '"remMap"',
                 fit_inclusive = '"DrFARM.one"', precision = '"precM"')
      attribution <- data.frame(role = names(roles), function_name = unname(roles),
        inclusive_sample_seconds = vapply(roles, function(f) {
          if (f %in% rownames(profile$by.total)) profile$by.total[f, "total.time"] else NA_real_
        }, numeric(1)),
        scope = "Overlapping samples; NA means not sampled, not zero; native internal iteration cost unresolved")
      write.table(attribution, file.path(opt$output, paste0(run, "-major-stages.tsv")),
                  sep = "\t", row.names = FALSE, quote = TRUE)
    }
  }
  if (opt$kinship == "yes") {
    # Diagnostic only: compare the evaluator's rotated participant basis with
    # the historical whole-fit inline score's original participant basis.
    # No algorithm correction or p-value validity claim is made here.
    metadata$kinship <- list()
    for (kind in c("identity", "unequal_spectrum")) {
      set.seed(metadata$seed)
      metadata$kinship[[kind]] <- tryCatch(stage("kinship", kind, {
        env <- new.env(parent = emptyenv())
        data("drfarm.dat", package = "drfarm", envir = env)
        X <- env$drfarm.dat$X[seq_len(80L), , drop = FALSE]
        Y <- env$drfarm.dat$Y[seq_len(80L), , drop = FALSE]
        n <- nrow(X); q <- ncol(Y)
        K <- if (kind == "identity") diag(n) else {
          toeplitz(0.4 ^ (0:(n - 1L))) + diag(seq(0.1, 0.9, length.out = n))
        }
        init <- remMap.whole(X, Y, n.lambda = 2)
        precision <- precM(X)
        args <- list(X = X, Y = Y, Theta0 = init$Theta0, precM = precision,
                     k = 2, lambda1 = init$lambda1.opt, lambda2 = init$lambda2.opt, K = K)
        if ("max.iter" %in% names(formals(DrFARM.one))) args$max.iter <- 1000L
        fit <- do.call(DrFARM.one, args)
        Xs <- scale(X); Ys <- scale(Y)
        eg <- eigen(K); U <- eg$vectors
        score <- function(x, y) {
          E <- y - x %*% t(fit$Theta) - fit$E.Z %*% t(fit$B)
          # Independent integer-combination reference for small q, including
          # the historical log(sum choose) penalty rather than sum(log choose).
          sum(colSums(E ^ 2) / fit$diag.Psi) + n * sum(log(fit$diag.Psi)) +
            log(n) * sum(fit$Theta != 0) +
            2 * log(sum(choose(q, rowSums(t(fit$Theta) != 0))))
        }
        rotated <- score(crossprod(U, Xs), crossprod(U, Ys))
        unrotated <- score(Xs, Ys)
        api <- DrFARM.EBIC(X, Y, fit$Theta, fit$B, fit$E.Z, fit$diag.Psi, K = K)
        stopifnot(abs(api - rotated) <= 1e-10 * (1 + abs(rotated)))
        out <- list(n = n, p = ncol(X), q = q, K = K,
                    eigenvalue_range = range(eg$values), fit = fit,
                    evaluator = api, independently_rotated = rotated,
                    historical_whole_inline_unrotated = unrotated,
                    unrotated_minus_rotated = unrotated - rotated,
                    scope = "One fixed fit: evaluator matches explicit rotation; a score difference is not evidence that the selected grid index changes.")
        saveRDS(out, file.path(opt$output, paste0("kinship-", kind, ".rds")))
        out[c("n", "p", "q", "eigenvalue_range", "evaluator", "independently_rotated",
              "historical_whole_inline_unrotated", "unrotated_minus_rotated", "scope")]
      }), error = function(e) list(status = "FAIL", error = conditionMessage(e)))
      if (identical(metadata$kinship[[kind]]$status, "FAIL")) status <- 1L
    }
  }
  metadata$gc_high_water <- gc()
  metadata$linux_process_memory <- if (file.exists("/proc/self/status")) {
    grep("^Vm(HWM|RSS):", readLines("/proc/self/status"), value = TRUE)
  } else "NOT RUN: Linux /proc memory counters unavailable"
  metadata$fit_status <- fit_status
  metadata$counts <- counts
  if (status == 0L) {
    invariant <- function(x) list(Theta = x$one$Theta,
      factor_product = x$one$E.Z %*% t(x$one$B),
      cov = tcrossprod(x$one$B) + diag(as.numeric(x$one$diag.Psi)),
      selection = x$selection, inference = x$inference)
    metadata$cold_warm_agreement <- all.equal(invariant(results$cold), invariant(results$warm), tolerance = 1e-10)
    if (!isTRUE(metadata$cold_warm_agreement)) status <- 1L
  }
  metadata$status <- if (status == 0L) "PASS" else "FAIL"
  metadata$whole_path_seconds <- lapply(results, `[[`, "whole_path_seconds")
  saveRDS(metadata, file.path(opt$output, "receipt.rds"))
  writeLines(capture.output(str(metadata)), file.path(opt$output, "receipt.txt"))
  print(do.call(rbind, rows))
  status
}
quit(status = main(), save = "no")
