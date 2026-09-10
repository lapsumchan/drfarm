#!/usr/bin/env Rscript
# Run against an explicitly installed package, never source() a mutable R file.
# Full mode preserves the original 100-cell remMap / 25-cell DrFARM example.
# Quick mode is a smoke example only; it is not the historical tuning workflow.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(output = "reproduction", mode = "full", library = NULL,
              compare = NULL, label = "unspecified", tolerance = "1e-10")
  if (length(args) %% 2L) stop("Arguments must be --name value pairs")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  stopifnot(opt$mode %in% c("quick", "full"))
  tolerance <- as.numeric(opt$tolerance)
  stopifnot(length(tolerance) == 1L, is.finite(tolerance), tolerance > 0)
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  opt$output <- normalizePath(opt$output)
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library), .libPaths()))
  log <- file(file.path(opt$output, "console.log"), "wt")
  sink(log, split = TRUE)
  on.exit({ sink(); close(log) }, add = TRUE)
  seed <- 20260909L
  set.seed(seed)
  rows <- list()
  receipt <- list(status = "RUNNING", label = opt$label, mode = opt$mode,
                  seed = seed, RNGkind = RNGkind(), command = commandArgs(),
                  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                  session = sessionInfo(), environment = Sys.info(),
                  threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                         "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")),
                  source_commit = Sys.getenv("DRFARM_SOURCE_COMMIT", unset = NA_character_),
                  interpretation = "Local reproducibility, not statistical validation")
  outputs <- list()
  flush_receipt <- function() {
    receipt$stages <- if (length(rows)) do.call(rbind, rows) else data.frame()
    saveRDS(receipt, file.path(opt$output, "receipt.rds"))
    if (length(rows)) write.table(receipt$stages, file.path(opt$output, "stages.tsv"),
                                  sep = "\t", row.names = FALSE, quote = TRUE)
    saveRDS(outputs, file.path(opt$output, "results.rds"))
  }
  stage <- function(name, expr) {
    cat("START ", name, "\n", sep = "")
    flush.console()
    elapsed <- proc.time()
    warnings <- character()
    error <- NULL
    result <- tryCatch(withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      cat("WARNING ", name, ": ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
    elapsed <- proc.time() - elapsed
    rows[[length(rows) + 1L]] <<- data.frame(
      stage = name, status = if (is.null(error)) "PASS" else "FAIL",
      elapsed_seconds = unname(elapsed["elapsed"]),
      user_seconds = unname(elapsed["user.self"]),
      system_seconds = unname(elapsed["sys.self"]),
      warning_count = length(warnings), warnings = paste(warnings, collapse = " | "),
      error = if (is.null(error)) "" else error)
    flush_receipt()
    cat(tail(rows, 1L)[[1]]$status, " ", name, " ",
        sprintf("%.3fs", elapsed["elapsed"]), "\n", sep = "")
    if (!is.null(error)) stop(error, call. = FALSE)
    result
  }
  status <- tryCatch({
    stage("package_load", library("drfarm", character.only = TRUE))
    receipt$session <- sessionInfo()
    receipt$package_path <- find.package("drfarm")
    receipt$package_version <- as.character(packageVersion("drfarm"))
    writeLines(capture.output(sessionInfo()), file.path(opt$output, "sessionInfo.txt"))
    installed <- installed.packages()
    write.table(installed[, c("Package", "Version", "LibPath")],
                file.path(opt$output, "installed-packages.tsv"),
                sep = "\t", row.names = FALSE, quote = TRUE)
    inputs <- stage("bundled_data", {
      env <- new.env(parent = emptyenv())
      data("drfarm.dat", package = "drfarm", envir = env)
      env$drfarm.dat
    })
    X <- inputs$X; Y <- inputs$Y
    stopifnot(identical(dim(inputs$Theta.t), c(ncol(X), ncol(Y))))
    receipt$dimensions <- c(n = nrow(X), p = ncol(X), q = ncol(Y), k = 2L)
    receipt$options <- list(standardize = TRUE, thres = 1e-4, K = NULL,
                            rotate = "none", scores = "regression", fm = "ml")
    outputs$inputs <- inputs
    saveRDS(inputs, file.path(opt$output, "inputs.rds"))
    receipt$input_md5 <- unname(tools::md5sum(file.path(opt$output, "inputs.rds")))
    outputs$preprocessing <- list(X_center = colMeans(X), X_scale = apply(X, 2, sd),
                                  Y_center = colMeans(Y), Y_scale = apply(Y, 2, sd))
    outputs$remMap <- stage("remMap_initialization", {
      if (opt$mode == "full") remMap.whole(X, Y)
      else remMap.whole(X, Y, n.lambda = 2)
    })
    outputs$precision <- stage("precision_glasso", precM(X))
    r <- outputs$remMap
    outputs$fit <- stage(if (opt$mode == "full") "DrFARM_grid_25" else "DrFARM_one", {
      if (opt$mode == "full") {
        DrFARM.whole(X, Y, r$Theta0, outputs$precision, 2,
                     r$lambda1.opt, r$lambda2.opt)
      } else {
        grid <- DrFARM.grid(X, Y, r$Theta0, outputs$precision, 2,
                            r$lambda1.opt, r$lambda2.opt)
        DrFARM.one(X, Y, r$Theta0, outputs$precision, 2, grid[1, 1], grid[1, 2])
      }
    })
    fit <- outputs$fit
    outputs$entry <- stage("entry_inference", entry.pvalue(X, Y, fit$Theta, fit$B,
                                                           fit$E.Z, outputs$precision))
    outputs$pleio <- stage("historical_predictor_inference", pleio.pvalue(
      X, Y, fit$Theta, fit$B, fit$E.Z, outputs$precision))
    outputs$invariants <- list(
      Theta0 = r$Theta0, precision = outputs$precision, Theta = fit$Theta,
      factor_product = fit$E.Z %*% t(fit$B),
      response_covariance = tcrossprod(fit$B) + diag(as.numeric(fit$diag.Psi)),
      diag_Psi = as.numeric(fit$diag.Psi), entry_p = outputs$entry,
      predictor_p = outputs$pleio,
      remMap_lambda = c(r$lambda1.opt, r$lambda2.opt),
      DrFARM_lambda = c(fit$lambda1.opt, fit$lambda2.opt))
    stage("shape_and_finiteness", {
      stopifnot(identical(dim(fit$Theta), c(ncol(Y), ncol(X))),
                identical(dim(fit$E.Z), c(nrow(X), 2L)),
                identical(dim(outputs$entry), c(ncol(Y), ncol(X))),
                length(outputs$pleio) == ncol(X), all(fit$diag.Psi > 0),
                all(vapply(outputs$invariants, function(z) all(is.finite(z)), logical(1))),
                all(outputs$entry >= 0 & outputs$entry <= 1),
                all(outputs$pleio >= 0 & outputs$pleio <= 1))
      TRUE
    })
    receipt$convergence <- if (is.null(fit$diagnostics)) {
      "Historical interface did not expose a convergence status; returned fit is not proof of convergence."
    } else fit$diagnostics
    if (!is.null(opt$compare)) {
      reference <- readRDS(opt$compare)
      comparison <- stage("baseline_comparison", {
        stopifnot(identical(reference$inputs, outputs$inputs))
        pieces <- lapply(names(outputs$invariants), function(name) {
          a <- outputs$invariants[[name]]; b <- reference$invariants[[name]]
          same_shape <- identical(dim(a), dim(b)) && length(a) == length(b)
          finite <- same_shape && all(is.finite(a)) && all(is.finite(b))
          delta <- if (finite) abs(as.numeric(a) - as.numeric(b)) else Inf
          absolute <- if (length(delta)) max(delta) else 0
          relative <- if (finite && length(delta)) max(delta / pmax(abs(as.numeric(b)), 1e-300)) else if (finite) 0 else Inf
          accepted <- finite && all(delta <= tolerance * (1 + abs(as.numeric(b))))
          data.frame(quantity = name, status = if (accepted) "PASS" else "FAIL",
                     max_abs = absolute, max_rel = relative, tolerance = tolerance)
        })
        tab <- do.call(rbind, pieces)
        write.table(tab, file.path(opt$output, "comparison.tsv"), sep = "\t",
                    row.names = FALSE, quote = TRUE)
        print(tab)
        if (any(tab$status != "PASS")) stop("Numerical compatibility failed; see comparison.tsv")
        tab
      })
      receipt$comparison <- comparison
      receipt$comparison_rule <- "abs(candidate-reference) <= tolerance * (1 + abs(reference)); factor products/covariance, not signs"
    }
    print(outputs$entry)
    print(outputs$pleio)
    receipt$status <- "PASS"
    0L
  }, error = function(e) {
    receipt$status <<- "FAIL"
    receipt$error <<- conditionMessage(e)
    cat("FAIL: ", conditionMessage(e), "\n", sep = "")
    1L
  })
  receipt$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  receipt$linux_process_memory <- if (file.exists("/proc/self/status")) {
    grep("^Vm(HWM|RSS):", readLines("/proc/self/status"), value = TRUE)
  } else "NOT RUN: Linux /proc memory counters unavailable"
  flush_receipt()
  status
}
quit(status = main(), save = "no")
