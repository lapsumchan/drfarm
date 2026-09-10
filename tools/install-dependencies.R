#!/usr/bin/env Rscript
# Version-pinned source installer. Tested local numerical dependencies were
# Ubuntu binaries; rebuilding these source versions is a separate CI check.
args <- commandArgs(trailingOnly = TRUE)
allowed <- args %in% c("--verify-only", "--plan") | startsWith(args, "--library=")
if (any(!allowed)) stop("Unknown arguments: ", paste(args[!allowed], collapse = ", "))
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "tools"
pins <- read.csv(file.path(script_dir, "dependencies.csv"), stringsAsFactors = FALSE)
if (anyDuplicated(pins$Package) || !identical(pins$Order, seq_len(nrow(pins)))) {
  stop("Invalid dependency order or duplicate package in dependencies.csv")
}
lib_arg <- sub("^--library=", "", args[startsWith(args, "--library=")])
if (length(lib_arg) > 1L) stop("Use at most one --library= argument")
user_lib <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1]]
lib <- if (length(lib_arg)) lib_arg else if (length(user_lib) && nzchar(user_lib[1]) && !grepl("%", user_lib[1], fixed = TRUE)) user_lib[1] else .libPaths()[1]
lib <- path.expand(lib)
if ("--plan" %in% args) {
  print(pins[, c("Order", "Package", "Version", "Role")], row.names = FALSE)
  cat("Destination:", lib, "\n")
  quit(status = 0L)
}
if (!dir.exists(lib) && !dir.create(lib, recursive = TRUE)) stop("Cannot create library: ", lib)
.libPaths(c(lib, .libPaths()))
installed_version <- function(p) {
  loc <- find.package(p, quiet = TRUE)
  if (!length(loc)) return(NA_character_)
  unname(read.dcf(file.path(loc, "DESCRIPTION"), fields = "Version")[1, 1])
}
receipt <- list()
for (i in seq_len(nrow(pins))) {
  p <- pins$Package[i]
  v <- pins$Version[i]
  found <- installed_version(p)
  if (identical(found, v)) {
    cat("MATCH", p, v, "\n")
    next
  }
  if ("--verify-only" %in% args) {
    stop(p, ": expected ", v, "; found ", if (is.na(found)) "absent" else found)
  }
  archive <- tempfile(paste0(p, "_"), fileext = ".tar.gz")
  urls <- c(pins$SourceURL[i], pins$ArchiveURL[i])
  downloaded <- FALSE
  errors <- character()
  for (url in urls) {
    status <- tryCatch(
      utils::download.file(url, archive, mode = "wb", quiet = FALSE),
      error = function(e) { errors <<- c(errors, conditionMessage(e)); 1L }
    )
    if (identical(status, 0L) && file.exists(archive) && file.info(archive)$size > 0) {
      downloaded <- TRUE
      break
    }
  }
  if (!downloaded) stop("Exact source unavailable for ", p, " ", v, ": ", paste(errors, collapse = "; "))
  unpack <- tempfile(paste0(p, "-metadata-"))
  dir.create(unpack)
  description <- paste0(p, "/DESCRIPTION")
  if (!description %in% utils::untar(archive, list = TRUE)) stop("Invalid source archive for ", p)
  utils::untar(archive, files = description, exdir = unpack)
  metadata <- read.dcf(file.path(unpack, description), fields = c("Package", "Version"))
  if (!identical(unname(metadata[1, ]), c(p, v))) stop("Downloaded source identity mismatch for ", p)
  utils::install.packages(archive, repos = NULL, type = "source", lib = lib, dependencies = FALSE)
  if (!identical(installed_version(p), v)) stop("Installation failed or wrong version for ", p, " ", v)
  receipt[[length(receipt) + 1L]] <- data.frame(Package = p, Version = v, URL = url,
                                             MD5 = unname(tools::md5sum(archive)))
  utils::write.csv(do.call(rbind, receipt), file.path(lib, "drfarm-source-install-receipt.csv"), row.names = FALSE)
  unlink(c(archive, unpack), recursive = TRUE)
}
cat("Verified", nrow(pins), "exact dependency versions under R", as.character(getRversion()), "\n")
if (getRversion() != "4.3.3") warning("The local reference environment used R 4.3.3; this R version requires its own validation.")
