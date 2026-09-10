#!/usr/bin/env Rscript
# Build only public files. Private project context is never an input to pkgdown.
started <- proc.time()[["elapsed"]]
if (length(commandArgs(trailingOnly = TRUE))) stop("Run without arguments from a source checkout")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
stopifnot(requireNamespace("pkgdown", quietly = TRUE))
if (as.character(packageVersion("pkgdown")) != "2.0.7") stop("Use tools/install-site-dependencies.R: this build pins pkgdown 2.0.7")
work <- file.path(root, "site-build")
dir.create(work, showWarnings = FALSE)
stage <- file.path(work, "source")
if (dir.exists(stage)) unlink(stage, recursive = TRUE)
dir.create(stage)
copy <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(from)) {
    dir.create(to, recursive = TRUE, showWarnings = FALSE)
    children <- list.files(from, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    if (length(children) && !all(file.copy(children, to, recursive = TRUE, overwrite = TRUE))) stop("Copy failed: ", from)
  } else if (!file.copy(from, to, overwrite = TRUE)) stop("Copy failed: ", from)
}
public <- c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "LICENSE.md",
            "CONTRIBUTING.md", ".Rbuildignore", "_pkgdown.yml", "R", "src", "man",
            "data", "inst", "tests", "vignettes", "pkgdown", "tools")
for (name in public) copy(file.path(root, name), file.path(stage, name))
unlink(list.files(file.path(stage, "src"), pattern = "\\.(o|so|dll)$", full.names = TRUE))
run <- function(args, log) {
  cat("R", paste(args, collapse = " "), "\n")
  status <- system2(file.path(R.home("bin"), "R"), shQuote(args), stdout = log, stderr = log)
  if (status != 0L) stop("Command failed; read ", log)
}
original_dir <- setwd(work)
run(c("CMD", "build", stage), file.path(work, "package-build.log"))
archive <- file.path(work, paste0("drfarm_", read.dcf(file.path(stage, "DESCRIPTION"))[1, "Version"], ".tar.gz"))
lib <- file.path(work, "library")
# This directory is owned by this builder; always install into an empty library.
if (dir.exists(lib)) unlink(lib, recursive = TRUE)
dir.create(lib)
run(c("CMD", "INSTALL", paste0("--library=", lib), archive), file.path(work, "package-install.log"))
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
writeLines(sort(getNamespaceExports("drfarm")), file.path(work, "installed-exports.txt"))
# Rewrite source-navigation links for the generated site, keeping GitHub's README native.
readme <- readLines(file.path(stage, "README.md"), warn = FALSE)
links <- c("vignettes/getting-started.Rmd" = "articles/getting-started.html",
           "vignettes/articles/methods-status.Rmd" = "articles/methods-status.html",
           "inst/CITATION" = "authors.html#citation", "LICENSE.md" = "LICENSE.html",
           "NEWS.md" = "news/index.html", "CONTRIBUTING.md" = "CONTRIBUTING.html",
           "docs/WEBSITE.md" = "WEBSITE.html", "docs/KNOWN_ISSUES.md" = "KNOWN_ISSUES.html")
for (old in names(links)) readme <- gsub(paste0("](", old, ")"), paste0("](", links[[old]], ")"), readme, fixed = TRUE)
writeLines(readme, file.path(stage, "README.md"))
assets <- file.path(stage, "pkgdown/assets")
copy(archive, file.path(assets, "downloads", basename(archive)))
for (name in list.files(file.path(root, "inst/examples")))
  copy(file.path(root, "inst/examples", name), file.path(assets, "downloads", name))
copy(file.path(root, "tools/reproduce.R"), file.path(assets, "downloads/reproduce.R"))
for (name in c("WEBSITE.md", "KNOWN_ISSUES.md", "GAUSSIAN_ECM_REFERENCE.md", "WEIGHTED_UPDATE.md", "OUTER_GAUSSIAN.md")) {
  copy(file.path(root, "docs", name), file.path(assets, "technical", name))
  # pkgdown renders these public technical notes with the same navigation/theme.
  # Add after R CMD build: they are website pages, not new package top-level files.
  copy(file.path(root, "docs", name), file.path(stage, name))
}
method_path <- file.path(stage, "vignettes/articles/methods-status.Rmd")
writeLines(gsub("../technical/", "../", readLines(method_path), fixed = TRUE), method_path)
capture.output(sessionInfo(), file = file.path(work, "sessionInfo.txt"))
# Optional cross-package link metadata must not cause network access while
# rendering. downlit catches unavailable metadata and uses its documented fallback.
# Keep this guard in this fresh build process only; package code is untouched.
options(pkgdown.internet = FALSE)
trace("download.file", where = asNamespace("utils"), print = FALSE,
      tracer = quote(stop("External metadata lookup disabled for local site rendering")))
pkgdown::build_site(stage, new_process = FALSE, install = FALSE, examples = TRUE, preview = FALSE)
untrace("download.file", where = asNamespace("utils"))
setwd(original_dir)
dest <- file.path(root, "_site")
if (dir.exists(dest)) unlink(dest, recursive = TRUE)
copy(file.path(stage, "_site"), dest)
cat("Site:", dest, "\nCandidate:", archive, "\nElapsed seconds:", proc.time()[["elapsed"]] - started, "\n")
