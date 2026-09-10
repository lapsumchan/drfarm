#!/usr/bin/env Rscript
# Reuse the existing exact-version installer, adding only the documentation pins.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
core <- read.csv(file.path(root, "tools/dependencies.csv"), stringsAsFactors = FALSE)
site <- read.csv(file.path(root, "tools/site-dependencies.csv"), stringsAsFactors = FALSE)
site$Role <- "website"
pins <- rbind(core, site[, names(core)])
pins$Order <- seq_len(nrow(pins))
stage <- tempfile("drfarm-site-pins-")
dir.create(stage)
write.csv(pins, file.path(stage, "dependencies.csv"), row.names = FALSE)
file.copy(file.path(root, "tools/install-dependencies.R"), stage)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(file.path(stage, "install-dependencies.R")),
    shQuote(commandArgs(trailingOnly = TRUE))))
unlink(stage, recursive = TRUE)
quit(status = status)
