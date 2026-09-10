# Contributing to drfarm

Keep changes small enough to explain and reproduce. Preserve original authors,
GPL notices, existing third-party attribution, and the public API unless a
behavior change is explicitly documented in `NEWS.md`.

## Local development

Work in an isolated checkout and R library. For example, from the repository:

```sh
mkdir -p .local-R-library
R_LIBS_USER="$PWD/.local-R-library" Rscript --vanilla tools/install-dependencies.R
R_LIBS_USER="$PWD/.local-R-library" R CMD INSTALL .
R_LIBS_USER="$PWD/.local-R-library" Rscript --vanilla inst/examples/quickstart.R
R_LIBS_USER="$PWD/.local-R-library" R CMD build .
R_LIBS_USER="$PWD/.local-R-library" R CMD check --no-manual drfarm_0.1.0.9000.tar.gz
```

The installer uses the recorded package versions in `tools/dependencies.csv`.
Retain its download/install log and record `sessionInfo()` alongside any claimed
result. Exact R package versions do not pin the compiler, BLAS/LAPACK, or system
libraries; the local execution evidence records those separately. A full manual additionally requires a
working LaTeX installation. The vignette requires knitr, rmarkdown, and Pandoc;
if any are unavailable, record the omitted checks explicitly.

Use `.local-R-library` only for temporary local dependencies. Do not add it to
Git. The package's executable behavior tests run during `R CMD check`. The CI
workflow runs installation, the quickstart, source build, and package checks;
it retains logs on failure. Check logs rather than inferring success from the
existence of a workflow or test file.

## Numerical changes

Read `AGENTS.md` and the compatible guidance in `docs/agent/` when those project
files are available. For each numerical rewrite, record the target, old/new
operations, justification, finite stopping/precision limits, and output meaning.
Choose the package and plumbing checks triggered by that change. A `PASS` needs
an executed command and its result; `NOT RUN` and failures remain visible.

Useful references test behavior independently: unequal predictor/outcome
shapes; a hand-derived scale transformation; a least-squares or sparse-group
subproblem with matched penalty normalization; factor products under rotations;
and an estimator together with its matching covariance. A numerical test does
not certify a statistical theorem. Preserve finite-budget and failure outcomes,
not only successful fits.

For performance changes, profile the complete path and major stages using
matched inputs, seed, tolerances, thread counts, dependency versions, and
cache state. Record wall/CPU time, memory, numerical agreement, and failures.
Measure a bottleneck before changing it. Preserve tuning order and tie-breaking;
keep expensive simulation studies outside routine CI.

## Reporting an issue or proposing a change

Include the source commit/package version, `sessionInfo()`, a small synthetic
reproducer, exact options, expected behavior, actual output, and complete error
or warning text. Avoid including participant-level or unpublished research data.
For generalized-response extensions, identify the supported family/link,
coefficient estimand, latent assumptions, objective, fit/predict behavior, and
inferential status separately from this continuous-response reference.

Prepare local commits and reviewable evidence before requesting release or
publication. Repository pushes, PR publication, tags, registries, and unpublished
research assets follow the project's actual authorization; this guide is not
permission to publish them.
