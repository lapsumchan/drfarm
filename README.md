# drfarm

Debiased-regularized factor analysis regression for multiple continuous outcomes.

DrFARM estimates the associations of observed predictors while representing residual
outcome dependence with latent factors:

$$
Y = X\Theta^T + ZB^T + E.
$$

The sparse-group penalty encourages both predictor-level and individual
coefficient sparsity. The factors model the residual response component;
`Theta` is the coefficient matrix for the observed predictors. The package
contains the original continuous-response implementation in R and Rcpp.
Generalized-response work is separate and is not implemented by this API.

## Install

From R 4.3.0 or later, with a C++ compilation toolchain available:

```r
install.packages("remotes")
remotes::install_github("lapsumchan/drfarm")
```

The historical `0.1.0` source is pinned at commit
`be6d52ee796161e732f398da5eadfc3d40812f34`:

```r
remotes::install_github("lapsumchan/drfarm@be6d52ee796161e732f398da5eadfc3d40812f34")
```

This checkout is the `0.1.0.9000` maintenance candidate. To install it locally,
run `Rscript --vanilla tools/install-dependencies.R` for the recorded dependency
versions, then `R CMD INSTALL .` from the checkout. The dependency installer
includes the vignette tools. A development version
in this checkout does not imply that GitHub or a package registry has released it.

Compilation uses Rcpp. A missing compiler or dependency is an installation
failure, not evidence about the statistical method. On a build failure, retain
the full installation log and include `sessionInfo()` when reporting it. The CI
definition targets R 4.3.3 on Linux; a workflow file alone is not evidence of a
successful run or support on other environments.

## A small installation-to-result example

This example fits a deliberately small 2 by 2 remMap grid, then one DrFARM fit.
It demonstrates the interface; it does not perform a full tuning analysis.

```r
library(drfarm)
data("drfarm.dat", package = "drfarm")
X <- drfarm.dat$X
Y <- drfarm.dat$Y

initial <- remMap.whole(X, Y, n.lambda = 2)
precision <- precM(X)
fit <- DrFARM.one(
  X, Y, initial$Theta0, precision, k = 2,
  lambda1 = initial$lambda1.opt, lambda2 = initial$lambda2.opt,
  max.iter = 1000
)
fit$diagnostics
fit$Theta
```

The finite `max.iter` argument and diagnostics require this maintenance
candidate. A full script, including scale conversion and dimension checks, is
installed with the package:

```r
source(system.file("examples", "quickstart.R", package = "drfarm"))
```

Use `vignette("getting-started", package = "drfarm")` when vignettes were built.
Function reference pages are available through `help(package = "drfarm")`.

## Shapes, scale, and fitted components

| Object | Shape | Meaning |
|---|---|---|
| `X` | n by p | Observed predictors |
| `Y` | n by q | Continuous outcomes |
| `Theta`, `Theta0` | q by p | Outcome rows, predictor columns |
| `B` | q by k | Factor loadings |
| `E.Z` | n by k | Estimated latent factor scores |
| `precM(X)` | p by p | Predictor precision estimate |
| `drfarm.dat$Theta.t` | p by q | Bundled generating coefficients, transposed relative to the fitted API |

By default, fitting and precision estimation use `scale()` to center each column
and divide by its sample standard deviation. `Theta` therefore acts on
standardized predictors and yields the standardized observed-predictor
component. The package does not return an intercept or automatically transform
coefficients back. Use the same preprocessing when fitting, scoring, and doing
inference; `standardize = FALSE` requires the caller to supply the intended
scale and compatible initial coefficients and precision matrix.

For a default fit, recover raw-scale coefficients and their intercept as follows:

```r
sx <- apply(X, 2, sd)
sy <- apply(Y, 2, sd)
Theta.raw <- sweep(sweep(fit$Theta, 1, sy, "*"), 2, sx, "/")
intercept <- colMeans(Y) - drop(Theta.raw %*% colMeans(X))
Y.predictor <- sweep(X %*% t(Theta.raw), 2, intercept, "+")
```

This computes the observed-predictor component. For training participants with
`K = NULL`, `E.Z %*% t(B)` is the additional fitted latent component on the
standardized outcome scale. Those fitted scores use the observed outcomes;
they are not automatically available for a new participant. There is no
separate `predict()` method. Reject missing, nonfinite, or constant columns
before standardizing. Retain row and column identities with your input data.

The bundled example contains 500 participants, 10 predictors, and 5 outcomes.
Its package documentation describes a rank-2 simulated factor model; that is why
the example uses `k = 2`. The original simulation script and RNG seed are not
bundled. The stored data reproduce this demonstration, not the entire original
simulation experiment. Comparing `Theta` directly with `t(Theta.t)` also requires
matching the coefficient scale.

## Termination and inference

`DrFARM.one()` and `DrFARM.whole()` retain the historical default
`max.iter = Inf`; specify a finite budget in new workflows. Inspect diagnostics
before interpreting a fit. `converged = TRUE` requires the historical outer
loss-change criterion and the inner coefficient-change criterion. It does not certify parameter stability, a KKT condition,
or a global optimum. Loss increases and exhausted iteration budgets are
reported explicitly. The historical trial-return behavior on a loss increase is
preserved. `remMap.one(..., diagnostics = TRUE)` exposes the native iteration
report; its default return remains the coefficient matrix.

`entry.pvalue()` retains the original outer-debiasing implementation. Its
variance calculation requires positive residual variance, positive residual
degrees of freedom, and a compatible precision estimate. A successful software
check does not establish inferential calibration for a new dataset.
`pleio.pvalue()` retains the historical two-sided Cauchy transform
`2 * pcauchy(-abs(mean(1 / tan(pi * p))))` across each predictor's entry p-values.
This is not the usual one-sided ACAT tail; its intended test and boundary
behavior remain under review. It has not been silently replaced.

The optional kinship matrix `K` rotates participants into an eigenbasis in some
paths. A fixed-fit diagnostic reproduced a mismatch between the rotated evaluator
and the historical whole-grid score; the quickstart uses `K = NULL`. See [NEWS.md](NEWS.md) for compatibility
changes, [known issues](docs/KNOWN_ISSUES.md) for executed diagnostic scope,
and [CONTRIBUTING.md](CONTRIBUTING.md) for reproducible checks.

## Citation and authorship

Use `citation("drfarm")` for the package citation and version. Original authors:
Lap Sum Chan, Gen Li, and Peter X.K. Song. The historical README cites:

> Chan, Lap Sum, et al. “DrFARM: Identification and inference for pleiotropic gene
> in GWAS.” bioRxiv (2022).

Publication metadata beyond that historical record has not been reconciled in
this maintenance slice. DrFARM retains its GPL (>= 3) license and existing
third-party notices.
