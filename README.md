# DrFARM

**Debiased-regularized factor analysis regression model**

Which predictors are associated with several continuous outcomes—and which
associations remain after accounting for shared residual variation?

DrFARM combines sparse multivariate regression with latent response factors.
It was developed for pleiotropy in multi-trait GWAS; observed predictors can
also be non-genetic. The implementation is **R with Rcpp**.

**Source versions:** `main` preserves the historical **0.1.0** package. The
quickstart below installs the checked **0.1.0.9002 development version** from a
pinned GitHub commit, including its explicitly labelled method options.

Start with this README, then use the development
[installation example](vignettes/getting-started.Rmd) and
[methods and status](vignettes/articles/methods-status.Rmd). Those same development
docs and R help generate the optional package website.

$$
\underbrace{Y}_{\text{outcomes}} =
\underbrace{X\Theta^\top}_{\text{observed predictors}} +
\underbrace{ZB^\top}_{\text{latent factors}} +
\underbrace{E}_{\text{noise}}.
$$

**Read it left to right:** observed predictors explain the mean; latent factors
explain residual dependence between outcomes. The coefficient matrix
$\Theta$ belongs to the observed predictors. It is **not** constrained to be
low rank.

| Data or parameter | Shape | One row represents |
|---|---|---|
| `X` / `Y` | n × p / n × q | A participant |
| `Theta` | q × p | An outcome; columns are predictors |
| `B` / `E.Z` | q × k / n × k | An outcome / a participant |

## Choose a method

| Path | Use it for | Status and interpretation |
|---|---|---|
| **Historical DrFARM** · `DrFARM.one()` / `DrFARM.whole()` | The original fitting and debiasing procedure | Default path; inspect outer and inner stopping diagnostics before inference |
| **Weighted coefficient option** · `coefficient.update = "weighted"` | A specified unequal-variance coefficient subproblem within DrFARM | Opt-in; inner debiasing remains; weighted inference is unvalidated |
| **Gaussian ECM reference** · `gaussian.ecm.reference()` | Optimization of an explicit Gaussian penalized likelihood, with independent rows | Different estimator; no inner debiasing and no inference supplied |

All three paths use continuous Gaussian-response machinery. The ECM reference
does not establish a generalized-response extension. There is no Python API.
See [methods and status](vignettes/articles/methods-status.Rmd) for the objectives,
assumptions and limits.

## Install the development version

This README documents **0.1.0.9002**. With R ≥ 4.3.0 and a C++ compiler,
install the checked development source from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "lapsumchan/drfarm@3ccf15176783ba7d0128c15403a2191e5491fd39",
  upgrade = "never"
)
```

The commit pins the source used for the local package checks and examples.
For the recorded dependency versions, follow the source-checkout instructions
in [Get started](vignettes/getting-started.Rmd). This development version has
no CRAN release or new Zenodo release DOI.

For the publicly available **historical 0.1.0** source only:

```r
install.packages("remotes")
remotes::install_github("lapsumchan/drfarm@be6d52ee796161e732f398da5eadfc3d40812f34")
```

That historical install does **not** contain `max.iter`, the new diagnostics,
the weighted option or `gaussian.ecm.reference()` shown in this candidate's help.
See [installation and troubleshooting](vignettes/getting-started.Rmd).

## Run the bundled example

After installing this candidate, this quickstart uses the original bundled
data: **500 participants, 10 predictors and 5 outcomes**. It fits a small remMap
initialization grid and one historical DrFARM model with two factors.

```r
library(drfarm)
data("drfarm.dat", package = "drfarm")
X <- drfarm.dat$X
Y <- drfarm.dat$Y

set.seed(20260909)
initial <- remMap.whole(X, Y, n.lambda = 2)
precision <- precM(X)
fit <- DrFARM.one(
  X, Y, initial$Theta0, precision, k = 2,
  lambda1 = initial$lambda1.opt, lambda2 = initial$lambda2.opt,
  standardize = TRUE, thres = 1e-4, max.iter = 1000
)
fit$diagnostics
dim(fit$Theta)  # 5 outcomes × 10 predictors
```

`Theta[r, j]` describes predictor j's association with outcome r on the
**standardized scale**. By default, X and Y are centered and divided by their
column sample standard deviations. Factor scores use the observed training
outcomes; they are not predictions for new participants.

The recorded quickstart returns `loss_increase`, not convergence. The example
is useful for learning the interface; its successful execution does not validate
the fit for inference. [Get started](vignettes/getting-started.Rmd) explains the
status, scale conversion, downloadable scripts and full original demonstration.

```r
source(system.file("examples", "quickstart.R", package = "drfarm"))
help(package = "drfarm")
?DrFARM.one
citation("drfarm")
```

## Interpretation and limits

Use matching preprocessing, coefficient scale and predictor precision when
fitting or evaluating inference. The original p-value routines are retained;
the historical predictor combination uses a **two-sided Cauchy tail**, and an
optional-kinship model-selection basis mismatch has been demonstrated. The
quickstart uses `K = NULL`. Lower objective values are not a global-optimum or
inferential guarantee. Read [methods and limitations](vignettes/articles/methods-status.Rmd)
before a scientific analysis.

## Paper, archived code and software citation

**Method paper:** Chan, L. S., Li, G., Fauman, E. B., Yin, X., Laakso, M.,
Boehnke, M. & Song, P. X. K. (2025). *DrFARM: identification of pleiotropic
genetic variants in genome-wide association studies.* Nature Communications,
16, 5789. [doi:10.1038/s41467-025-60439-4](https://doi.org/10.1038/s41467-025-60439-4).

**Historical code cited by the paper:** [DrFARM 0.1.0 on Zenodo](https://doi.org/10.5281/zenodo.15252156),
archived from this repository's `0.1.0` release. It does not contain the newer
development APIs documented here.

For a reproducible analysis, cite the method **and the exact software version**
used. A paper's historical code archive and this evolving development candidate
are different records. See the [citation and release policy](docs/CITING.md).

Package authors: **Lap Sum Chan, Gen Li and Peter X.K. Song**.
Use `citation("drfarm")` to cite the installed software and version. The
[citation record](inst/CITATION), [GPL ≥ 3 license and third-party notices](LICENSE.md),
[changelog](NEWS.md), and [contribution guide](CONTRIBUTING.md) are included.
Maintainer contact details are preserved in `DESCRIPTION`.

For the optional documentation build and proposed Pages publication path, see
[building the package website](docs/WEBSITE.md).

*Latent factors. Explicit assumptions.*
