# drfarm 0.1.0.9002

* Add a pkgdown review site, native-math README, curated reference and explicit
  method/status article. The getting-started vignette runs the original bundled
  workflow and preserves stopping warnings. Numerical R/Rcpp source is unchanged.
* Add pinned site-build dependencies, downloadable matching source/examples,
  and a Pages workflow whose deployment requires manual dispatch on main.


* Add the separately named `gaussian.ecm.reference()` optimization baseline for
  an explicit Gaussian observed-data likelihood with sparse-group penalties.
  It requires a supplied starting tuple and working scale, supports independent
  rows only, and uses the existing weighted coefficient solver.
* Use one coefficient tuple throughout each frozen-posterior coefficient,
  loading and variance cycle. Retain the full posterior residual-variance
  contribution and recompute the observed objective and stationarity at the
  complete candidate tuple.
* Check all three conditional-objective changes and the fresh observed
  objective. Reject failed complete trials, retain the last accepted tuple and
  report finite budgets, boundary failures and observed stationarity. An
  optional positive variance lower bound defines an explicit constrained target.
* Add an installed small example and `docs/GAUSSIAN_ECM_REFERENCE.md` with the
  mathematical target, controls, output meaning and comparison protocol.

Migration: this reference omits DrFARM's inner debiasing and therefore changes
the estimation procedure. It supplies no inference and inherits no DrFARM
inferential validation. The historical and weighted DrFARM paths, their default
behavior and their existing p-value functions remain unchanged. Gaussian ECM
stationarity is not a global-optimum guarantee or a generalized DrFARM extension.

# drfarm 0.1.0.9001

* Add `remMap.weighted()` for an explicitly specified, fixed-variance weighted
  sparse-group coefficient objective. It operates on the supplied working
  scale, takes variances in `sigma`, and does not automatically standardize or
  add an intercept. `C = 0` excludes an entry; `C = 2` exempts it from both the
  entry penalty and the predictor group norm.
* Use a monotone scalar root for unequal-curvature predictor blocks and report
  full-objective KKT residuals, the objective trace, and termination information.
  Finite coefficient/root budgets and unsuccessful solves remain explicit.
* Add opt-in `coefficient.update = "weighted"` with `weighted.control` to
  `DrFARM.one()` and `DrFARM.whole()`. The default `"historical"` retains the
  original coefficient update. A nonconverged weighted inner solve stops its
  outer caller; the standalone solver returns its last finite iterate with a
  warning and status.
* Add an installed analytic coefficient example and a developer derivation in
  `docs/WEIGHTED_UPDATE.md`. These describe a coefficient-subproblem correction;
  they do not establish convergence of the full DrFARM algorithm, repair
  optional-K scoring or the historical predictor Cauchy formula, or introduce
  generalized-response inference.

Migration: weighted mode changes the coefficient objective's numerical update
and can change estimates and tuning results. Select it explicitly, retain its
controls and statuses with results, and keep comparisons with the historical
mode identifiable. No historical numerical default is replaced.

# drfarm 0.1.0.9000

Development maintenance candidate based on public `0.1.0`, commit
`be6d52ee796161e732f398da5eadfc3d40812f34`.

* Add a short installation-to-result example, executable getting-started
  vignette, package citation, contribution guide, and CI definition.
* Clarify predictor/outcome coefficient orientation, default standardization,
  raw-scale conversion, and the distinction between the observed-predictor
  component and fitted latent scores.
* Repair the README code fence and source/help inconsistencies, including the
  `q` by `p` initial-coefficient shape and the `print.iter` argument.
* Add an appended `max.iter` argument to `DrFARM.one()` and `DrFARM.whole()`.
  Its default `Inf` preserves the historical unbounded budget. A finite budget
  can stop before the historical loss-change rule and is reported as such.
* Expose termination diagnostics and warnings for a loss increase or exhausted
  budget. Reported convergence requires the historical outer loss-change and inner
  coefficient-change rules only.
  Preserve the historical trial-return behavior on a loss increase; this change
  does not replace the objective or establish a global-optimization guarantee.
* Add optional `remMap.one(..., diagnostics = TRUE)` to return `Theta0` together
  with the native iteration diagnostic. The default still returns a matrix.
* Preserve the historical predictor-level two-sided Cauchy formula while its
  intended inferential interpretation is reviewed. Do not identify it as
  standard one-sided ACAT. The optional-`K` basis review also remains open.

Migration: existing calls retain their coefficient orientation and
standardization defaults. Code that requires an exact set of list names should
allow the added diagnostic fields. New workflows should specify a finite
`max.iter` and inspect the returned status. No generalized-response API is
introduced by this maintenance candidate.

A defined CI job is not a recorded GitHub Actions run. Local execution evidence
and known failures belong in the maintenance evidence record, with exact source
and dependency identities.

Known issues reproduced locally on the maintenance fixtures: all 25 bundled
historical-grid fits stop on loss increase; optional-K whole-grid scoring can
mix participant bases; and unequal-variance group shrinkage can meet its
coefficient-change rule without satisfying the natural weighted-objective KKT
condition. The numerical updates and historical selection remain unchanged.
See `docs/KNOWN_ISSUES.md` for conditions and reproducible commands.

The declared R minimum now matches the historical README requirement (4.3.0).
The executed reference in this slice is R 4.3.3 on Linux; other supported-range
versions and platforms require their own check receipts.
