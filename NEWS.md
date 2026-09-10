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
