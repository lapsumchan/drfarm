# Known numerical issues and compatibility scope

Reference source: `be6d52ee796161e732f398da5eadfc3d40812f34` (0.1.0).
Development candidate: 0.1.0.9001. Local diagnostic environment: R 4.3.3,
glasso 1.11, glmnet 4.1-8, psych 2.4.1 and Rcpp 1.0.12. These are finite
software/algorithm checks, not a new assessment of published inferential theory.

## A returned historical fit need not have met a convergence criterion

On the bundled 500-by-10 design and 500-by-5 response, all 25 historical
DrFARM grid candidates stopped on an increase in the monitored loss. The
selected row was 22; its attempted outer iteration count was 2. The candidate
reports that state explicitly. It preserves the historical trial return and
first-EBIC-minimum selection, including eligibility of fits whose stopping
criterion was not met. Historical and candidate coefficients, factor products,
implied covariance, both p-value arrays and selected penalties agreed exactly
in this recorded environment. Agreement does not establish a valid optimizer.

Reproduce with `tools/run-checks.sh --baseline-source PATH --output-dir NEW_PATH
--mode full` in the recorded dependency environment. Supply an untouched
checkout of the reference commit as PATH. Each program is externally bounded
by the runner's timeout; an optional finite `max.iter` bounds attempted outer
updates but not time spent inside an individual native update.

## Optional K: the score must use the same participant basis

For a single fitted model with K, `DrFARM.one` returns scores in the eigenbasis
of K. `DrFARM.EBIC(..., K=K)` rotates X and Y consistently. The historical
inline score in `DrFARM.whole` instead combines original-basis X/Y with those
rotated scores; its grid initialization is also a distinct path.

A fixed-fit diagnostic on the first 80 bundled participants produced:

| K | Correctly rotated score | Historical inline unrotated score | Fit termination |
|---|---:|---:|---|
| Identity | -254.7427 | 25325.82 | loss_increase |
| Toeplitz plus unequal positive diagonal | -134.0761 | 6173.401 | max_iter (1000) |

The evaluator agreed with an independently written rotated score. This
establishes a basis mismatch at these finite returned fits. It does not
establish convergence, a changed winning grid row, or valid K-adjusted
p-values. Eigenvectors of repeated eigenvalues may vary across environments.
Run `Rscript --vanilla tools/profile.R --library PATH --output NEW_PATH
--kinship yes` to retain the data construction, spectra, fitted objects and
score calculations. The quickstart uses K=NULL.

## Unequal uniquenesses: coefficient stability is not weighted stationarity

Consider one predictor with squared norm 1, ordinary least-squares coefficients
a=(3,4), variance vector sigma=(1,2), entry penalty 0, and group penalty 1.
The historical kernel, still the candidate default, returns b=(2.4,2.4). For the natural
weighted objective

`f(b) = sum((b-a)^2 / (2*sigma)) + sqrt(sum(b^2))`,

the gradient there is `(0.1071067812, -0.0928932188)`, not zero. The candidate
native coefficient-change diagnostic nonetheless reports convergence with
final delta zero. The radial shrinkage identity applies to equal response
curvature; it does not justify this unequal-curvature step. This is a
counterexample to that objective interpretation, not a claim about every
possible objective.

Run `Rscript --vanilla tools/diagnose-weighted-group.R PATH`, where PATH is an
installed baseline or candidate R library. The script verifies the historical
output and the nonzero analytic gradient. It is separate from the passing
contract tests because the default preserves the historical update.

The separately named `remMap.weighted()` solver in 0.1.0.9001 targets this
explicit objective. On the same fixture it returns approximately
b=(2.323269,2.527538), lowers the objective from 4.214113 to 4.204097, and
has independently computed gradient infinity norm 2.22e-16 in the recorded
R environment. It agrees with a separate scalar-root reference to 8.89e-16.
Run `Rscript --vanilla tools/compare-weighted-update.R --library PATH --output
NEW_PATH` to reproduce the counterexample and a matched small-grid comparison.
See [WEIGHTED_UPDATE.md](WEIGHTED_UPDATE.md) for the precise target and limits.

All 20 weighted inner coefficient solves in that matched comparison met the
scaled KKT threshold 1e-8. All 20 outer fits across both update modes and two
passes still stopped on `loss_increase`; both four-cell grids selected row 3.
This demonstrates a coefficient-subproblem correction, not full-fit convergence.
The inherited outer monitored loss also includes C=2 entries in its penalties,
whereas the new coefficient objective excludes them; it is not the same target.

## Historical predictor combination has two tails

For one entry probability p strictly between zero and one, the historical
combination equals `2*min(p,1-p)`. Thus probabilities close to one may produce
small combined values. This differs from the conventional upper-tail Cauchy
combination. The contract tests preserve and expose this behavior. Boundary
roundoff and statistical calibration require separate treatment; this patch
does not introduce a corrected inference option.
