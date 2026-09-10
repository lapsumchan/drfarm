# Weighted coefficient update: objective and implementation contract

Author: Student W. This developer note describes a separately named coefficient
solver. Original package authorship and historical numerical behavior remain
unchanged. It does not establish a new DrFARM inferential result.

## Target and API

`remMap.weighted(X, Y, lambda1, lambda2, sigma, C, Theta0, control)` targets the
convex, fixed-variance coefficient objective

$$
F(\Theta)=\sum_{r=1}^{q}\frac{\|Y_{\cdot r}-X\Theta_{r\cdot}^{T}\|_2^2}
 {2\sigma_r}
 +\lambda_1\sum_{(r,j):C_{rj}=1}|\Theta_{rj}|
 +\lambda_2\sum_{j=1}^{p}\|\Theta_{P_j,j}\|_2,
\qquad P_j=\{r:C_{rj}=1\},
$$

subject to `Theta[C == 0] = 0`. Entries with `C == 2` are unpenalized and are
**outside both penalties**, including the predictor group norm. Empty groups
contribute zero. This explicitly chosen objective is the contract for the new
solver; the unequal-variance historical update is not assumed to minimize it.

| Input/output | Contract |
|---|---|
| `X`, `Y` | Finite numeric matrices, n by p and n by q, with n >= 2 and positive finite squared X-column norms, on the supplied working scale |
| `sigma` | q positive response **variances**, not standard deviations; default ones |
| `lambda1`, `lambda2` | Nonnegative finite entry/group penalties; no division by n is implicit |
| `C` | q by p mask: 0 excluded, 1 penalized, 2 unpenalized; default all 1 |
| `Theta0` | Optional q by p starting coefficients on that same scale |
| `control` | `tol = 1e-8`, `max.sweeps = 1000L`, `root.tol = 1e-14`, `root.maxit = 200L` |
| `Theta0` in result | q by p coefficient estimate for this objective |
| `diagnostics` | Objective trace, full-objective KKT residual, convergence and termination information |

Zero predictor columns are rejected. Nonzero constant design columns are
allowed on the supplied working scale. There is no automatic centering, scaling,
or intercept. Preprocess explicitly
and retain the transformation; pass an explicit design column if an intercept
is intended. Scaling the loss while leaving the penalties fixed changes the
target. Finite inputs alone do not rule out ill-conditioning or overflow.

The historical counterexample and its exact source/command are retained in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md). The reusable asset here is its demonstrated
failure of unequal-curvature radial shrinkage; it is not a valid replacement
kernel. This note derives the new block calculation directly from the chosen
objective. Applicable checks: DF-01/02/06, P05/06/09, N01/02. Runtime results and
source bindings belong in the corresponding change record; a derivation or an
example file does not establish that the implementation passed.

## One predictor block

For predictor j, let `b` be its q response coefficients and form the partial
residual by adding that predictor's current contribution back to `Y - X Theta'`.
Write

$$
d_j=X_{\cdot j}^{T}X_{\cdot j},\qquad
h_r=d_j/\sigma_r,\qquad
s_r=X_{\cdot j}^{T}R_{\cdot r,-j}/\sigma_r.
$$

For `d_j > 0`, each `h_r > 0`. Ignoring a constant, the block objective is

$$
\tfrac12\sum_r h_r b_r^2-\sum_r s_r b_r
 +\lambda_1\sum_{r\in P_j}|b_r|
 +\lambda_2\|b_{P_j}\|_2.
$$

Excluded entries stay zero. Unpenalized entries have `b_r = s_r / h_r` and do
not affect the penalized group norm. On the penalized coordinates, put
`u = soft(s, lambda1)`, with signed soft thresholding applied componentwise.

* If `||u|| <= lambda2`, the penalized block is zero. This is precisely the
  zero-block subgradient condition: the distance from s to the entry-penalty
  box is at most the group-penalty radius.
* If `lambda2 == 0`, use `b = u / h` directly.
* Otherwise, the nonzero block satisfies

$$
b_r=\frac{u_r}{h_r+t},\qquad t=\frac{\lambda_2}{\|b\|_2}>0,
\qquad
\left\|\frac{t u}{h+t}\right\|_2=\lambda_2.
$$

The final scalar function is continuous and strictly increasing in t when u is
nonzero and every h is positive. It starts at zero and tends to `||u||`.
Consequently the nonzero branch has one positive finite root. A bracketed root
solve supplies the block update. The solution preserves the signs of u; the
first-order conditions give the expression above. This is a derivation for
a positive diagonal quadratic metric, not a general rule allowing arbitrary
proximal operators to be composed.

When every h is equal, the formula reduces to ordinary radial group shrinkage
of the soft-thresholded vector. Unequal response variances generally give
unequal h, so that shortcut no longer supplies the same block solution. In floating-point computation the implementation transforms the root to
`alpha = t / (H + t)`, where `H = max(h)`, and solves on `[0, 1]` to avoid an
unbounded bracket. `root.tol` is an absolute tolerance on alpha, not a
coefficient-error bound. `max.block.root.residual` records a normalized scalar
equation residual. Root error, unsupported curvature ratios, and an exhausted
root budget remain distinct from an exact algebraic identity.

## Full-objective stopping

Predictor blocks are revisited in cyclic order. Changing one block changes the
partial residuals for correlated predictors, so a solved scalar root for one
block does not certify the full coefficient matrix.

Let `G = -t(crossprod(X, Y - X %*% t(Theta))) / sigma`, with division by sigma
along the outcome rows. At a solution:

* excluded coefficients meet the zero constraint;
* unpenalized coordinates have G equal to zero;
* a nonzero penalized group satisfies the entrywise lasso subgradient condition
  after adding `lambda2 * b / ||b||` to its smooth gradient;
* a zero penalized group has `||soft(G, lambda1)|| <= lambda2`.

The implementation recomputes residuals from X, Y and the current coefficients
after each sweep, and checks these conditions over the full matrix.
`kkt.residual` is the maximum per-predictor Euclidean distance to the admissible
subgradient set, including the unpenalized-coordinate gradients.
`kkt.scaled` divides each predictor residual by its fixed scale and then takes
the maximum. That scale is the maximum of 1, the norm of the weighted zero-fit
`X'Y` score over included coordinates, `lambda1 * sqrt(number of penalized
entries)`, and `lambda2`; penalty terms are omitted for entirely unpenalized
blocks. `control$tol` is the threshold on this **scaled** residual. The returned
`objective.trace` contains the initial objective and evaluated sweep iterates. In exact arithmetic the stated convex objective makes KKT
conditions sufficient for a global coefficient minimum. A finite computed
residual is approximate evidence; it is not an exact certificate, a bound on
coefficient error without conditioning assumptions, or an outer DrFARM result.
Rank deficiency can make coefficients nonunique.

On an exhausted budget or a solve failure leaving an evaluable finite iterate,
the standalone solver returns that iterate, status, and a warning. If a trial
cannot be evaluated, it returns the previous fully evaluated sweep with
`termination = "numerical_failure"` and the error in `failure`. Invalid inputs
or unrepresentable initialization calculations raise an error before iteration. A
DrFARM caller using this option stops
when the weighted coefficient solver fails to converge. Retain those failures
in summaries. `max.sweeps` limits coefficient sweeps and `root.maxit` limits
individual scalar solves; `max.iter` in DrFARM is a separate outer budget.
`sweeps` counts attempted sweeps, including a partial or rolled-back final attempt.

## Analytic fixture

Take a single predictor `X = (1, -1)' / sqrt(2)`, response slopes `(8, 5.5)`,
`sigma = (1, 0.25)`, `lambda1 = 2`, and `lambda2 = 5`. The block curvature is
`h = (1, 4)` and `u = (6, 20)`. The root `t = 1` gives `b = (3, 4)`, whose
norm is 5. Its gradient including both penalties is

$$
(1,4)\mathbin{\odot}(3,4)-(8,22)+(2,2)+(3,4)=(0,0).
$$

The objective is `17 + 14 + 25 = 56`. Positive curvature makes this fixture's
coefficient minimizer unique. The installed example checks these coefficients
against their analytic value with absolute tolerance `1e-7`; it also computes
the objective and gradient independently of the solver diagnostics:

```r
source(system.file("examples", "weighted-update.R", package = "drfarm"))
```

This case tests unequal variance and both penalties together. It cannot by
itself validate correlated-predictor cycling, all mask patterns, zero penalties,
rank deficiency, convergence failures, or the complete DrFARM procedure.

## Opt-in DrFARM integration and limits

`DrFARM.one()` and `DrFARM.whole()` accept
`coefficient.update = "weighted"` and `weighted.control = list(...)`.
The default `coefficient.update = "historical"` retains the original numerical
update. Weighted mode passes the current augmented response and uniqueness
variances to the solver above, on the working scale already selected by the
outer fit. Reuse an initial estimate only with matching data order and scale.

The option changes the coefficient subproblem. It does not repair or replace
the outer monitored loss, inner debiasing, factor/variance updates, optional-K
basis handling, model-selection scoring, or historical Cauchy combination.
The inherited outer monitored loss also includes C=2 coefficients in its
penalties, although the weighted coefficient objective excludes them. Outer
convergence remains a separate question. In particular, a small weighted
KKT residual does not justify inference from the full returned fit. This is a
continuous-response coefficient option, not a generalized-response model or a
transfer of Gaussian debiasing theory to new families.
