# Gaussian ECM optimization reference

Author: Student W. This is a separately named optimization baseline for a
specified Gaussian penalized likelihood. It is a changed estimation procedure;
it does not replace DrFARM or inherit DrFARM's inferential claims.

## Target and scope

`gaussian.ecm.reference()` fits the model

$$
y_i=\Theta x_i+Bz_i+\epsilon_i,\qquad
z_i\sim N_k(0,I),\qquad \epsilon_i\sim N_q(0,\Psi),
\qquad \Psi=\operatorname{diag}(\psi)>0.
$$

Rows, latent factors and errors are independent; the latent factors have mean
zero conditional on X. Theta is the q by p coefficient matrix for the observed
predictor mean. Conditional on X and Z, the mean is `X %*% t(Theta) + Z %*% t(B)`;
after integrating over Z it is `X %*% t(Theta)`. The latent structure belongs to
the response covariance, not to a low-rank factorization of Theta. These mean
associations have no automatic causal interpretation.

With `E = Y - X %*% t(Theta)` and `Sigma = B %*% t(B) + diag(psi)`, the objective
is the observed Gaussian negative log-likelihood plus a sparse-group penalty,
omitting the constant `n * q * log(2*pi) / 2`:

$$
L(\Theta,B,\psi)=\frac n2\log\det\Sigma+
\frac12\operatorname{tr}(\Sigma^{-1}E^TE)+\mathcal P_C(\Theta),
$$
$$
\mathcal P_C(\Theta)=\lambda_1\sum_{j}\sum_{r\in P_j}|\Theta_{rj}|
+\lambda_2\sum_j\|\Theta_{P_j,j}\|_2,
\qquad P_j=\{r:C_{rj}=1\}.
$$

The constraint is `Theta[C == 0] = 0`; C=2 entries are outside both penalties,
including the group norm. The likelihood and penalty are summed, with no
implicit division by n. Scaling the data, likelihood or penalties changes the
target unless the corresponding transformations are made explicitly.

Only independent rows (the `K = NULL` model) are implemented; there is no K
argument, and supplying one raises an unused-argument error. There is no automatic
centering, standardization, intercept, precision estimation, factor-rank
selection or tuning grid. The caller supplies the initial coefficient,
loading and variance tuple on the intended scale. An explicit constant design
column can represent an intercept; its mask determines its penalty. No
generalized-response family or p-value procedure is supplied.

Orthogonal rotations of B leave Sigma unchanged. Compare implied covariance,
observed-predictor means and fitted factor products, rather than requiring
arbitrary loading coordinates to agree. The problem is nonconvex in B and psi;
stationarity can occur at a saddle and does not establish global optimality or
unique factor identification.

## One frozen-posterior cycle

At the accepted tuple, compute

$$
S=(I+B^T\Psi^{-1}B)^{-1},\qquad
M=E\Psi^{-1}BS,\qquad T=M^TM+nS.
$$

M contains posterior mean factor scores. S is their common conditional
covariance; T includes both the mean cross-product and posterior uncertainty.
Freeze M, S and T for all three updates in this cycle. With these moments,
the relevant expected complete-data objective is

$$
Q(\Theta,B,\psi\mid old)=\frac n2\sum_r\log\psi_r+
\frac12\sum_r
\frac{\|(Y-X\Theta^T-MB^T)_{\cdot r}\|_2^2+n\,b_r^TSb_r}{\psi_r}
+\mathcal P_C(\Theta),
$$

where b_r is the transpose of loading row r. The omitted latent-prior term is
constant only while this posterior is frozen.

1. **Coefficients:** call the existing `remMap.weighted()` on
   `Y - M %*% t(B.old)`, with `sigma = psi.old`, C and the supplied penalties.
   Require its stated KKT stopping condition. Keep this new Theta for both
   following updates; there is no inner-debiasing operation.
2. **Loadings:** form `E.new = Y - X %*% t(Theta.new)` and solve
   `B.new = t(E.new) %*% M %*% solve(T)`.
3. **Variances:** use the full posterior expected residual sum of squares,
   `psi.new[r] = (sum((E.new - M %*% t(B.new))[, r]^2) +
   n * t(b.new[r]) %*% S %*% b.new[r]) / n`, subject to the explicit lower
   bound described below.

These are conditional minimizations of the same frozen Q. The loading update
does not depend on the positive diagonal variance weights because each response
row has a common weight. Its normal equations also justify the compressed
variance identity derived in [the outer reconciliation](OUTER_GAUSSIAN.md).
The reference evaluates the full residual-plus-posterior-variance expression so
the target of that calculation remains visible.

Q is evaluated before and after each conditional update. Only same-posterior Q
differences are meaningful here. Fresh posterior moments and the observed
objective L are then computed at the complete candidate tuple. The matching
Gaussian identity gives

$$
L(new)-L(old)=Q(new\mid old)-Q(old\mid old)
-\mathrm{KL}\{p_{old}(Z\mid Y)\Vert p_{new}(Z\mid Y)\}.
$$

Thus exact nonincrease of the frozen Q implies nonincrease of L. The finite
implementation checks the three Q changes and the freshly evaluated L directly;
it does not infer success from small coefficient movement or a stale posterior.

## API, scale and controls

```r
fit <- gaussian.ecm.reference(
  X, Y, Theta0, B0, psi0, lambda1, lambda2,
  C = NULL,
  control = list(max.iter = 500L, objective.tol = 1e-8,
                 score.tol = 1e-6, psi.min = 0),
  coefficient.control = list(tol = 1e-8, max.sweeps = 1000L,
                             root.tol = 1e-14, root.maxit = 200L)
)
```

| Quantity | Contract |
|---|---|
| X, Y | Finite numeric n by p and n by q matrices on the supplied scale, with n >= 2 and q >= 2; zero predictor columns are rejected |
| Theta0 | Finite q by p starting coefficients; nonzero C=0 entries are rejected |
| B0 | Finite q by k starting loadings; defines the fixed rank with `1 <= k < q` |
| psi0 | q strictly positive starting uniqueness variances, not standard deviations |
| C | q by p mask; default all 1 |
| `max.iter` | Positive integer outer trial budget; distinct from coefficient sweeps and scalar-root iterations |
| `objective.tol` | Nonnegative tolerance below 1 for the relative change in fresh observed penalized objective |
| `score.tol` | Positive threshold for the observed stationarity diagnostic described below |
| `psi.min` | Nonnegative scalar or q-vector of explicit variance lower bounds |
| `coefficient.control` | Controls of the existing weighted coefficient solver |

`psi.min = 0` means a strictly positive variance domain with no hidden floor.
A computed zero expected residual sum of squares at a zero bound reports
`variance_boundary`. Negative or nonfinite values are numerical failures, as is
underflow when a positive computed sum is divided by n. Underflow can also
affect the sum-of-squares calculation itself.
A positive bound changes
the optimization domain to `psi >= psi.min`: the conditional variance solution
is the larger of the full posterior expected residual variance and that bound.
Starting variances must belong to the requested domain. Report a chosen bound
with the result; it is part of the target, not merely a numerical tolerance.

## Monitoring, stationarity and failures

Every trial must keep its frozen-Q changes and fresh observed-objective change
nonpositive, allowing only `128 * .Machine$double.eps` times the largest of 1
and the compared objective magnitudes. A failed coefficient solve, variance
boundary, unevaluable trial or failed descent check rejects the **whole trial**.
The returned parameter tuple remains the last accepted tuple, with a warning.
Invalid inputs or unrepresentable initial calculations raise an error before
iteration. Rejected-trial
diagnostics are retained; the stopping flag never converts a rejected trial
into successful convergence.

For the fresh observed objective, write A for `solve(Sigma)` and
`G.Sigma = (n*A - A %*% t(E) %*% E %*% A) / 2`. Its smooth gradients are

$$
G_\Theta=-A E^TX,\qquad G_B=2G_\Sigma B,\qquad
g_\psi=\operatorname{diag}(G_\Sigma).
$$

The coefficient diagnostic is the maximum predictor-block distance of
`G.Theta` to the admissible negative sparse-group subgradient, including the
free-coordinate gradients. C=0 constraints are enforced separately. The loading diagnostic is
the Frobenius norm of G_B. At an active positive variance lower bound, the
variance diagnostic retains only the negative part of its gradient; interior
coordinates retain the full gradient. The aggregate score is the maximum of
these three residuals divided by n. **These are working-scale gradient units,
not a dimensionless error bound.** Record preprocessing and tolerances before
comparing scores across scales.

An initial tuple meeting the score criterion may return with zero cycles.
Otherwise successful stopping requires both the observed objective-change and
fresh observed stationarity criteria. The relative change is
`abs(L.new - L.old) / max(1, abs(L.old))`. A small loss change alone is insufficient.
Stationarity and finite nonincrease checks support an optimization comparison
at the recorded start; neither certifies a global minimum, estimator calibration
or a DrFARM inference theorem. Zero loadings can be stationary even when a
different loading start would reach a better objective.

The result has class `gaussian_ecm_reference` and returns Theta, B, `diag.Psi`,
fresh `E.Z` and `posterior.covariance`, plus diagnostics containing the method,
objective, stationarity, history containing the initial and accepted states,
and all attempted-trial statuses. `iterations` counts accepted cycles;
`attempted` also counts a rejected trial. `objective.change` is the last
accepted relative change, or zero before any accepted cycle. Its
`inference.status` is `"not_provided"`. Posterior factor scores use observed Y;
they are not new-participant predictors. The observed-predictor mean can be
computed as `X.new %*% t(fit$Theta)` only on a matching design and working scale.
The diagnostics also return the fresh unpenalized observed gradients and
`variance.bound.active`; `stationarity$raw` contains the three residual
magnitudes before division by n, rather than the gradient arrays.

## Comparison with preserved DrFARM

| Behavior | DrFARM, historical or weighted coefficient option | Gaussian ECM reference |
|---|---|---|
| Public entry point | `DrFARM.one()` / `DrFARM.whole()` | `gaussian.ecm.reference()` |
| Coefficient step | Historical native update by default; weighted update is opt-in | Weighted conditional minimization |
| Inner debiasing | Preserved, using the supplied predictor precision | Absent |
| Covariance tuple | Uses internally debiased coefficients, returns sparse coefficients | Uses and returns the same new coefficients |
| Outer monitor | Historical diagonal-noise loss | Fresh integrated Gaussian penalized likelihood |
| Loss-increase return | Historical trial-return behavior | Reject complete trial, return last accepted tuple |
| Initialization | Existing DrFARM initialization and preprocessing | Explicit caller-supplied tuple and working scale |
| Inference | Existing procedures and their stated limitations | None supplied or validated |

Migration is an explicit choice of estimator. Removing inner debiasing changes
the procedure, its finite estimates and any inferential argument attached to
those estimates. Do not substitute this result into `entry.pvalue()` or
`pleio.pvalue()` and describe the output as validated ECM inference. Both
preserved DrFARM paths remain available with their existing behavior.

The installed deterministic example starts from the asymmetric fixture used in
the earlier outer reconciliation:

```r
source(system.file("examples", "gaussian-ecm-reference.R", package = "drfarm"))
```

The installed development candidate's example ran in a clean R 4.3.3 process
on Linux. From that explicit start it accepted 65 cycles, returned
`stationarity_and_objective`, and reached L approximately 6.408187 with maximum
observed score/n approximately `9.54e-7`. The example independently recomputed
L from the returned covariance and residuals and checked agreement within
`1e-8`. This describes one starting state and the declared tolerances.

For an optimization comparison, use matching X, Y, working scale, C, penalties,
rank and explicit starting tuple; compute the same observed L for every returned
tuple. The older diagnostic's injected DrFARM initialization is a comparison
device, not a new public fitting argument. Preserve failures and iteration
budgets, and compare means/covariances rather than arbitrary factor signs.
The small fixtures establish their tested optimization behavior; they do not
estimate typical failure rates or comparative statistical performance.

The comparison used `max.iter=500`, observed `score.tol=1e-6`,
`objective.tol=1e-8`, `psi.min=0`, and coefficient tolerance `1e-10`.
The preserved procedures used their own monitored-loss tolerance `1e-8`;
its meaning differs from the ECM score condition.

| Fixed start | Initial L | Final DrFARM L (both coefficient options) | Final ECM L | ECM status |
|---|---:|---:|---:|---|
| Zero loadings, Gaussian stationary | 5.912023 | 6.500000 | 5.912023 | Initial stationary; 0 cycles |
| Nonzero loadings, Gaussian stationary | 9.197225 | 9.288467 | 9.197225 | Initial stationary; 0 cycles |
| Asymmetric nonstationary | 9.696191 | 7.347153 | 6.408187 | Stationarity and objective; 65 cycles |

The first two preserved fits report `loss_tolerance` after 3 and 55 attempts;
the third reports `loss_increase` after 2. Those are the preserved monitored
loss statuses, not Gaussian stationarity statements. They do not coincide with
the shorter one-cycle diagnostic of the earlier reconciliation.

A separate four-cell bundled-data comparison reuses one explicit initial tuple
for all methods, with n=500, p=10, q=5, k=2 and the recorded working scale.
The Gaussian objective is lower in each ECM cell than in the corresponding
preserved cells, but all four reach `max_iter=500`. Their maximum observed
scores/n are respectively `2.82e-5`, `1.02e-5`, `1.82e-5`, and `4.97e-6`.
The variance residual dominates each. No larger budget, favorable replacement
start, cross-method winner selection or inference was used for this result.

Reproduce the fixed comparison from the repository root after installing:

```sh
Rscript --vanilla tools/reconcile-outer.R --output outer-fixtures
Rscript --vanilla tools/compare-gaussian-ecm.R \
  --fixtures outer-fixtures/fixtures.rds --output ecm-comparison
```

Use fresh output directories to retain earlier receipts. The optional bundled
profile additionally accepts `--inputs` and `--preparation` paths to the saved
weighted-comparison assets; its receipt binds exact inputs and controls.
