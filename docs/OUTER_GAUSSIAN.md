# What the current outer DrFARM algorithm optimizes

Author: Student W. Source under study: `44e6a012e64d08040bd01276b79913183fef5299`
(development version 0.1.0.9001). This is a reconciliation of existing operations;
the historical and weighted fitting paths are preserved.

**The complete implemented cycle is a hybrid estimating procedure.** The weighted
coefficient step minimizes an explicit conditional Gaussian penalized objective.
The loading and variance steps minimize the Gaussian conditional objective at a
different, internally debiased coefficient matrix. The monitored loss omits the
factor contribution. These facts do not establish one Gaussian objective that
the full returned parameter tuple minimizes. The small fixtures below distinguish
these statements from optimizer or inferential guarantees.

## Question, estimand, and data

Write `beta = t(Theta)`, so X is n by p, Y is n by q, beta is p by q, and the
observed-predictor mean is `X beta`. Theta remains the q by p coefficient matrix
for observed predictors. B is q by k; the low-rank structure belongs to the latent
response component. All formulas below use the supplied working scale, after
the existing standardization and any K eigenbasis rotation. No intercept is
silently added. The examples turn standardization off.

The explicit Gaussian comparison model is

$$
y_i = \beta^T x_i + Bz_i + \epsilon_i,\qquad
z_i\sim N_k(0,d_iI_k),\qquad \epsilon_i\sim N_q(0,\Psi),
\quad \Psi=\operatorname{diag}(\psi_1,\ldots,\psi_q)>0.
$$

Latent factors and errors are independent, and rows are independent in the
working basis. For `K=NULL`, all $d_i$ are 1. For positive-definite K, let
$K=U\operatorname{diag}(d)U^T$; use $U^TX$ and $U^TY$. This corresponds to latent row covariance K in
the original participant basis and independent noise across participants.
Conditional on X, the factors have mean zero and are independent of X. Thus the
Gaussian model's marginal and factor-conditional mean coefficient is the same
Theta; neither this model statement nor source compatibility establishes a
causal interpretation or generalized-response inference.
Orthogonal changes $B\mapsto BR$ leave the implied covariance unchanged; factor
coordinates are not uniquely identified. Covariances and fitted factor products
are the appropriate comparisons. This invariance is not a claim that every
other parameter is identifiable on these small data fixtures.

For $P_j=\{r:C_{rj}=1\}$, define the existing weighted solver's penalty

$$
\mathcal P_C(\beta)=\lambda_1\sum_j\sum_{r\in P_j}|\beta_{jr}|
+\lambda_2\sum_j\|\beta_{j,P_j}\|_2,
\qquad \beta_{jr}=0\text{ when }C_{rj}=0.
$$

C=2 entries enter neither penalty and are outside the group norm. Losses are
summed; there is no division by n in the coefficient objective.

## Three distinct scalar objectives

With $e_i=y_i-\beta^T x_i$ and $V_i=\Psi+d_i BB^T$, the observed-data negative
log-likelihood plus penalty, omitting `n*q*log(2*pi)/2`, is

$$
L(\beta,B,\Psi)=\frac12\sum_i
\{\log\det V_i+e_i^TV_i^{-1}e_i\}+\mathcal P_C(\beta).
$$

At an old parameter tuple, let $m_i=E[z_i\mid y_i]$ and
$S_i=\operatorname{Var}(z_i\mid y_i)$, with the posterior frozen throughout the next conditional
updates. Direct Gaussian conditioning gives

$$
S_i=(d_i^{-1}I+B_{old}^T\Psi_{old}^{-1}B_{old})^{-1},\qquad
m_i=S_iB_{old}^T\Psi_{old}^{-1}(y_i-\beta_{old}^Tx_i).
$$

Put $T_i=S_i+m_i m_i^T$. The part of the expected complete-data negative
log-likelihood depending on the candidate parameters is

$$
Q(\beta,B,\Psi\mid old)=\frac n2\log\det\Psi+
\frac12\sum_i\left\{
(e_i-Bm_i)^T\Psi^{-1}(e_i-Bm_i)
+\operatorname{tr}(\Psi^{-1}BS_iB^T)\right\}
+\mathcal P_C(\beta).
$$

The expected latent-prior term is constant in these candidate parameters
because d is fixed; it can be omitted for same-posterior Q differences. Q values
computed under different old posteriors must not be compared as though their
omitted constants were identical.

For the matching Gaussian densities and the same penalty, the EM identity is

$$
L(\eta)-L(\eta_{old})=
Q(\eta\mid old)-Q(\eta_{old}\mid old)
-\mathrm{KL}\{p_{old}(Z\mid Y)\Vert p_{\eta}(Z\mid Y)\},
$$

where $\eta$ denotes $(\beta,B,\Psi)$. Nonincrease of one frozen Q therefore implies
nonincrease of L. This is the condition that must apply to the complete candidate
tuple; minimizing different blocks at different coefficient tuples does not
establish it.

The source monitors a third quantity:

$$
H(\beta,\Psi)=\frac12\sum_{i,r}\frac{(Y-X\beta)_{ir}^2}{\psi_r}
+\frac n2\sum_r\log\psi_r
+\lambda_1\sum_{j,r}|\beta_{jr}|
+\lambda_2\sum_j\|\beta_{j,\cdot}\|_2.
$$

H is a diagonal, no-factor Gaussian score with penalties on every entry. It
contains neither B nor posterior means/variances. It differs from L in general,
and from Q even if a posterior-mean residual is substituted without its variance
term. For mixed C, its penalties also differ from the coefficient solver's
penalties. Agreement in the special case B=0 and C=1 does not make the complete
cycle an optimizer of H or L.

## Reconcile each implemented operation

Let $M$ be the $p\times p$ supplied `precM`; let $\bar Z$ be the $n\times k$
matrix with rows $m_i^T$. The E-step in `DrFARM.one()` implements the posterior
above by a Woodbury formula. Its `E.zzt` is $T=\sum_i T_i$, not merely
$\bar Z^T\bar Z$. The subsequent source updates are

$$
\begin{aligned}
Y_{aug}&=Y-\bar ZB_{old}^T,\\
\beta_{db}&=\beta_s+MX^T(Y_{aug}-X\beta_s)/n,
\qquad E_{db}=Y-X\beta_{db},\\
B_{new}&=E_{db}^T\bar ZT^{-1},\\
\psi_{new}&=\operatorname{diag}(E_{db}^TE_{db}-B_{new}\bar Z^TE_{db})/n.
\end{aligned}
$$

| Source operation | Gaussian interpretation and limit |
|---|---|
| Weighted coefficients $\beta_s$ | Conditional minimization of Q over $\beta$ at fixed $B_{old},\Psi_{old}$; the posterior-variance trace is constant in $\beta$. The finite KKT result applies to this subproblem. |
| Historical coefficient update | Preserved native path. Unequal-variance group shrinkage need not minimize that weighted objective; the previous counterexample remains applicable. |
| Inner debiasing $\beta_{db}$ | A one-step score correction. It does not include sparse-group penalties or enforce C=0. If $M=(X^TX/n)^{-1}$, it gives the unpenalized least-squares coefficients for $Y_{aug}$. |
| Loading update $B_{new}$ | Exact loading conditional minimizer of Q at $\beta_{db}$, subject to invertibility of T. It is not computed at the sparse $\beta_s$ returned by the fit. |
| Variance update $\psi_{new}$ | Gaussian variance conditional minimizer at $\beta_{db}$ and $B_{new}$, by the loading normal equations; requires positive resulting variances. |
| H at $\beta_s,\psi_{new}$ | Stopping score, not the Gaussian objective used in those conditional updates. Its decrease is not a KKT or likelihood certificate. |

The compressed variance update **does retain posterior factor uncertainty**.
For outcome r, the loading normal equation gives

$$
\frac1n\left\{\|E_{db,\cdot r}-\bar Zb_r\|_2^2
+\sum_i b_r^TS_i b_r\right\}
=\frac1n\left\{\|E_{db,\cdot r}\|_2^2
-b_r^T\bar Z^TE_{db,\cdot r}\right\},
$$

where $b_r$ is row r of $B_{new}$ written as a column. The equality uses
$T b_r=\bar Z^T E_{db,\cdot r}$. Dropping the positive trace term and retaining only the
posterior-mean residual would be a different variance estimate. Conversely,
this correct algebra does not make the variance estimate appropriate at $\beta_s$:
the normal equations were formed using $\beta_{db}$.

For fixed frozen posterior, cyclically minimizing Q in $\beta,B,\Psi$ using the
**same candidate beta** is a coherent Gaussian conditional-maximization
construction. The current implementation includes the additional $\beta_{db}$ step
and returns $\beta_s$. This note does not replace that implementation with an ECM
algorithm or claim that every possible scalar objective has been ruled out.

## Exact failure fixture

Take $n=2,p=1,q=2,k=1$, $X=(1,-1)^T$, and

$$
\begin{aligned}
Y&=\begin{pmatrix}3&4\\-1&-2\end{pmatrix},\quad
\beta_0=(0,0),\quad B_0=0,\quad \Psi_0=\operatorname{diag}(5,10),\\
M&=1,\quad \lambda_1=1,\quad\lambda_2=0,\quad C=1.
\end{aligned}
$$

This is an explicitly injected valid Gaussian initial state, not a claim about
what `psych::fa()` estimates from these two rows. Both coefficient modes agree
here because the group penalty is zero. The posterior has $m_i=0$ and $S_i=1$.

The first weighted minimizer is $\beta_s=(0,0)$; inner debiasing gives $\beta_{db}=(2,3)$.
Then $E_{db}$ has both rows $(1,1)$, so $B_{new}=0$ and $\psi_{new}=(1,1)$. The returned sparse
coefficient remains zero. Consequently

$$
L(\beta_0,B_0,\Psi_0)=\log 50+2\approx5.912023,
\qquad L(\beta_s,B_{new},\Psi_{new})=15.
$$

The full hybrid step increases the explicit Gaussian objective even though its
weighted coefficient update is exact. The source initializes the monitored
previous loss to `1e300`, so it accepts this first trial without comparing it
with the true initial H.

With further steps, the map reaches $\beta_s=(1.5,2.5)$, $\beta_{db}=(2,3)$, $B=0$,
$\psi=(1,1)$, and $H=L=6.5$. Repeating the map leaves that tuple unchanged. Yet at
the returned tuple,

$$
\frac{\partial L}{\partial\psi_r}
=\frac{n}{2\psi_r}-\frac{\|Y_{\cdot r}-X\beta_{\cdot r}\|_2^2}{2\psi_r^2}
=1-\frac{2.5}{2}=-0.25,\quad r=1,2.
$$

Thus even a fixed point of this small hybrid map need not be a stationary point
of the explicit Gaussian penalized likelihood. This example neither quantifies
frequency on scientific datasets nor addresses inferential calibration.
Holding its beta and B fixed while changing both variances to 1.25 lowers L to
`6 + 2*log(1.25) = 6.446287103`. This is an explicit descent direction.

## Nonzero-factor check and the variance identity

A second, fully specified fixture demonstrates that the issue is not confined
to zero loadings. Let $x=(-1,-1,1,1)^T$, $z=(1,-1,1,-1)^T$, and
$u=(1,-1,-1,1)^T$. Set

$$
X=x,\quad
Y=x(7/4,11/4)+z(\sqrt{15}/4,\sqrt{15}/4)+u(1,-1)/\sqrt2,
$$

with $\beta_{old}=(1,2)$, $B_{old}=(1,1)^T$, $\psi_{old}=(1,1)$, $M=1$, all $d_i=1$,
$\lambda_1=1$, $\lambda_2=0$, and C=1. The three design vectors are pairwise orthogonal
with squared norm 4. The initial residual covariance is exactly
$\left(\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right)=B_{old}B_{old}^T+\Psi_{old}$.

The weighted coefficient minimizer stays at (1,2). A matching factor/variance
conditional update at that sparse coefficient leaves B and Psi unchanged. Inner
debiasing instead gives (5/4,9/4); the source factor/variance formulas then yield
$B_{new}=(7/8,7/8)^T$ and $\psi_{new}=(59/64,59/64)$.

| Candidate, same posterior | Q | Observed L | Monitored H |
|---|---:|---:|---:|
| Initial / unchanged sparse coefficient | 7.000000000 | 9.197224577 | 11.000000000 |
| Debiased coefficient, old B and Psi | 7.250000000 | 9.280557911 | 10.250000000 |
| Sparse coefficient, $B_{new},\Psi_{new}$ | 7.081397103 | 9.247385563 | 11.352583544 |

Here the isolated debiasing step increases both Q and L while decreasing H.
The complete hybrid step also increases Q and L. Constants omitted from Q and L
differ; compare changes within a column, not magnitudes across columns.

For $B_{new}$ and $\beta_{db}$, the expected residual variance is 59/64 in both outcomes.
The posterior-variance contribution alone is 49/192, which is strictly positive
and is already included in the compressed source expression. At the **returned
sparse beta** and the same $B_{new}$, the correct full expected residual variance is
65/64 instead. Substituting that sparse residual into the compressed expression
without re-solving B would give 9/8; the loading normal equation no longer
justifies that compression. These exact fractions distinguish a valid algebraic
compression from using it at a mismatched parameter tuple.

An additional asymmetric nonzero-factor fixture in the execution receipts shows
the complementary failure: a coherent same-beta factor/variance conditional
update lowers Q from 6.603439690 to 5.736290581 and L from 9.003513317 to
8.038573117, while H rises from 10.850844534 to 11.437449741. The monitor can
therefore reject genuine Gaussian-objective progress as well as obscure a
debiasing step that moves in the wrong direction.

## What a returned stopping flag means

`loss_tolerance` combines the historical loss-change criterion with the inner
solver's criterion. It does not check the Gaussian likelihood gradient or
factor/variance normal equations at the returned sparse coefficient. On
`loss_increase`, the source breaks before restoring saved parameters, so the
increasing trial is returned. This behavior is preserved.

Returned `E.Z` normally comes from the E-step preceding the returned parameter
updates. It need not be the posterior mean at the returned tuple. An exact
zero-loss-change tie is a special case: the code restores the previous parameter
tuple, for which that E-step is current. Freshly recomputing posterior moments
is necessary when evaluating a returned tuple's conditional moments; it does
not retroactively change the implemented optimization path.

## Execution evidence and reproduction

The final actual-production trace ran on 9 September 2026 with isolated R 4.3.3,
Ubuntu 24.04.3 LTS, drfarm 0.1.0.9001, and one computational thread. The three
fixtures are deterministic explicit matrices; there is no random data generation
or simulation selection. They exercise the unchanged fitting operations under
fixed injected initial states, not a new public initialization API. The observer
and uninstrumented fixed-initializer closure returned exactly matching output
and warnings in both coefficient modes on all three fixtures.

| Question | Executed result and scope |
|---|---|
| Gaussian calculations and preservation | **PASS: 78 diagnostic assertions** in both modes; absolute tolerance 1e-10, or 1e-9 for the K basis check. |
| Separate reference implementation | **PASS: 32 checks** at 1e-12 absolute tolerance; no package numerical functions called. |
| Complete outer cycle always decreases L | **FAIL:** both counterexamples increase L. Tests pass by detecting these failures. |
| Loss convergence implies Gaussian stationarity | **FAIL:** the three-step fixture gives variance derivative (-0.25,-0.25). |
| Package check and inference in this slice | **NOT RUN anew.** Prior 0.1.0.9001 package-check evidence is reused; no inference performed. |

The 78 new diagnostic assertions are distinct from the 78 package-test outcomes
reported in the preceding maintenance slice. The complete R process took
1.411 seconds, with peak child RSS 214472 KiB; the diagnostic body took about
0.18 seconds. These are resource receipts for this tiny execution, not a speed
comparison or scaling claim.

The exact recorded invocation, run from the repository root after activating
the isolated runtime, was:

```sh
Rscript --vanilla tools/reconcile-outer.R \
  --library ../weighted-receipts/candidate-final-build/drfarm.Rcheck \
  --output ../outer-receipts/production-trace-r3
```

Use a new output directory when rerunning; the script refuses to overwrite
previous receipts. On another installation, omit `--library` to use the normal
R library search path, or point it to the isolated library containing the
matching package version. `tools/gaussian-objective.R` must accompany
`tools/reconcile-outer.R`; they evaluate L, fixed-posterior Q, and H independently
of the package monitor and expose their values at intermediate steps.

Receipts, relative to the workspace root one level above the checkout:
`outer-receipts/production-trace-r3/` contains `console.log`, `stages.tsv`,
`receipt.rds`, and `receipt.txt`; resource measurements are in
`outer-receipts/actual-production-trace-r3.resources.json`. The separate reference,
`outer-receipts/independent-fixtures.R`, has adjacent checks/values CSV, log, and
RDS outputs. Keep source and receipts together. No inference or generalized-
response validation is claimed.

Applicable checks: DF-01/02/03/04/06; P01/P05/P06/P09; N01/N02/N03.
Original attribution, both fitting paths, and their prior evidence remain intact.
