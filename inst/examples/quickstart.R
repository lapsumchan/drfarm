# A bounded demonstration using the installed maintenance candidate.
# Run in a clean session: Rscript --vanilla inst/examples/quickstart.R
# This is interface evidence, not a simulation study or full tuning analysis.
library(drfarm)
data("drfarm.dat", package = "drfarm")
set.seed(20260909)

X <- drfarm.dat$X
Y <- drfarm.dat$Y
stopifnot(
  is.matrix(X), is.matrix(Y), nrow(X) == nrow(Y),
  all(is.finite(X)), all(is.finite(Y)),
  all(apply(X, 2, sd) > 0), all(apply(Y, 2, sd) > 0)
)

# Four remMap candidates, in the package's original grid order and EBIC rule.
initial <- remMap.whole(X, Y, n.lambda = 2, standardize = TRUE)
precision <- precM(X, standardize = TRUE)

# One pair, selected by remMap: this does not tune the DrFARM grid.
fit <- DrFARM.one(
  X, Y, initial$Theta0, precision, k = 2,
  lambda1 = initial$lambda1.opt, lambda2 = initial$lambda2.opt,
  standardize = TRUE, thres = 1e-4, max.iter = 1000
)
print(fit$diagnostics)
print(fit$Theta)
stopifnot(
  identical(dim(fit$Theta), c(ncol(Y), ncol(X))),
  identical(dim(fit$B), c(ncol(Y), 2L)),
  identical(dim(fit$E.Z), c(nrow(X), 2L))
)

# Convert only the observed-predictor component back to the input scale.
# The fitted factor scores depend on observed outcomes and are not predictions
# of new participants' latent factors.
sx <- apply(X, 2, sd)
sy <- apply(Y, 2, sd)
Theta.raw <- sweep(sweep(fit$Theta, 1, sy, "*"), 2, sx, "/")
intercept <- colMeans(Y) - drop(Theta.raw %*% colMeans(X))
Y.predictor <- sweep(X %*% t(Theta.raw), 2, intercept, "+")
Y.predictor.from.standardized <- sweep(
  sweep(scale(X) %*% t(fit$Theta), 2, sy, "*"),
  2, colMeans(Y), "+"
)
stopifnot(isTRUE(all.equal(
  unname(Y.predictor), unname(Y.predictor.from.standardized),
  tolerance = 1e-10, check.attributes = FALSE
)))
print(head(Y.predictor))

# Inspect fit$diagnostics before any scientific interpretation or inference.
# No p-values are emitted by this demonstration; historical inference functions
# and their current scope are documented in the function help and README.
sessionInfo()
