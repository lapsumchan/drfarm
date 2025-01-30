remMap <- function(X.m, Y.m, lamL1, lamL2, phi0 = NULL, C.m = NULL, sigma = NULL) {
  # Basic dimensions
  n <- nrow(X.m)
  p <- ncol(X.m)
  q <- ncol(Y.m)

  # If no C.m given, default to all penalized
  if (is.null(C.m)) {
    C.m <- matrix(1L, nrow = p, ncol = q)
  }

  # If sigma is NULL, default to rep(1, q)
  if (is.null(sigma)) {
    sigma <- rep(1, q)
  } else {
    # optional check: length(sigma) == q
    if (length(sigma) != q) {
      stop("sigma must be a numeric vector of length equal to the number of responses (q).")
    }
  }

  lambda1 <- lamL2
  lambda2 <- lamL1

  out <- MultiRegGroupLasso_sigma_unified(
    X_m = X.m,
    Y_m = Y.m,
    C_m = C.m,
    lam1 = lambda1,
    lam2 = lambda2,
    sigma = sigma,
    Phi_initial = phi0  # can be NULL or a numeric matrix
  )

  phi.result <- out$Phi_output
  E.matrix   <- out$E_debug
  rss.v      <- apply(E.matrix^2, 2, sum)

  list(phi = phi.result, rss.v = rss.v)
}

logsumchoose <- function(q, vec) {
  logchoose <- lgamma(q + 1) - lgamma(q - vec + 1) - lgamma(vec + 1)
  maxlogchoose <- max(logchoose)
  est <- maxlogchoose + log(sum(exp(logchoose - maxlogchoose)))
  return(est)
}

SoftThreshold <- function(x, lambda) {
  if (x > lambda) {
    return (x - lambda)
  } else {
    if (x < (-lambda)) {
      return (x + lambda)
    } else {
      return(0)
    }
  }
}

InverseLinftyOneRow <- function(sigma, i, mu,
                                maxiter = 50, threshold = 1e-2) {
  p <- nrow(sigma)
  rho <- max(abs(sigma[i,-i])) / sigma[i,i]
  mu0 <- rho / (1 + rho)
  beta <- rep(0,p)

  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i]
    returnlist <- list("optsol" = beta, "iter" = 0)
    return(returnlist)
  }

  diff.norm2 <- 1
  last.norm2 <- 1
  iter <- 1
  iter.old <- 1
  beta[i] <- (1 - mu0) / sigma[i,i]
  beta.old <- beta
  sigma.tilde <- sigma
  diag(sigma.tilde) <- 0
  vs <- -sigma.tilde %*% beta

  while ((iter <= maxiter) && (diff.norm2 >= threshold * last.norm2)){

    for (j in 1:p){
      oldval <- beta[j]
      v <- vs[j]
      if (j==i) v <- v + 1
      beta[j] <- SoftThreshold(v, mu)/sigma[j,j]
      if (oldval != beta[j]){
        vs <- vs + (oldval - beta[j]) * sigma.tilde[,j]
      }
    }

    iter <- iter + 1
    if (iter == 2 * iter.old){
      d <- beta - beta.old
      diff.norm2 <- sqrt(sum(d*d))
      last.norm2 <- sqrt(sum(beta*beta))
      iter.old <- iter
      beta.old <- beta
      if (iter > 10) vs <- -sigma.tilde %*% beta
    }
  }

  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol = 1.5,
                          mu = NULL, maxiter = 50, threshold = 1e-2,
                          verbose = TRUE) {
  isgiven <- 1
  if (is.null(mu)){
    isgiven <- 0
  }

  p <- nrow(sigma)
  M <- matrix(0, p, p)
  xperc <- 0
  xp <- round(p/10)
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc <- xperc + 10
      if (verbose) {
        print(paste(xperc,"% done",sep=""))
      }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1 - (0.1/(p^2)))
    }
    mu.stop <- 0
    try.no <- 1
    incr <- 0
    beta <- rep(0, p)   # initialize beta (missing in the original code snippet)
    while ((mu.stop != 1) && (try.no < 10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      } else {
        if (try.no == 1){
          if (iter == (maxiter+1)){
            incr <- 1
            mu <- mu * resol
          } else {
            incr <- 0
            mu <- mu / resol
          }
        }
        if (try.no > 1){
          if ((incr == 1) && (iter == (maxiter+1))){
            mu <- mu * resol
          }
          if ((incr == 1) && (iter < (maxiter+1))){
            mu.stop <- 1
          }
          if ((incr == 0) && (iter < (maxiter+1))){
            mu <- mu / resol
          }
          if ((incr == 0) && (iter == (maxiter+1))){
            mu <- mu * resol
            beta <- last.beta
            mu.stop <- 1
          }
        }
      }
      try.no <- try.no + 1
    }
    M[i,] <- beta
  }
  return(M)
}

#' Generate tuning parameter grid for \code{remMap}
#'
#' Generates a 2D grid of tuning parameters (\code{lambda1}, \code{lambda2}) for \code{remMap}
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param n.lambda The number of lambda points in a single dimension. Default is \code{10}
#' @param lambda.min.ratio The smallest value for lambda as a fraction of \code{lambda.max}. Default is \code{0.01}
#'
#' @return A data frame with two columns: \code{lambda1} and \code{lambda2}
#'
#' @export
remMap.grid <- function(X, Y, standardize = TRUE,
                        n.lambda = 10, lambda.min.ratio = 0.01) {

  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  lambda.vec <- rep(NA, q)
  for (j in 1:q) {
    lambda.vec[j] <- max(abs(colSums(X * Y[,j])))
  }
  lambda.max <- max(lambda.vec)

  lambda.seq <- round(exp(seq(log(lambda.max), log(lambda.max*lambda.min.ratio),
                              length.out = n.lambda)), digits = 10)
  remMap.lambda.grid <- expand.grid(lambda.seq, lambda.seq)
  colnames(remMap.lambda.grid) <- c("lambda1", "lambda2")

  return(remMap.lambda.grid)
}

#' Fit one \code{remMap} model
#'
#' Fit \code{remMap} for a single pair (\code{lambda1}, \code{lambda2}), serving as a helper function for
#' \code{remMap.all} to facilitate for parallelization.
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param lambda1 A user supplied value for the lasso tuning parameter
#' @param lambda2 A user supplied value for the group lasso tuning parameter
#' @param C A size \eqn{q \times p} integer matrix specifying model inclusion and penalty rules:
#'   \itemize{
#'     \item 0 = excluded from model
#'     \item 1 = penalized
#'     \item 2 = unpenalized
#'   }
#'   If \code{null}, all entries are penalized (1) by default
#'
#' @return A \eqn{q \times p} coefficient matrix, \code{Theta0}
#'
#' @export
remMap.one <- function(X, Y, standardize = TRUE,
                       lambda1, lambda2, C = NULL) {

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  if (is.null(C)) {
    Theta0 <- t(remMap(X, Y, lambda1, lambda2)$phi)
  } else {
    Theta0 <- t(remMap(X, Y, lambda1, lambda2, C.m = t(C))$phi)
  }

  return(Theta0)
}

#' Compute EBIC for \code{remMap}
#'
#' Computes the Extended Bayesian Information Criterion (EBIC) for a given \eqn{\Theta0}.
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta0 A \eqn{q \times q} coefficient matrix
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{1}
#'
#' @return A numeric EBIC value
#'
#' @export
remMap.EBIC <- function(X, Y, Theta0, standardize = TRUE,
                        gamma = 1) {

  n <- dim(Y)[1]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  E <- Y - X %*% t(Theta0)
  RSS <- diag(t(E) %*% E)
  EBIC <- n * sum(log(RSS)) +
    log(n) * sum(Theta0 != 0) +
    2 * gamma * logsumchoose(q, rowSums(t(Theta0) != 0))

  return(EBIC)
}

#' Fit all \code{remMap} models at once
#'
#' Searches over a 2D grid of (\code{lambda1}, \code{lambda2}) pairs, computes \code{EBIC} for each, and
#' returns the \eqn{q x p} coefficient matrix that minimizes the EBIC
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param n.lambda The number of lambda points in a single dimension. Default is \code{10}
#' @param lambda.min.ratio The smallest value for lambda as a fraction of \code{lambda.max}. Default is \code{0.01}
#' @param C A size \eqn{q \times p} integer matrix specifying model inclusion and penalty rules:
#'   \itemize{
#'     \item 0 = excluded from model
#'     \item 1 = penalized
#'     \item 2 = unpenalized
#'   }
#'   If \code{null}, all entries are penalized (1) by default
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{1}
#'
#' @return A list containing:
#' \item{Theta0}{The EBIC-chosen \eqn{q \times p} coefficient matrix}
#' \item{lambda1.opt}{Optimal \eqn{\lambda_1} value}
#' \item{lambda2.opt}{Optimal \eqn{\lambda_2} value}
#'
#' @export
remMap.whole <- function(X, Y, standardize = TRUE,
                         n.lambda = 10, lambda.min.ratio = 0.01,
                         C = NULL, gamma = 1) {

  n <- dim(Y)[1]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  lambda.vec <- rep(NA, q)
  for (j in 1:q) {
    lambda.vec[j] <- max(abs(colSums(X * Y[,j])))
  }
  lambda.max <- max(lambda.vec)

  lambda.seq <- round(exp(seq(log(lambda.max), log(lambda.max*lambda.min.ratio),
                              length.out = n.lambda)), digits = 10)
  remMap.lambda.grid <- expand.grid(lambda.seq, lambda.seq)
  colnames(remMap.lambda.grid) <- c("lambda1", "lambda2")

  n.lambda.sq <- dim(remMap.lambda.grid)[1]
  ls <- vector("list", n.lambda.sq)

  for (i in 1:n.lambda.sq) {
    if (is.null(C)) {
      ls[[i]] <- t(remMap(X, Y, remMap.lambda.grid[i,1], remMap.lambda.grid[i,2])$phi)
    } else {
      ls[[i]] <- t(remMap(X, Y, remMap.lambda.grid[i,1], remMap.lambda.grid[i,2], C.m = t(C))$phi)
    }
  }

  EBIC.vec <- rep(NA, n.lambda.sq)
  for (i in 1:n.lambda.sq) {
    E <- Y - X %*% t(ls[[i]])
    RSS <- diag(t(E) %*% E)
    EBIC.vec[i] <- n * sum(log(RSS)) +
      log(n) * sum(ls[[i]] != 0) +
      2 * gamma * logsumchoose(q, rowSums(t(ls[[i]]) != 0))
  }

  opt.idx <- which.min(EBIC.vec)

  return(list(Theta0 = ls[[opt.idx]],
              lambda1.opt = remMap.lambda.grid[opt.idx,1],
              lambda2.opt = remMap.lambda.grid[opt.idx,2]))
}

#' Estimate precision matrix via \code{glasso}
#'
#' Searches over a grid of tuning parameters in \code{rholist}, computes \code{EBIC} for each, and
#' returns the \eqn{p x p} precision matrix that minimizes the EBIC
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} by column. Default is \code{TRUE}
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{0.5}
#' @param rholist A list of nonnegative regularization parameters for the lasso
#'   (penalizing off-diagonal entries in the \code{glasso} framework).
#' @param thr,maxit Numerics specifying convergence threshold and maximum iterations for \code{glasso},
#'   respectively. Defaults are \code{1e-4} and \code{1e4}.
#'
#' @return A size \eqn{p \times p} precision matrix estimate
#'
#' @export
precM.glasso <- function(X, standardize = TRUE,
                         gamma = 0.5, rholist = NULL,
                         thr = 1e-4, maxit = 1e4) {

  n <- dim(X)[1]
  p <- dim(X)[2]

  if (standardize == TRUE) {
    X <- scale(X)
  }

  CovM <- t(X) %*% (X) / n
  a <- glassopath(CovM, rholist, thr, maxit)
  n.grid <- dim(a$wi)[3]
  EBIC.vec <- rep(NA, n.grid)

  for (i in 1:n.grid) {
    L <- n * (tr(CovM %*% a$wi[,,i]) - determinant(a$wi[,,i])$modulus[1]) / 2
    E <- sum(a$wi[,,i] != 0)
    EBIC.vec[i] <- 2 * L + E * log(n) + 4 * gamma * E * log(p)
  }

  opt.idx <- which.min(EBIC.vec)
  return(a$wi[,,opt.idx])
}

#' Estimate one row of precision matrix via Nodewise Lasso
#'
#' Estimates one row of the precision matrix via Nodewise Lasso, serving as a helper function
#' for \code{precM.NL.all} to facilitate parallelization.
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param row.idx An integer from 1 to \eqn{p} indicating which row of the precision matrix to estimate.
#' @param method A character, either \code{"EBIC"} or \code{"CV"}, determining how \code{lambda} is chosen
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{0.5}
#' @param n.lambda An integer giving the number of \code{lambda} points for \code{glmnet}. Default is \code{100}
#' @param lambda A nonnegative numeric vector (optional user-supplied \code{lambda} grid). Default is \code{NULL}
#'
#' @return A list with the following elements:
#' \item{beta}{A length-\eqn{(p-1)} vector of coefficients for the selected row of the precision matrix}
#' \item{tau.sq}{A numeric value representing the residual variance plus penalty contribution}
#'
#' @export
precM.NL.one <- function(X, row.idx, method = "EBIC", gamma = 0.5,
                         n.lambda = 100, lambda = NULL) {

  n <- dim(X)[1]
  p <- dim(X)[2]

  if (method == "CV") {
    cv.fit <- cv.glmnet(x = X[,-row.idx], y = X[,row.idx],
                        family = "gaussian", nlambda = n.lambda,
                        lambda = lambda, standardize = FALSE,
                        intercept = FALSE)
    lambda.opt <- cv.fit$lambda[cv.fit$index[2]]
    beta <- cv.fit$glmnet.fit$beta[, cv.fit$index[2]]
    SSE <- colSums((X[,row.idx, drop=FALSE] - X[,-row.idx] %*% beta)^2)
    tau.sq <- SSE / n + lambda.opt * sum(abs(beta))
  } else {
    fit <- glmnet(x = X[,-row.idx], y = X[,row.idx],
                  family = "gaussian", nlambda = n.lambda,
                  lambda = lambda, standardize = FALSE,
                  intercept = FALSE)
    lambda <- fit$lambda
    n.grid <- length(lambda)
    SSE.vec <- colSums((matrix(X[,row.idx], n, n.grid)
                        - X[,-row.idx] %*% fit$beta)^2)
    df.vec <- colSums(fit$beta != 0)
    EBIC.vec <- SSE.vec + df.vec * log(n) + 2 * gamma * lchoose(p, df.vec)
    opt.idx <- which.min(EBIC.vec)
    beta <- fit$beta[,opt.idx]
    tau.sq <- SSE.vec[opt.idx] / n + lambda[opt.idx] * sum(abs(fit$beta[,opt.idx]))
  }
  ls <- list(beta = beta, tau.sq = tau.sq)
  return(ls)
}

#' Estimate whole precision matrix via Nodewise Lasso
#'
#' Estimates the entire precision matrix via Nodewise Lasso by minimizing \code{EBIC}
#' or via cross-validation (\code{CV}).
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} by column. Default is \code{TRUE}
#' @param method A character, either \code{"EBIC"} or \code{"CV"}, determining how \code{lambda} is chosen
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{0.5}
#' @param n.lambda An integer giving the number of \code{lambda} points for \code{glmnet}. Default is \code{100}
#'
#' @return A \eqn{p \times p} precision matrix estimate
#'
#' @export
precM.NL.whole <- function(X, standardize = TRUE, method = "EBIC",
                           gamma = 0.5, n.lambda = 100) {

  p <- dim(X)[2]

  if (standardize == TRUE) {
    X <- scale(X)
  }

  C.hat <- diag(1, p, p)
  tau.sq.vec <- rep(NA, p)

  for (i in 1:p) {
    ls <- precM.NL.one(X, i, method, gamma, n.lambda)
    C.hat[i,-i] <- -ls$beta
    tau.sq.vec[i] <- ls$tau.sq
  }

  precM <- diag(1/tau.sq.vec) %*% C.hat
  return(precM)
}

#' Estimate precision matrix via Quadratic Optimization
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param standardize A logical indicating whether to standardize \code{X} by column. Default is \code{TRUE}
#'
#' @return A size \eqn{p \times p} precision matrix
#'
#' @export
precM.QO <- function(X, standardize = TRUE) {

  n <- dim(X)[1]

  if (standardize == TRUE) {
    X <- scale(X)
  }

  CovM <- t(X) %*% (X) / n
  precM <- InverseLinfty(CovM, n, resol = 1.3)
  return(precM)
}

#' Precision Matrix Estimation Wrapper
#'
#' @param X Data matrix \eqn{n \times p}.
#' @param method One of \code{"glasso"}, \code{"NL"}, or \code{"QO"}.
#' @param standardize Whether to standardize \code{X}.
#'
#' @return A \eqn{p \times p} precision matrix estimate.
#'
#' @export
precM <- function(X, method = "glasso", standardize = TRUE) {

  if (method == "NL") {
    pm <- precM.NL.whole(X, standardize)
  } else if (method == "QO") {
    pm <- precM.QO(X, standardize)
  } else {
    pm <- precM.glasso(X, standardize)
  }
  return(pm)
}

#' Generate tuning parameter grid for \code{DrFARM}
#'
#' Generates a 2D grid of tuning parameters (\code{lambda1}, \code{lambda2}) for \code{DrFARM}
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta0 A size \eqn{q \times p} initial coefficient matrix
#' @param precM A size \eqn{p \times p} precision matrix
#' @param k A positive integer specifying the number of latent factors used in DrFARM
#' @param lambda1.opt The chosen lasso tuning parameter from \code{remMap.whole}
#' @param lambda2.opt The chosen group-lasso tuning parameter from \code{remMap.whole}
#' @param K An optional size \eqn{n \times n} kinship matrix. Default is \code{NULL}
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#'
#' @return A data frame with two columns: \code{lambda1} and \code{lambda2}
#'
#' @export
DrFARM.grid <- function(X, Y, Theta0, precM, k,
                        lambda1.opt, lambda2.opt,
                        K = NULL,
                        standardize = TRUE) {

  n <- dim(X)[1]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  if (!is.null(K)) {
    eigen.res <- eigen(K)
    U <- eigen.res$vectors
    d <- eigen.res$values
    X <- t(U) %*% X
    Y <- t(U) %*% Y
  } else {
    d <- rep(1, n)
  }

  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  E.star <- Y - X %*% Theta0.db.t
  fa.res <- fa(E.star, nfactors = k, rotate = rotate, scores = scores, fm = fm, covar = TRUE)
  diag.Psi <- fa.res$uniquenesses

  Psi.range <- unname(summary(diag.Psi)[-4])
  DrFARM.lambda.grid <- expand.grid(lambda1.opt / Psi.range, lambda2.opt / Psi.range)
  colnames(DrFARM.lambda.grid) <- c("lambda1", "lambda2")

  return(DrFARM.lambda.grid)
}

#' Fit one \code{DrFARM} model
#'
#' Fit \code{DrFARM} for a single pair (\code{lambda1}, \code{lambda2}), serving as a helper function for
#' \code{DrFARM.all} to facilitate for parallelization.
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta0 A size \eqn{q \times p} initial coefficient matrix
#' @param precM A size \eqn{p \times p} precision matrix
#' @param k A positive integer specifying the number of latent factors used in DrFARM
#' @param lambda1 A user supplied value for the lasso tuning parameter
#' @param lambda2 A user supplied value for the group lasso tuning parameter
#' @param K An optional size \eqn{n \times n} kinship matrix. Default is \code{NULL}
#' @param C A size \eqn{q \times p} integer matrix specifying model inclusion and penalty rules:
#'   \itemize{
#'     \item 0 = excluded from model
#'     \item 1 = penalized
#'     \item 2 = unpenalized
#'   }
#'   If \code{null}, all entries are penalized (1) by default
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param thres A numeric value for the convergence threshold. Default is \eqn{1 \times 10^{-4}}
#' @param rotate Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"none"})
#' @param scores Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"regression"})
#' @param fm Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"ml"})
#'
#' @return A list containing:
#' \item{Theta}{A size \eqn{q \times p} coefficient matrix}
#' \item{B}{A size \eqn{q \times k} factor loading matrix}
#' \item{E.Z}{A size \eqn{n \times k} latent factor score matrix}
#' \item{diag.Psi}{A length-\eqn{q} vector of uniquenesses (diagonal of \eqn{\Psi}).}
#'
#' @export
DrFARM.one <- function(X, Y, Theta0, precM, k,
                       lambda1, lambda2, K = NULL,
                       C = NULL,
                       standardize = TRUE,
                       thres = 1e-4,
                       rotate = "none",
                       scores = "regression",
                       fm = "ml") {

  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  if (!is.null(K)) {
    eigen.res <- eigen(K)
    U <- eigen.res$vectors
    d <- eigen.res$values
    X <- t(U) %*% X
    Y <- t(U) %*% Y
  } else {
    d <- rep(1, n)
  }

  # Debias the sparse estimator used for initial value
  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  E.star <- Y - X %*% Theta0.db.t

  # Factor analysis (initial)
  fa.res <- fa(E.star, nfactors = k, rotate = rotate, scores = scores, fm = fm, covar = TRUE)
  B <- fa.res$loadings
  diag.Psi <- fa.res$uniquenesses
  Psi.Inv <- diag(1/diag.Psi)
  PsiInv.B <- Psi.Inv %*% B
  Bt.PsiInv.B <- t(B) %*% Psi.Inv %*% B

  InvList <- list()
  for (i in 1:n) {
    sandInv <- solve(diag(1/d[i], k, k) + Bt.PsiInv.B)
    InvList[[i]] <- Psi.Inv - PsiInv.B %*% sandInv %*% t(PsiInv.B)
  }

  iter <- 0
  prev.loss <- 1e300
  diff <- 1e300

  Theta.t <- t(Theta0)
  prev.Theta.t <- Theta.t
  prev.B <- B
  prev.Psi <- diag.Psi

  while (diff > thres) {
    iter <- iter + 1
    e <- Y - X %*% Theta.t

    WList <- list()
    for (i in 1:n) {
      WList[[i]] <- d[i] * t(B) %*% InvList[[i]]
    }
    E.zt <- matrix(NA, k, n)
    for (i in 1:n) {
      E.zt[,i] <- WList[[i]] %*% e[i,]
    }

    E.zzt.List <- list()
    for (i in 1:n) {
      E.zzt.List[[i]] <- d[i] * (diag(1, k, k) - WList[[i]] %*% B) + E.zt[,i] %*% t(E.zt[,i])
    }

    E.zzt <- Reduce("+", E.zzt.List)
    E.zzt.inv <- solve(E.zzt)

    Y.aug <- Y - t(B %*% E.zt)
    if (is.null(C)) {
      Theta.t <- remMap(X.m = X, Y.m = Y.aug, lamL1 = lambda1, lamL2 = lambda2,
                        phi0 = prev.Theta.t, C.m = NULL, sigma = diag.Psi)$phi
    } else {
      Theta.t <- remMap(X.m = X, Y.m = Y.aug, lamL1 = lambda1, lamL2 = lambda2,
                        phi0 = prev.Theta.t, C.m = t(C), sigma = diag.Psi)$phi
    }
    e <- Y - X %*% Theta.t
    Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n
    E.star <- Y - X %*% Theta.db.t

    B <- t(E.zt %*% E.star) %*% E.zzt.inv
    EET.star <- t(E.star) %*% E.star

    diag.Psi <- diag(EET.star - B %*% E.zt %*% E.star)/n
    Psi.Inv <- diag(1/diag.Psi)
    PsiInv.B <- Psi.Inv %*% B
    Bt.PsiInv.B <- t(B) %*% Psi.Inv %*% B

    InvList <- list()
    for (i in 1:n) {
      sandInv <- solve(diag(1/d[i], k, k) + Bt.PsiInv.B)
      InvList[[i]] <- Psi.Inv - PsiInv.B %*% sandInv %*% t(PsiInv.B)
    }

    loss <- sum(e^2 * matrix(1/diag.Psi, n, q, byrow = TRUE)) / 2 +
      lambda1 * sum(abs(Theta.t)) +
      lambda2 * sum(sqrt(rowSums(Theta.t^2))) +
      n * sum(log(diag.Psi)) / 2

    print(c("M-step", "iter", iter, loss))

    diff <- prev.loss - loss
    prev.loss <- loss
    if (diff > 0) {
      prev.Theta.t <- Theta.t
      prev.B <- B
      prev.Psi <- diag.Psi
    } else if (diff < 0) {
      iter <- iter - 1
      break
    }
    Theta.t <- prev.Theta.t
    B <- prev.B
    diag.Psi <- prev.Psi
  }

  return(list(Theta = t(Theta.t), B = B, E.Z = t(E.zt), diag.Psi = diag.Psi))
}

#' Compute EBIC for \code{DrFARM}
#'
#' Computes the Extended Bayesian Information Criterion (EBIC) for a given DrFARM coefficient matrix \eqn{\Theta}
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta A \eqn{q \times p} coefficient matrix
#' @param B A size \eqn{q \times k} factor loading matrix
#' @param E.Z A size \eqn{n \times k} latent factor score matrix
#' @param diag.Psi A length-\eqn{q} vector of uniquenesses (diagonal of \eqn{\Psi})
#' @param K An optional size \eqn{n \times n} kinship matrix. Default is \code{NULL}
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{1}
#'
#' @return A numeric EBIC value
#'
#' @export
DrFARM.EBIC <- function(X, Y, Theta, B, E.Z, diag.Psi, K = NULL,
                        standardize = TRUE, gamma = 1) {

  n <- dim(Y)[1]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  if (!is.null(K)) {
    eigen.res <- eigen(K)
    U <- eigen.res$vectors
    d <- eigen.res$values
    X <- t(U) %*% X
    Y <- t(U) %*% Y
  }

  E <- Y - X %*% t(Theta) - E.Z %*% t(B)
  ratio <- diag(t(E) %*% E) / diag.Psi
  EBIC <- sum(ratio) + n * sum(log(diag.Psi)) +
    log(n) * sum(Theta != 0) +
    2 * gamma * logsumchoose(q, rowSums(t(Theta) != 0))

  return(EBIC)
}

#' Fit all \code{DrFARM} models at once
#'
#' Searches over a 2D grid of (\code{lambda1}, \code{lambda2}) pairs, computes \code{EBIC} for each, and
#' returns the \eqn{q x p} coefficient matrix that minimizes the EBIC, along with corresponding
#' factor loadings, factor scores, and uniquenesses
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta0 A size \eqn{q \times p} initial coefficient matrix
#' @param precM A size \eqn{p \times p} precision matrix
#' @param k A positive integer specifying the number of latent factors used in DrFARM
#' @param lambda1.opt The chosen lasso tuning parameter from \code{remMap.whole}
#' @param lambda2.opt The chosen group-lasso tuning parameter from \code{remMap.whole}
#' @param K An optional size \eqn{n \times n} kinship matrix. Default is \code{NULL}
#' @param C A size \eqn{q \times p} integer matrix specifying model inclusion and penalty rules:
#'   \itemize{
#'     \item 0 = excluded from model
#'     \item 1 = penalized
#'     \item 2 = unpenalized
#'   }
#'   If \code{null}, all entries are penalized (1) by default
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#' @param gamma A numeric value in \([0, 1]\) for the EBIC hyperparameter. When \code{gamma = 0}, EBIC reduces to the ordinary BIC. Default is \code{1}
#' @param thres A numeric value for the convergence threshold. Default is \eqn{1 \times 10^{-4}}
#' @param rotate Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"none"})
#' @param scores Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"regression"})
#' @param fm Character string specifying factor analysis options, consistent with the \pkg{psych} package (e.g., \code{"ml"})
#'
#' @return A list containing:
#' \item{Theta}{The EBIC-chosen \eqn{q \times p} coefficient matrix}
#' \item{B}{A size \eqn{q \times k} factor loading matrix}
#' \item{E.Z}{A size \eqn{n \times k} latent factor score matrix}
#' \item{diag.Psi}{A length-\eqn{q} vector of uniquenesses (diagonal of \eqn{\Psi}).}
#' \item{lambda1.opt}{Optimal \eqn{\lambda_1} value}
#' \item{lambda2.opt}{Optimal \eqn{\lambda_2} value}
#'
#' @export
DrFARM.whole <- function(X, Y, Theta0, precM, k,
                         lambda1.opt, lambda2.opt,
                         K = NULL, C = NULL,
                         standardize = TRUE,
                         gamma = 1,
                         thres = 1e-4,
                         rotate = "none",
                         scores = "regression",
                         fm = "ml") {

  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  E.star <- Y - X %*% Theta0.db.t
  fa.res <- fa(E.star, nfactors = k, rotate = rotate, scores = scores, fm = fm, covar = TRUE)
  diag.Psi <- fa.res$uniquenesses

  Psi.range <- unname(summary(diag.Psi)[-4])
  DrFARM.lambda.grid <- expand.grid(lambda1.opt / Psi.range, lambda2.opt / Psi.range)
  colnames(DrFARM.lambda.grid) <- c("lambda1", "lambda2")

  n.lambda.sq <- dim(DrFARM.lambda.grid)[1]
  ls <- vector("list", n.lambda.sq)

  for (i in 1:n.lambda.sq) {
    ls[[i]] <- DrFARM.one(X, Y, Theta0, precM, k,
                          DrFARM.lambda.grid[i,1], DrFARM.lambda.grid[i,2],
                          K = K, C = C, standardize = FALSE,
                          thres = thres, rotate = rotate,
                          scores = scores, fm = fm)
  }

  EBIC.vec <- rep(NA, n.lambda.sq)
  for (i in 1:n.lambda.sq) {
    E <- Y - X %*% t(ls[[i]]$Theta) - ls[[i]]$E.Z %*% t(ls[[i]]$B)
    ratio <- diag(t(E) %*% E) / ls[[i]]$diag.Psi
    EBIC.vec[i] <- sum(ratio) + n * sum(log(ls[[i]]$diag.Psi)) +
      log(n) * sum(ls[[i]]$Theta != 0) +
      2 * gamma * logsumchoose(q, rowSums(t(ls[[i]]$Theta) != 0))
  }

  opt.idx <- which.min(EBIC.vec)

  return(list(Theta = ls[[opt.idx]]$Theta,
              B = ls[[opt.idx]]$B,
              E.Z = ls[[opt.idx]]$E.Z,
              diag.Psi = ls[[opt.idx]]$diag.Psi,
              lambda1.opt = DrFARM.lambda.grid[opt.idx,1],
              lambda2.opt = DrFARM.lambda.grid[opt.idx,2]))
}

#' Compute entrywise p-values for DrFARM
#'
#' Computes entrywise p-values for DrFARM derived from outer-debiasing
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta A size \eqn{q \times p} initial coefficient matrix
#' @param B A size \eqn{q \times k} factor loading matrix
#' @param E.Z A size \eqn{n \times k} latent factor score matrix
#' @param precM A size \eqn{p \times p} precision matrix
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#'
#' @return A size \eqn{q \times p} matrix of p-values
#'
#' @export
entry.pvalue <- function(X, Y, Theta, B, E.Z, precM,
                         standardize = TRUE) {

  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  Y.aug <- Y - E.Z %*% t(B)

  Theta.t <- t(Theta)
  Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n

  s <- colSums(Theta.t != 0)
  sses <- diag(t(Y - X %*% Theta.t - E.Z %*% t(B)) %*% (Y - X %*% Theta.t - E.Z %*% t(B)))
  Psi.star <- sses / (n - s)

  CovM <- t(X) %*% X / n
  sqrtPhi <- sqrt(diag(precM %*% CovM %*% t(precM)))

  Z <- matrix(NA, p, q)
  for (i in 1:q) {
    Z[,i] <- sqrt(n) * Theta.db.t[,i] / (sqrt(Psi.star[i]) * sqrtPhi)
  }

  pval <- 2 * pnorm(-abs(Z))
  return(t(pval))
}

#' Compute group p-values for DrFARM
#'
#' Computes group-level p-values via the Cauchy combination test, providing a single p-value per predictor
#' across all outcomes.
#'
#' @param X A size \eqn{n \times p} matrix of predictors (e.g., genetic variants). Missing values are not allowed
#' @param Y A size \eqn{n \times q} matrix of outcomes (e.g., continuous omics traits). Missing values are not allowed
#' @param Theta A size \eqn{q \times p} initial coefficient matrix
#' @param B A size \eqn{q \times k} factor loading matrix
#' @param E.Z A size \eqn{n \times k} latent factor score matrix
#' @param precM A size \eqn{p \times p} precision matrix
#' @param standardize A logical indicating whether to standardize \code{X} and \code{Y} by column. Default is \code{TRUE}
#'
#' @return A length-\eqn{p} vector of p-values, where each p-value corresponds to one predictor across all outcomes
#'
#' @export
pleio.pvalue <- function(X, Y, Theta, B, E.Z, precM,
                         standardize = TRUE) {

  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]

  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }

  Y.aug <- Y - E.Z %*% t(B)

  Theta.t <- t(Theta)
  Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n

  s <- colSums(Theta.t != 0)
  sses <- diag(t(Y - X %*% Theta.t - E.Z %*% t(B)) %*% (Y - X %*% Theta.t - E.Z %*% t(B)))
  Psi.star <- sses / (n - s)

  CovM <- t(X) %*% X / n
  sqrtPhi <- sqrt(diag(precM %*% CovM %*% t(precM)))

  Z <- matrix(NA, p, q)
  for (i in 1:q) {
    Z[,i] <- sqrt(n) * Theta.db.t[,i] / (sqrt(Psi.star[i]) * sqrtPhi)
  }

  pval <- 2 * pnorm(-abs(Z))

  T <- rowSums(1 / tan(pi * pval)) / q
  pval2 <- 2 * pcauchy(-abs(T))
  return(pval2)
}
