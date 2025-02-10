# Debiased-regularized factor analysis regression model (DrFARM)

### Overview
This repository provides a demonstration on how to use the `R` package `drfarm`

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using `drfarm` package, users should have `R` version 4.3.0 or higher.

### Installation  

First, install `drfarm` from GitHub using `devtools`:  

    # Install devtools if not already installed
    # install.packages("devtools") 
    devtools::install_github("lapsumchan/drfarm")
    
Installation should complete within a couple of minutes on a standard machine.

# Demo

Below is a toy example demonstrating the end-to-end workflow of DrFARM. For reproducibility, we use a small simulated dataset (`drfarm.dat`) included in the package.

To get started, load the necessary packages:

```
library(drfarm)
library(glmnet)
library(glasso)
library(psych)
```

The lazy-loaded `drfarm.dat` contains
- `X`: a `n` x `p` matrix of predictors (or variants)
- `Y`: a `n` x `q` matrix of outcomes (or traits)
- `Theta.t`: a `p` x `q` matrix of true (simulated) coefficients
where `n = 500`, `p = 10` and `q = 5`

```
X <- drfarm.dat$X
Y <- drfarm.dat$Y
Theta.t <- drfarm.dat$Theta.t
```

You can inspect it directly:
```
> t(Theta.t)
     [,1] [,2]      [,3] [,4] [,5] [,6] [,7]     [,8] [,9]     [,10]
[1,]    0    0 -4.383062    0    0    0    0 0.000000    0 -3.746484
[2,]    0    0  1.556498    0    0    0    0 0.000000    0 -4.612953
[3,]    0    0  0.000000    0    0    0    0 5.198627    0  0.000000
[4,]    0    0 -1.654552    0    0    0    0 0.000000    0  4.947609
[5,]    0    0  0.000000    0    0    0    0 1.166915    0  3.515037
```
which shows there are 3 "pleiotropic variants" (variant #3, #8 and #10).

To use DrFARM, we first obtain an initial sparse estimate (`Theta0`) using remMap<sup>[1]</sup>, which assumes standardization by default:
```
remMap.res <- remMap.whole(X, Y)
Theta0 <- remMap.res$Theta0

> Theta0
     [,1] [,2]         [,3] [,4] [,5] [,6] [,7]       [,8] [,9]      [,10]
[1,]    0    0 -0.208809944    0    0    0    0 0.00000000    0 -0.1635109
[2,]    0    0  0.110882820    0    0    0    0 0.00000000    0 -0.2878153
[3,]    0    0  0.005000023    0    0    0    0 0.23695892    0  0.0000000
[4,]    0    0 -0.103653563    0    0    0    0 0.00000000    0  0.2709923
[5,]    0    0  0.000000000    0    0    0    0 0.04195427    0  0.1935116
```
**Note: remMap estimates often differ from the true scale since it standardizes `X` and `Y` by default (likewise for DrFARM).**

Next, we need to estimate the precision matrix for the design matrix `X`. By default, `precM` uses glasso, which our paper recommends for this purpose.
```
precM <- precM(X)
```

With that, we have all the input needed for running DrFARM. We assume the number of latent factors (`k = 2`) is known, as used in generating the simulated data (For real data, you might determine `k` using exploratory graph analysis or other criteria).
```
k <- 2
DrFARM.res <- DrFARM.whole(X, Y, Theta0, precM, k, 
                           remMap.res$lambda1.opt, 
                           remMap.res$lambda2.opt)
Theta <- DrFARM.res$Theta;
B <- DrFARM.res$B; 
E.Z <- DrFARM.res$E.Z;
```

Once we have the `q` x `p` DrFARM coefficient matrix `Theta`, `q` x `k` loading matrix and `n` x `k` expected latent factors, we can compute:
1. Entrywise *p*-values for each variant-trait pair
2. Pleiotropy (group-level) *p*-values for each variant across all traits
```
pval1 <- entry.pvalue(X, Y, Theta, B, E.Z, precM)
pval2 <- pleio.pvalue(X, Y, Theta, B, E.Z, precM)
```

# Output
```
pval1
           [,1]       [,2]         [,3]      [,4]       [,5]      [,6]      [,7]         [,8]      [,9]        [,10]
[1,] 0.74585396 0.87490709 4.194237e-19 0.3924213 0.07035777 0.3057141 0.3449264 6.824145e-01 0.6802427 6.335812e-13
[2,] 0.07777188 0.54971679 2.629018e-07 0.7766920 0.64737382 0.1144596 0.7923100 5.357060e-01 0.5760429 5.770091e-26
[3,] 0.83586036 0.14528803 1.145534e-01 0.4878057 0.89476283 0.3889823 0.2641217 3.884086e-25 0.5107359 3.232748e-01
[4,] 0.71172948 0.09200038 1.539533e-07 0.9771724 0.28471456 0.8645441 0.2219014 1.185907e-01 0.6266486 7.943486e-35
[5,] 0.49088373 0.07114813 7.201191e-01 0.3354013 0.00654572 0.1026620 0.9432392 1.154125e-02 0.5977810 4.135523e-51

pval2
[1] 9.344327e-01 3.852209e-01 4.194237e-18 2.163034e-01 6.287407e-02 5.328321e-01 5.538712e-01 3.884086e-24 7.970215e-01 4.135523e-50
```

## Additional Notes

### Fitting `remMap` with Helper Functions

The functions `remMap.one()` (and similarly `DrFARM.one()`) fit a single (`\lambda_1`,`\lambda_2`) pair. A typical workflow is:

1. **Generate a tuning grid**  
   ```
   # Suppose this yields 10 x 10 = 100 grid cells
   remMap.lambda.grid <- remMap.grid(X, Y, standardize = TRUE)
   ```
2. **Fit a model for each grid cell**
   ```
   # Example using the i-th row of the grid:
   i <- 36
   Theta0.cand <- remMap.one(X, Y, lambda1 = remMap.lambda.grid[i,1], lambda2 = remMap.lambda.grid[i,2])
   ```
   In practice, you can loop over all grid rows (or use parallelization) to obtain 100 candidate solutions, each returning a `Theta0` matrix
3. **Select the best candidate via EBIC**
   ```
   EBIC <- remMap.EBIC(X, Y, Theta0.cand, standardize = TRUE)
   ```
   By default, `remMap.EBIC` sets `gamma = 1`, encouraging stronger sparsity. If you prefer a standard BIC (less sparse solution), set `gamma = 0`:
   ```
   BIC <- remMap.EBIC(X, Y, Theta0.cand, gamma = 0, standardize = TRUE)
   ```
   In this toy example, `i = 36` gives the smallest EBIC, so we choose that corresponding `Theta0` as our final remMap estimate.

### Fitting `DrFARM` with Helper Functions

After selecting the optimal `Theta0` from remMap (and corresponding (`lambda1.star` and `lambda2.star`), we can perform a similar procedure for DrFARM:

1. **Obtain a precision matrix**
   
   DrFARM requires estimating the precision matrix of `X`. By default, the function `precM()` uses the graphical lasso (as recommended in our paper):
   ```
   precM <- precM(X, method = "glasso", standardize = TRUE)
   ```
2. **Generate a tuning grid**
   
   `DrFARM.grid()` build a 5 x 5 grid around the chosen lasso (`lambda1.star`) and group-lasso (`lambda2.star`) parameters from remMap. It also requires specifying the number of latent factors `k`:
   ```
   Theta0 <- Theta0.cand
   lambda1.star <- remMap.lambda.grid[i,1]
   lambda2.star <- remMap.lambda.grid[i,2]
   
   DrFARM.lambda.grid <- DrFARM.grid(X, Y, Theta0, precM, k = 2, lambda1.opt = lambda1.star, lambda2.opt = lambda2.star)
3. **Fit a model for each grid cell**
   ```
   # Example with the i-th row of the grid:
   i <- 22
   DrFARM.one.res <- DrFARM.one(X, Y, Theta0, precM, k = 2, lambda1 = DrFARM.lambda.grid[i,1], lambda2 = DrFARM.lambda.grid[i,2])
   ```
4. **Select the best solution via EBIC**
   ```
   Theta <- DrFARM.one.res$Theta;
   B <- DrFARM.one.res$B;
   E.Z <- DrFARM.one.res$E.Z;
   diag.Psi <- DrFARM.one.res$diag.Psi;

   EBIC <- DrFARM.EBIC(X, Y, Theta, B, E.Z, diag.Psi)
   ```
   In this example, `i = 22` yields the smallest EBIC, so `Theta` here should be identical to that found by `DrFARM.whole()`.
   
### Special Case: (`lambda2 = 0`)
In remMap, setting `lambda2 = 0` removes the group-lasso penalty and yields a pure multivariate lasso model. Below is a minimal example:
```
# Generate a 2D grid (100 x 100), but only keep the unique lambda1 values
remMap.lambda.grid <- remMap.grid(X, Y, n.lambda = 100)
lasso.grid <- sort(unique(remMap.lambda.grid[, 1]))

# i = 63 yields the smallest EBIC
i <- 63
mlasso <- remMap.one(X, Y, lambda1 = lambda1.grid[i], lambda2 = 0)
```
A similar approach works for `DrFARM.one()` if you wish to consider a lasso-only penalty in the DrFARM setting as well.

# Citation

If you find `drfarm` useful, please cite:
> Chan, Lap Sum, et al. "DrFARM: Identification and inference for pleiotropic gene in GWAS." bioRxiv (2022): 2022-11.
