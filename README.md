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

To get started, load the necessary packages:

```
library(glmnet)
library(glasso)
library(psych)
```

**Note: To ensure consistency, a version of the simulated data (`drfarm.dat`) has been lazy-loaded within the `drfarm` package. Using this pre-loaded data ensures reproducibility across different systems.**

`drfarm.dat` contains a small simulated toy example with sample size `n = 500`, `p = 10` variants and `q = 5` traits, with 3 pleiotropic variants (variant #3, #8 and #10): `X` (`n x p` variants matrix) and `Y` (`n x q` trait matrix), as well as the ground truth `p x q` coefficient matrix `Theta.t:

```
X <- drfarm.dat$X
Y <- drfarm.dat$Y
Theta.t <- drfarm.dat$Theta.t
```

We estimate the initial value `Theta0` using remMap<sup>[1]</sup>:
```
remMap.res <- remMap.whole(X, Y)
Theta0 <- remMap.res$Theta0
```

Next, we estimate the precision matrix using `precM`:
```
precM <- precM(X)
```
By default, `precM` estimates the precision matrix using glasso (the recommended approach in our paper).

In this example, we assume the number of latent factors is known (`k = 2`), which is same as that of the number of latent factors used to generate the simulated data. Using `DrFARM.whole`:
```
k <- 2
DrFARM.res <- DrFARM.whole(X, Y, Theta0, precM, k, 
                           remMap.res$lambda1.opt, 
                           remMap.res$lambda2.opt)
Theta <- DrFARM.res$Theta;
B <- DrFARM.res$B; 
E.Z <- DrFARM.res$E.Z;
```
we obtain the estimated `q x p` sparse coefficient matrix `Theta`, `q x k` loading matrix `B` and `n x k` expected latent factors. These output are essential for the final step of calculating the entrywise *p*-values as well as the pleiotropic *p*-values.

For statistical inference, the `q x p` entrywise (`pval1`) and length `p` pleiotropic (`pval2`) *p*-values can simply be obtained using:
```
pval1 <- entry.pvalue(X, Y, Theta, B, E.Z, precM)
pval2 <- pleio.pvalue(X, Y, Theta, B, E.Z, precM)
```
