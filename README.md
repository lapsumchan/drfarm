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
where `n = 500`, `p = 10` and `q = 5`.

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

Theta0
> Theta0
     [,1] [,2]         [,3] [,4] [,5] [,6] [,7]       [,8] [,9]      [,10]
[1,]    0    0 -0.208809944    0    0    0    0 0.00000000    0 -0.1635109
[2,]    0    0  0.110882820    0    0    0    0 0.00000000    0 -0.2878153
[3,]    0    0  0.005000023    0    0    0    0 0.23695892    0  0.0000000
[4,]    0    0 -0.103653563    0    0    0    0 0.00000000    0  0.2709923
[5,]    0    0  0.000000000    0    0    0    0 0.04195427    0  0.1935116
```
As a result, `Theta0` estimated by remMap is not all the same scale as the true coefficient matrix:
```
> Theta0
     [,1] [,2]         [,3] [,4] [,5] [,6] [,7]       [,8] [,9]      [,10]
[1,]    0    0 -0.208809944    0    0    0    0 0.00000000    0 -0.1635109
[2,]    0    0  0.110882820    0    0    0    0 0.00000000    0 -0.2878153
[3,]    0    0  0.005000023    0    0    0    0 0.23695892    0  0.0000000
[4,]    0    0 -0.103653563    0    0    0    0 0.00000000    0  0.2709923
[5,]    0    0  0.000000000    0    0    0    0 0.04195427    0  0.1935116
```

Next, we estimate the precision matrix using `precM`:
```
precM <- precM(X)
```
and default, `precM` estimates the precision matrix using glasso (the recommended approach in our paper).

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

The estimated coefficient matrix is 
Finally, the `q x p` entrywise (`pval1`) and length `p` pleiotropic (`pval2`) *p*-values can simply be obtained using:
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
