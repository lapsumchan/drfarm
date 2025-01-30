# Debiased-regularized factor analysis regression model (DrFARM)

### Overview
This repository provides a demonstration on how to use the `R` package `drfarm`

# System Requirements

## Software Requirements

### OS Requirements

The source code has been tested on Microsoft's Windows 10 operating system and Linux (Ubuntu 18.04). The source code should be compatible with Windows, Mac, and Linux operating systems.

Before using `drfarm` package, users should have `R` version 4.3.0 or higher.

### Installation  

First, install `camp` from GitHub using `devtools`:  

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

`drfarm.dat` contains a small simulated toy example with sample size `n = 500`, `p = 10` variants and `q = 5` traits, with 3 pleiotropic variants (variant #3, #8 and #10): `X` (`n x p` variants matrix) and `Y` (`n x q` trait matrix), as well as the ground truth `p x q` coefficient matrix `Theta.t`. Notice that there are `k = 2` underlying latent factors that contribute to the 5 traits.
