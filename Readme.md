# neftest
R package "neftest" for the goodness of fit tests based on zero regression characterizations of Tweedie, Bar-Lev and Enis class of distributions. Provide p-values of test statistics according to several specified distributions, which is based on the method introduced by Authors (2021).

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/neftest")

# Usage

   - [x] [neftest-manual.pdf](https://github.com/xliusufe/neftest/blob/master/inst/neftest-manual.pdf) ---------- Details of the usage of the package.
# Example

    ## Poisson
    library(neftest)
    n   <- 100
    distr <- "Poisson"
    x     <- rpois(n,lambda = 1)
    pval  = pvals(x, distr)

    pval

    ## Gamma
    library(neftest)
    n     <- 100
    distr <- "Gamma"
    x     <- rgamma(n, shape = 1, rate = 1)
    pval  = pvals(x, distr)

    pval


    ## Inverse Gaussian
    library(neftest)
    n     <- 100
    distr <- "Inverse Gaussian"

    x     <- rIGauss(n, mu = 1, lambda = 1)
    pval  = pvals(x, distr)

    pval

# References

Authors (2021). Goodness of fit tests based on zero regression characterizations of Tweedie, Bar-Lev and Enis class of distributions. Manuscript.

# Development
This R package is developed by Panpan Ren and Xu Liu (liu.xu@sufe.edu.cn).
