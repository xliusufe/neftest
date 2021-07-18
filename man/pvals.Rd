\name{pvals}
\alias{pvals}
\title{
  The p-value of the test based on the test statistic \eqn{T{nw}}
}
\description{
  Compute the p-value of the test based on the test statistic \eqn{T_{nw}}.
}

\usage{pvals(x, distr="Poisson", bootstrap = FALSE, B = 1000, signif = 0.05, 
                weight="normal", a = 1.0, max.iter = 100, tol = 1e-8)
}

\arguments{
  \item{x}{A length \eqn{n} vector of input data.}

  \item{distr}{The true distribution including Poisson distribution \code{distr = "Poisson"}, Gamma distribution \code{distr = "Gamma"} and Inverse Gaussian distribution \code{distr = "Inverse Gaussian"}. Default is \code{distr = "Poisson"}.}

  \item{bootstrap}{logical. Bootstrap method is used to compute the p-value if FALSE (default), and the maximum likelihood method otherwise.}

  \item{B}{Number of bootstrap samples. Default is 1000.}

  \item{signif}{The significance level of the test. Default is 0.05.}

  \item{weight}{The weight functions including normal distribution \code{weight = "normal"} and Laplace distribution \code{weight = "laplace"}. Default is \code{weight = "normal"}.}

  \item{a}{The parameter of the weight function. Default is 1.0. See details in the paper}

  \item{max.iter}{The maximum number of iterations in Newton method. Default is 100.}

  \item{tol}{The precision of the Newton method. Default is \code{1e-8}.}
}

\value{
  \item{pval}{The p-value of the test.}
}


\references{Authors (2021). Goodness of fit tests based on zero regression characterizations of Tweedie, Bar-Lev and Enis class of distributions. Manuscript.
}

\examples{ 
  # Poisson
  n   <- 100
  NS  <- 500
  distr <- "Poisson"
  pval  <- rep(NA, NS)
  set.seed(100861)
  for (i in 1:NS) {
    x       <- rpois(n,lambda = 1)
    pval[i] = pvals(x, distr, bootstrap = FALSE, B = 1000, signif = 0.05)
  }
  pvalue = mean(pval<0.05)


  # Poisson
  n   <- 100
  NS  <- 500
  B   <- 500
  distr <- "Poisson"
  pval  <- rep(NA, NS)
  for(i in 1:NS){
    x	      <- rpois(n,lambda = 1)
    pval[i] = pvals(x, distr, bootstrap = TRUE, B = B, signif = 0.05)
  }
  Pval = mean(pval<0.05)



  # Gamma
  n     <- 100
  NS    <- 500
  distr <- "Gamma"
  pval  <- rep(NA, NS)
  set.seed(100861)
  for (i in 1:NS) {
    x       <- rgamma(n, shape = 1, rate = 1)
    pval[i] = pvals(x, distr, bootstrap = FALSE, B = 1000, signif = 0.05)
  }
  pvalue = mean(pval<0.05)


  # Gamma
  n     <- 100
  NS    <- 500
  B     <- 500
  distr <- "Gamma"
  pval  <- rep(NA, NS)
  for(i in 1:NS){
    x	      <- rgamma(n,lambda = 1)
    pval[i] = pvals(x, distr, bootstrap = TRUE, B = B, signif = 0.05)
  }
  Pval = mean(pval<0.05)



  # Inverse Gaussian
  n     <- 100
  NS    <- 500
  distr <- "Inverse Gaussian"
  pval  <- rep(NA, NS)
  set.seed(100861)
  for (i in 1:NS) {
    x       <- rinvGauss(n, nu = 1, lambda = 1)
    pval[i] = pvals(x, distr, bootstrap = FALSE, B = 1000, signif = 0.05)
  }
  pvalue = mean(pval<0.05)



  # Inverse Gaussian
  n     <- 100
  NS    <- 500
  B     <- 500
  distr <- "Inverse Gaussian"
  pval  <- rep(NA, NS)
  for(i in 1:NS){
    x	      <- rinvGauss(n, nu = 1, lambda = 1)
    pval[i] = pvals(x, distr, bootstrap = TRUE, B = B, signif = 0.05)
  }
  Pval = mean(pval<0.05)
}