\name{pvals}
\alias{pvals}
\title{
  The p-value of the test based on the test statistic \eqn{T{nw}}
}
\description{
  Compute the p-value of the test based on the test statistic \eqn{T_{nw}}.
}

\usage{pvals(x, distr = "Poisson", bootstrap = FALSE, B = 1000, 
                weight = 1, a = 1.0, max.iter = 100, tol = 1e-8)
}

\arguments{
  \item{x}{A length \eqn{n} vector of input data.}

  \item{distr}{The true distribution including Poisson distribution \code{distr = "Poisson"}, Gamma distribution \code{distr = "Gamma"} and Inverse Gaussian distribution \code{distr = "Inverse Gaussian"}. Default is \code{distr = "Poisson"}.}

  \item{bootstrap}{logical. The bootstrap method is used to compute the p-value if FALSE (default), and the maximum likelihood method otherwise.}

  \item{B}{Number of bootstrap samples. Default is 1000.}

  \item{weight}{The weight functions including Laplace distribution \code{weight = 1} and normal distribution \code{weight = 2}. 
  \code{weight = 3} denotes the test statistic to be limit of "a" when it goes to \eqn{+\infty} in Remark 5 or the test statistic \eqn{ST_{n,w}} in Remark 4. 
  Default is \code{weight = 1}.}

  \item{a}{The parameter of the weight functions. Default is 1.0. See details in the paper}

  \item{max.iter}{The maximum number of iterations in Newton method. Default is 100.}

  \item{tol}{The precision of the Newton method. Default is \code{1e-8}.}
}

\value{
  \item{pval}{The p-value of the test statistic.}
}


\references{Authors (2021). Goodness of fit tests based on zero regression characterizations of Tweedie, Bar-Lev and Enis class of distributions. Manuscript.
}

\examples{
  # Poisson
  n     <- 100
  distr <- "Poisson"
  x     <- rpois(n,lambda = 1)
  pval  <- pvals(x, distr)

  pval


  # Poisson with Bootstrap
  n     <- 100
  B     <- 500
  distr <- "Poisson"
  x	    <- rpois(n,lambda = 1)
  pval  <- pvals(x, distr, bootstrap = TRUE, B = B)

  pval



  # Gamma
  n     <- 100
  distr <- "Gamma"
  x     <- rgamma(n, shape = 1, rate = 1)
  pval  <- pvals(x, distr)

  pval


  # Gamma with Bootstrap
  n     <- 100
  B     <- 500
  distr <- "Gamma"
  x	    <- rgamma(n, shape = 1, rate = 1)
  pval  <- pvals(x, distr, bootstrap = TRUE, B = B)

  pval



  # Inverse Gaussian
  n     <- 100
  distr <- "Inverse Gaussian"
  x     <- rIGauss(n, mu = 1, lambda = 1)
  pval  <- pvals(x, distr)

  pval



  # Inverse Gaussian with Bootstrap
  n     <- 100
  B     <- 500
  distr <- "Inverse Gaussian"
  x	    <- rIGauss(n, mu = 1, lambda = 1)
  pval  <- pvals(x, distr, bootstrap = TRUE, B = B)

  pval
}