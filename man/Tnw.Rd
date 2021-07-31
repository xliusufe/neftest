\name{Tnw}
\alias{Tnw}
\title{
  The test statistic for testing if a distribution is a TBE(\eqn{\gamma0})
}
\description{
  Compute the test statistic \eqn{T_{nw}}.
}

\usage{Tnw(x, gamma0 = 1, weight = 1, a = 1.0)
}

\arguments{
  \item{x}{A length \eqn{n} vector of input data.}

  \item{gamma0}{The power parameter in TBE(\eqn{\gamma_0}). Default is 1.}

  \item{weight}{The weight functions including Laplace distribution \code{weight = 1} and normal distribution \code{weight = 2}. 
  \code{weight = 3} denotes the test statistic to be limit of "a" when it goes to \eqn{+\infty} in Remark 5 or the test statistic \eqn{ST_{n,w}} in Remark 4. 
  Default is \code{weight = 1}.}

  \item{a}{The parameter of the weight function. Default is 1.0.}
}

\value{
  \item{Tn}{The test statistic \eqn{T_{nw}}.}
}

\references{Authors (2021). Goodness of fit tests based on zero regression characterizations of Tweedie, Bar-Lev and Enis class of distributions. Manuscript.
}

\examples{
  # Poisson
  n   <- 100
  x   <- rpois(n,lambda = 1)
  B   <- 1000
  Tn  <- Tnw(x, gamma0 = 1, weight = 1, a = 1.0)
  Tb 	= rep(NA, B)
  lambdahat 	= mean(x)
  for(b in 1:B){
    xb 		<- rpois(n, lambda = lambdahat)
    Tb[b] <- Tnw(xb, gamma0 = 1, weight = 2, a = 1.0)
  }
  pval = mean(Tb > Tn)

  # Gamma
  n   <- 100
  x   <- rgamma(n, shape = 1, rate = 1)
  Tn  <- Tnw(x, gamma0 = 1, a = 1.0)

  # Inverse Gaussian
  n   <- 100
  x   <- rIGauss(n, mu = 1, lambda = 1)
  B   <- 1000
  Tn  <- Tnw(x, gamma0 = 1, a = 1.0)
  Tb 	= rep(NA, B)
  nuhat 	= mean(x)
  lambdahat = (mean(1/x)-nuhat^(-1))^(-1)
  for(b in 1:B){
    xb 		<- rIGauss(n, mu = nuhat, lambda = lambdahat)
    Tb[b] <- Tnw(xb, gamma0 = 1, weight = 1, a = 1.0)
  }
  pval = mean(Tb > Tn)
}