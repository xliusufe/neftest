\name{rIGauss}
\alias{rIGauss}
\title{
  Generates random numbers from Inverse Gaussian distribution.
}
\description{
  Generates random numbers from Inverse Gaussian distribution.
}

\usage{rIGauss(n, mu = 1.0, lambda = 1.0)
}

\arguments{
  \item{n}{Number of random numbers to be generated.}

  \item{mu}{shape parameters \eqn{\mu}. Default is 1.}

  \item{lambda}{scale parameters \eqn{\lambda}. Default is 1.}
}

\details{
  The probability density function:

    \eqn{
      f(x; \mu, \lambda) = \sqrt{\frac{\lambda}{2\pi x^3}}\exp\left(-\frac{\lambda(x - \mu)^2}{2\mu^2x}\right),  x>0.
    }
}

\value{
  \item{x}{The random numbers from Inverse Gaussian distribution.}
}


\examples{
  x = rIGauss(n = 10, mu = 1.0, lambda = 1.0)
  x
}