\name{empirical.copula}
\alias{empirical.copula}
\title{ Constructor of Empirical Copulas }
\description{
  Construct an object of class empirical.copula most likely to be used as argument of the 
function fit.copula.
}
\usage{
empirical.copula(x, y = NA)
}
\arguments{
  \item{x}{ Random sample in [0,1] assumed to be uniformly distributed }
  \item{y}{ Random sample in [0,1] assumed to be uniformly distributed }
}

\value{
  Returns an object of class empirical.copula capturing the dependencies
between the uniformly distributed samples x and y
}
\references{ R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 3}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
U <- runif(1024) # Generates a uniform sample of size 1024
V <- runif(1024) # Generates a uniform sample of size 1024

EMPCOP <- empirical.copula(U,V)
COP <- fit.copula(EMPCOP,"gumbel")
}

\keyword{distribution}
\keyword{multivariate}
