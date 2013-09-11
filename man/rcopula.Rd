\name{rcopula}
\alias{rcopula}

\title{ Random generator for bivariate samples from a copula }
\description{
 Generate bivariate samples from a copula
}
\usage{
rcopula(copula, n)
}
\arguments{
  \item{copula}{ An object of class copula }
  \item{n}{ An integer giving the size of the desired random sample }
}
\value{
  A list (or data frame) comprising
  \item{x }{Numeric vector of length n giving samples from the first marginal}
  \item{y }{{Numeric vector of length n giving samples from the first marginal}}
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{ R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 3}
\examples{
N <- 5000
SD <- rcopula(gumbel.copula(1.5),N)
TAU1.sim <- qexp(0.5, SD$x)
TAU2.sim <- qexp(0.9, SD$y)
P1 <- mean(TAU1.sim < TAU2.sim)
TT <- 1
P2 <- mean ( (TAU1.sim > TT) & (TAU2.sim > TT))
}
\keyword{ distribution }
