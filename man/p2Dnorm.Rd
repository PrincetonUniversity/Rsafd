\name{p2Dnorm}
\alias{p2Dnorm}
\title{Cumulative Distribution Function of a Bivariate Normal Distribution}
\description{
Compute the value of the cumulative distribution function of a bivariate normal
distribution at a point of the plane.
}
\usage{
p2Dnorm(x, mean = rep(0, 2), std = rep(1, 2), rho = 0)
}
\arguments{
  \item{x}{numeric vector of length 2 giving the point at which the cdf is computed}
  \item{mean}{numeric vector of length 2 giving the means}
  \item{std}{numeric vector of length 2 giving the (positive) standard deviations}
  \item{rho}{number between -1 and +1 giving the correlation coefficient}
}
\details{
The cdf is computed by numerical integration (with \code{\link[stats]{integrate}}) of the
conditional cdf of the second component given the first.
}
\value{
A number giving the value of the cdf at \code{x}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\code{\link{dmvnorm}}, \code{\link{pmvnorm}}}
\examples{
p2Dnorm(c(0, 0), rho = 0.5)
}
\keyword{distribution}
\keyword{multivariate}
