\name{bivd}
\alias{bivd}
\title{ Manipulation of Bivariate Distributions }
\description{
  Creates an object of class \code{bivd} encapsulating the building blocks 
of a bivariate distribution
}
\usage{
bivd(cop, marginX = "unif", marginY = marginX, param.marginX = c(0, 1), param.marginY = param.marginX)
}
\arguments{
  \item{cop}{ Object of class copula }
  \item{marginX}{ Character string giving the name of a univariate distribution }
  \item{marginY}{ Character string giving the name of a univariate distribution }
  \item{param.marginX}{ Numeric vector giving the parameters of the marginal distribution of X }
  \item{param.marginY}{ Numeric vector giving the parameters of the marginal distribution of Y }
}
\value{
  Object of class bivd
}
\seealso{ \code{\link{contour.dbivd}}, \code{\link{persp.dbivd}},\code{\link{contour.pbivd}}, \code{\link{persp.pbivd}} }
\examples{
BIV <- bivd(gumbel.copula(1.4), "norm", "t", c(3,4), 3)
persp.dbivd(BIV)

par(mfrow=c(1,2))
contour.dbivd(BIV)
contour.pbivd(BIV)
par(mfrow=c(1,1))
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 3
}
\keyword{distribution}
