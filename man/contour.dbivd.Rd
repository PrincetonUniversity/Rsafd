\name{contour.dbivd}
\alias{contour.dbivd}
\alias{contour.pbivd}
\alias{contour.dcoupla}
\alias{contour.pcopula}
\title{
Function producing a contour plot of the level sets of the density of cdf of a bivariate distribution or a copula
}
\description{
Generic method for the classes bivd (bivariate distribution), copula and empirical.copula plotting the level sets of the density of cdf (contour.pbivd, contour.pcopula) or the density (contour.dbivd, contour.dcopula) of the distribution
}
\usage{
contour.dbivd(dist, n = 50, xlim = NA, ylim = NA, ...)
}
\arguments{
  	\item{dist}{a bivariate distribution (object of class "bid") or a copula (object of class "copula" or "empirical.copula")}
  	\item{n}{number of points on each axis forming the rectangular grid used in the discretization of the	cdf or density}
  	\item{xlim}{limits of the plot on the x-axis}
  	\item{slim}{limits of the plot on the y-axis}
}
\value{
Used for its side effects}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}

\seealso{
See Also  \code{\link{bivd}}, \code{\link{copula}},
}
\examples{
BIV1 <- bivd(gumbel.copula(1.4), "norm", "t", c(3,4), 3)
persp.dbivd(BIV1)

par(mfrow=c(1,2))
contour.dbivd(BIV1)
contour.pbivd(BIV1)
par(mfrow=c(1,1))
}
\keyword{datasets}
