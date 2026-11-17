\name{contour.dbivd}
\alias{contour.dbivd}
\alias{contour.pbivd}
\alias{contour.dcoupla}
\alias{contour.pcopula}
\alias{contour.dcopula}
\alias{contour.dbivd,bivd-method}
\alias{contour.dbivd,gpd.bivd-method}
\alias{contour.pbivd,bivd-method}
\alias{contour.pbivd,gpd.bivd-method}
\alias{contour.dcopula,copula-method}
\alias{contour.dcopula,empirical.copula-method}
\alias{contour.pcopula,copula-method}
\alias{contour.pcopula,empirical.copula-method}
\title{
Function producing a contour plot of the level sets of the density / cdf of a bivariate distribution or a copula
}
\description{
Generic method for the classes bivd (bivariate distribution), copula and empirical.copula plotting the level sets of the density of cdf (contour.pbivd, contour.pcopula) or the density (contour.dbivd, contour.dcopula) of the distribution
}
\usage{
\S4method{contour.dbivd}{bivd}(dist, n = 50, xlim = NA, ylim = NA, ...)
\S4method{contour.dbivd}{gpd.bivd}(dist, n = 50, xlim = NA, ylim = NA, ...)
\S4method{contour.pbivd}{bivd}(dist, n = 50, xlim = NA, ylim = NA, ...)
\S4method{contour.pbivd}{gpd.bivd}(dist, n = 50, xlim = NA, ylim = NA, ...)
\S4method{contour.dcopula}{copula}(cop, n = 50, ...)
\S4method{contour.dcopula}{empirical.copula}(cop, n = 50, ...)
\S4method{contour.pcopula}{copula}(cop, n = 50, ...)
\S4method{contour.pcopula}{empirical.copula}(cop, n = 50, ...)
}
\arguments{
  	\item{dist}{a bivariate distribution (object of class "bid") or a copula (object of class "copula" or "empirical.copula")}
  	\item{n}{number of points on each axis forming the rectangular grid used in the discretization of the	cdf or density}
  	\item{xlim}{limits of the plot on the x-axis}
  	\item{ylim}{limits of the plot on the y-axis}
  	\item{cop}{a copula (object of class "copula" or "empirical.copula")}
  	\item{\dots}{additional arguments}
}
\value{
Used for its side effects}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{
See Also  \code{\link{bivd}}, \linkS4class{copula}, \code{\link{persp.dbivd}}
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
