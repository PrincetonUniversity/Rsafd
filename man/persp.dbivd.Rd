\name{persp.dbivd}
\alias{persp.dbivd}
\alias{persp.pbivd}
\alias{persp.dcopula}
\alias{persp.pcopula}
\alias{persp.dbivd,bivd-method}
\alias{persp.dbivd,gpd.bivd-method}
\alias{persp.pbivd,bivd-method}
\alias{persp.pbivd,gpd.bivd-method}
\alias{persp.dcopula,copula-method}
\alias{persp.dcopula,empirical.copula-method}
\alias{persp.pcopula,copula-method}
\alias{persp.pcopula,empirical.copula-method}
\title{
Function producing a surface plot of the density or the cdf of a bivariate distribution
}
\description{
Produce a surface plot of the density or of the cdf of a bivariate distribution
(\code{persp.dbivd}, \code{persp.pbivd}) or of a copula (\code{persp.dcopula},
\code{persp.pcopula}).
}
\usage{
persp.dbivd(dist, n = 50, xlim = NA, ylim = NA, ...)
persp.pbivd(dist, n = 50, xlim = NA, ylim = NA, ...)
persp.dcopula(cop, n = 50, ...)
persp.pcopula(cop, n = 50, ...)
}
\arguments{
  \item{dist}{
bivariate distribution
}
  \item{n}{
integer giving the number of points in each dimension of the grid over which the plot is performed
}
  \item{xlim}{
limits of the x - variable}
  \item{ylim}{
limits of the y - variable}
  \item{cop}{
a copula (object of class "copula" or "empirical.copula")}
  \item{\dots}{
additional arguments}
}
\value{
Produce a surface plot of the density or the cdc of the bivariate distribution}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag
}
\seealso{
See Also as \code{\link{contour.dbivd}}
}
\examples{
BIV <- bivd(gumbel.copula(1.4), "norm", "t", c(3,4), 3)
persp.dbivd(BIV)}
