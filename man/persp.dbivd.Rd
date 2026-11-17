\name{persp.dbivd}
\alias{persp.dbivd}
\alias{persp.pbivd}
\title{
Function producing a surface plot of the density or the cdf of a bivariate distribution
}
\usage{
persp.dbivd(dist, n = 50, xlim = NA, ylim = NA, ...)
persp.pbivd(dist, n = 50, xlim = NA, ylim = NA, ...)
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
