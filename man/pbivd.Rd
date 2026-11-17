\name{pbivd}
\alias{pbivd}
\alias{dbivd}
\alias{pbivd,bivd-method}
\alias{pbivd,gpd.bivd-method}
\alias{dbivd,bivd-method}
\alias{dbivd,gpd.bivd-method}
\title{Cumulative Distribution Function and Density of a Bivariate Distribution}
\description{
Compute the values of the cumulative distribution function (\code{pbivd}) and of the
density (\code{dbivd}) of a bivariate distribution given by a copula and two marginal
distributions.
}
\usage{
pbivd(dist, x, y)
dbivd(dist, x, y)
}
\arguments{
  \item{dist}{object of class \code{bivd} (as created by \code{\link{bivd}}) or \code{gpd.bivd}}
  \item{x}{numeric vector of values of the first variable}
  \item{y}{numeric vector of values of the second variable}
}
\details{
The marginal cumulative distribution functions are evaluated at \code{x} and \code{y},
and the result is passed to \code{\link{pcopula}} (respectively \code{\link{dcopula}},
multiplied by the marginal densities) for the copula of \code{dist}.
}
\value{
Numeric vector of the values of the cdf (\code{pbivd}) or the density (\code{dbivd}).
}
\section{Methods}{
\describe{
\item{\code{signature(dist = "bivd")}}{marginals given by the names and parameters of univariate distributions}
\item{\code{signature(dist = "gpd.bivd")}}{marginals given by objects of class \code{gpd}}
}}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\code{\link{bivd}}, \code{\link{rbivd}}, \code{\link{contour.dbivd}}, \code{\link{persp.dbivd}}}
\examples{
BIV <- bivd(gumbel.copula(1.4), "norm", "t", c(3,4), 3)
pbivd(BIV, c(2, 3), c(0, 1))
dbivd(BIV, c(2, 3), c(0, 1))
}
\keyword{distribution}
\keyword{multivariate}
