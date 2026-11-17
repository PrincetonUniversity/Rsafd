\name{gpd.1p}
\alias{gpd.1p}
\alias{gpd.1q}
\alias{gpd.2p}
\alias{gpd.2q}
\title{Cumulative Distribution Function and Quantiles of Fitted GPD Tails}
\description{
Compute the cumulative distribution function (\code{gpd.1p}, \code{gpd.2p}) and the
quantile function (\code{gpd.1q}, \code{gpd.2q}) of a distribution estimated by
\code{\link{fit.gpd}}, with one upper tail (\code{gpd.1p}, \code{gpd.1q}) or two
tails (\code{gpd.2p}, \code{gpd.2q}) fitted by generalized Pareto distributions.
}
\usage{
gpd.1p(obj, x, linear = TRUE)
gpd.1q(obj, p, linear = TRUE)
gpd.2p(obj, x, linear = TRUE)
gpd.2q(obj, p, linear = TRUE)
}
\arguments{
  \item{obj}{object of class \code{gpd.upper.tail} (for \code{gpd.1p}, \code{gpd.1q})
    or \code{gpd.two.tails} (for \code{gpd.2p}, \code{gpd.2q}), as returned by
    \code{\link{fit.gpd}}}
  \item{x}{numeric vector (or \code{timeSeries}) of values at which the cdf is computed}
  \item{p}{numeric vector of probabilities at which the quantiles are computed}
  \item{linear}{logical, currently not used}
}
\details{
Between the thresholds, the empirical distribution of the data is used. Beyond the
thresholds, the generalized Pareto distributions with the estimated parameters are used.
}
\value{
Numeric vector of the same length as \code{x} (respectively \code{p}). The quantile
functions return \code{NA} for probabilities outside [0,1].
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 2}
\seealso{\code{\link{pgpd}}, \code{\link{qgpd}}, \code{\link{fit.gpd}}}
\examples{
SHAPE.XI <- TRUE
X <- rpareto(1024, xi = 0.3)
X.est <- fit.gpd(X, tail = "upper", plot = FALSE)
gpd.1p(X.est, c(1, 2, 5))
gpd.1q(X.est, c(0.5, 0.9, 0.99))
}
\keyword{distribution}
