\name{as.ts.timeSeries}
\alias{as.ts.timeSeries}
\alias{acf.timeSeries}
\alias{arima.timeSeries}
\title{Standard Time Series Functions Applied to timeSeries Objects}
\description{
\code{as.ts.timeSeries} converts a \code{timeSeries} object into an object of class
\code{ts}. \code{acf.timeSeries} computes and plots the auto-correlation function of
a \code{timeSeries} object. \code{arima.timeSeries} fits an ARIMA model to the data of
a \code{timeSeries} object.
}
\usage{
\method{as.ts}{timeSeries}(x, ...)
acf.timeSeries(x, y)
arima.timeSeries(x, order)
}
\arguments{
  \item{x}{object of class \code{timeSeries}}
  \item{y}{character string giving the type of auto-correlation function, passed as
    argument \code{type} to \code{\link[stats]{acf}}: \code{"correlation"},
    \code{"covariance"} or \code{"partial"}}
  \item{\dots}{additional arguments, currently not used}
  \item{order}{vector of three integers giving the order of the ARIMA model, passed to
    \code{\link[stats]{arima}}}
}
\value{
\code{as.ts.timeSeries} returns an object of class \code{ts} with the data of \code{x}
indexed by \code{1, 2, ...}. \code{acf.timeSeries} returns the object computed by
\code{\link[stats]{acf}} and plots it. \code{arima.timeSeries} returns the object
computed by \code{\link[stats]{arima}}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{\code{\link{timeSeries}}}
\examples{
LRET <- diff(log(CPN.ts))
X <- as.ts.timeSeries(LRET)
acf.timeSeries(LRET, "correlation")
arima.timeSeries(LRET, order = c(1, 0, 0))
}
\keyword{ts}
