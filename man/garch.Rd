\name{garch}
\alias{garch}
\title{Fit a GARCH Model to a Time Series or a Vector}
\description{
Fit a GARCH model to a numeric vector or a timeSeries object.
}
\usage{garch(x, order = c(1, 1), series = NULL, control = garch.control(...), ...)}
\arguments{
  \item{x}{Input time series, a numeric vector or a timeSeries object}
  \item{order}{A vector of two integers c(p,q) giving the AR and MA orders of the model}
  \item{series}{For compatibility with \code{garch} function of \code{tseries} package}
  \item{control}{For compatibility with \code{garch} function of \code{tseries} package}
  \item{\dots}{For compatibility with \code{garch} function of \code{tseries} package}
}
\details{
Simple wrapper around the function \code{garch} of the package \code{tseries} to accommodate timeSeries objects.
}
\value{
Returns a list comprising
  \item{order}{The order of the model, same as input}
  \item{coef}{The estimated coefficients of the model}
  \item{.likelin}{The negative of the likelihood for the estimated coefficients}
  \item{resid}{Residuals of the model, as a timeSeries object if so is the input}
  \item{fitted.values}{Fitted values of the model, as a timeSeries object if so is the input}
  \item{sigma.t}{Instantaneous conditional standard deviation of the model, as a timeSeries object if so is the input}
  \item{vcov}{Variance/covariance matrix of the estimated coefficients}
  \item{name}{Name of the input time series}
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{
\code{\link{sim.garch}}
}
\examples{
LRET <- diff(log(EnronIndex.ts))
LRET.ar <- ar(LRET) 
LRES.g <- garch(LRET.ar$resid[!is.na(LRET.ar$resid)])
plot(LRES.g$sigma.t)
 }
\keyword{univar}
\keyword{distribution}
