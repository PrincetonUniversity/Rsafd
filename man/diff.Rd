\name{diff}
\alias{diff}
\alias{diff,timeSeries-method}
\title{Differentiation of timeSeries objects }
\description{
Compute differences of a timeSeries and some of its lagged versions.
}
\usage{
\S4method{diff}{timeSeries}(x, lag = 1, diff = 1, ...)
}
\arguments{
  \item{x}{ \code{timeSeries} object }
  \item{lag}{ Integer giving the number of lags by which the \code{timeSeries} object should be lagged. }
  \item{diff}{ Integer giving the number of iterated differences to compute }
  \item{\dots}{ additional arguments, currently not used }
}
\value{
if Y <- diff(X,lag=k,diff=1), then Y[t]=X[t]-X[t-k] if k>=0, and Y[t]=X[t-k]-X[t] if k<0. If diff>1, then
the procedure is iterated diff times.}
\note{Essentially the same as the function \code{diff.timeSeries} of the package \code{fCalendar}
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{
    \code{\link[=plot]{plot.timeSeries}},\code{\link{timeSeries}}, \code{\link{dump.timeSeries}}, \code{\link{merge.timeSeries}}
}
\examples{
Y1 <- diff(Charlotte.ts)
Y5 <- diff(Charlotte.ts,lag=5)
plot(Y1[1:1000])
}
\keyword{manip}
