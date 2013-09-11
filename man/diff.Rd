\name{diff}
\alias{diff}
\title{Differentiation of timeSeries objects }
\description{
Compute differences of a timeSeries and some of its lagged versions.
}
\usage{
diff(x, lag = 1, diff = 1, trim = TRUE, pad = NA, ...)
}
\arguments{
  \item{x}{ \code{timeSeries} object }
  \item{lag}{ Integer giving the number of lags by which the \code{timeSeries} object should be lagged. }
  \item{diff}{ Integer giving the number of iterated differences to compute }
  \item{trim}{ Boolean to return a \code{timeSeries} objects with the same positions as \code{x}
when \code{FALSE}, or if \code{TRUE}, a \code{timeSeries} object without the positions for which would have carried an NA. }
  \item{pad}{ What should be added in lieu of entries which cannot be computed }
}
\details{
}
\value{
if Y <- diff(X,lag=k,diff=1), then Y[t]=X[t]-X[t-k] if k>=0, and Y[t]=X[t-k]-X[t] if k<0. If diff>1, then
the procedure is iterated diff times.}
\note{Essentially the same as the function \code{diff.timeSeries} of the package \code{fCalendar}
}
\seealso{
    \code{\link{plot.timeSeries}},\code{\link{timeSeries}}, \code{\link{dump.timeSeries}}, \code{\link{merge.timeSeries}}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
Y1 <- diff(Charlotte.ts)
Y5 <- diff(Charlotte.ts,lag=5)
plot(Y1[1:1000])
}
\keyword{manip}
