\name{c}
\alias{c}
\alias{concat}
\alias{concat,timeSeries,timeSeries-method}
\title{ Concatenate Two timeSeries Objects }
\description{
  One (simple) way to put together two timeSeries with the same variables.
}
\usage{
concat(x, y)
}
\arguments{
  \item{x}{ timeSeries object }
  \item{y}{ timeSeries object }
}
\details{
  Check that the dimensions of \code{x} and \code{y} are the same.
The current implementation assumes that the positions
of \code{x} and \code{y} do not overlap.
}
\value{
  A timeSeries object whose positions are given by the concatenation of the 
positions of \code{x} and \code{y}, and whose data is obtained by row-binding
the data of \code{x} and \code{y}. 
}
\author{ Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{  \code{\link{timeSeries}}}
\examples{
library(timeDate)
DAY1 <- timeDate("08/15/1998")
DAY2 <- timeDate("08/20/1998")
DAY3 <- timeDate("08/25/1998")
TMP1 <- SPsep98.ts[seriesPositions(SPsep98.ts) >= DAY1 & seriesPositions(SPsep98.ts) < DAY2]
TMP2 <- SPsep98.ts[seriesPositions(SPsep98.ts) >= DAY2 & seriesPositions(SPsep98.ts) < DAY3]
SDSP.ts <- concat(TMP1,TMP2)
plot(SDSP.ts)
}
\keyword{ts}
