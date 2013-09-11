\name{c}
\alias{c}

\title{ Concatenate Two timeSeries Objects }
\description{
  One (simple) way to put together two timeSeries with the same variables.
}
\usage{
c(X, Y)
}

\arguments{
  \item{X}{ timeSeries object }
  \item{Y}{ timeSeries object }
}
\details{
  Check that the dimensions of \code{X} and \code{Y} are the same.
The current implementation assumes that the positions
of \code{X} and \code{Y} do not overlap.
}
\value{
  A timeSeries object whose positions are given by the concatenation of the 
positions of \code{X} and \code{Y}, and whose data is obtained by row-binding
the data of \code{X} and \code{Y}. 
}

\author{ Rene Carmona, \email{rcarmona@princeton.edu}}

\seealso{  \code{\link{timeSeries}}}
\examples{
DAY1 <- timeDate("08/15/1998")
DAY2 <- timeDate("08/20/1998")
DAY3 <- timeDate("08/25/1998")
TMP1 <- SPsep98.ts[seriesPositions(SPsep98.ts) >= DAY1 & seriesPositions(SPsep98.ts) < DAY2]
TMP2 <- SPsep98.ts[seriesPositions(SPsep98.ts) >= DAY2 & seriesPositions(SPsep98.ts) < DAY3]
SDSP.ts <- c(TMP1,TMP2)
plot(SDSP.ts)
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\keyword{ts}
