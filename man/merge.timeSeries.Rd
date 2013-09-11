\name{merge}
\alias{merge}

\title{ Merge two objects of class timeSeries into a single one }
\description{
Given two objects of class timeSeries, compute the intersection of the set of positions, and 
use this intersection as the set of positions of the output timeSeries whose data is the union
of the columns of the input timeSeries restricted to these positions.
}
\usage{
merge.timeSeries(x, y, units = NULL, ...)
}
\arguments{
  \item{x}{ timeSeries object }
  \item{y}{ timeSeries object }
  \item{units}{ Desired column names }
}
\details{
  Works only for two timeSeries
}
\value{
A single timeSeries object as described above.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\seealso{ \code{\link{timeSeries}}, \code{\link{plot.timeSeries}},\code{\link{dump.timeSeries}}  }
\examples{
data(DesMoines.ts)
DM.ts <- timeSeries(positions=DesMoines.ts$pos,data=DesMoines.ts$data,units="DM")
data(Covington.ts)
CVG.ts <- timeSeries(positions=Covington.ts$pos,data=Covington.ts$data,units="CVG")
data(Portland.ts)
PDX.ts <- timeSeries(positions=Portland.ts$pos,data=Portland.ts$data,units="PDX")
TEMPS.ts <- merge(DM.ts,CVG.ts)
TEMPS.ts <- merge(TEMPS.ts,PDX.ts)
plot(TEMPS.ts[1:1000,])
}
\keyword{ ts }
