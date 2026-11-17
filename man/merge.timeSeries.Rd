\name{merge}
\alias{merge}
\alias{merge.timeSeries}
\alias{merge,timeSeries,timeSeries-method}
\title{ Merge two objects of class timeSeries into a single one }
\description{
Given two objects of class timeSeries, compute the intersection of the set of positions, and 
use this intersection as the set of positions of the output timeSeries whose data is the union
of the columns of the input timeSeries restricted to these positions.
}
\usage{
\method{merge}{timeSeries}(x, y, INTERSECT = FALSE, units = NULL, ...)
\S4method{merge}{timeSeries,timeSeries}(x, y)
}
\arguments{
  \item{x}{ timeSeries object }
  \item{y}{ timeSeries object }
  \item{INTERSECT}{ Logical. If \code{TRUE}, the positions of the output are the
intersection of the positions of \code{x} and \code{y} as described above. If \code{FALSE}
(the default, and the behavior of the \code{merge} method), the positions of the output
are the union of the positions of \code{x} and \code{y}, and missing values are set to \code{NA}. }
  \item{units}{ Desired column names }
  \item{\dots}{ additional arguments, currently not used }
}
\details{
  Works only for two timeSeries
}
\value{
A single timeSeries object as described above.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\seealso{ \code{\link{timeSeries}}, \code{\link[=plot]{plot.timeSeries}},\code{\link{dump.timeSeries}}  }
\examples{
data(DesMoines.ts)
DM.ts <- DesMoines.ts
data(Covington.ts)
CVG.ts <- Covington.ts
data(Portland.ts)
PDX.ts <- Portland.ts
TEMPS.ts <- merge(DM.ts,CVG.ts)
TEMPS.ts <- merge(TEMPS.ts,PDX.ts)
plot(TEMPS.ts[1:1000,])
}
\keyword{ ts }
