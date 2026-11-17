\name{plot}
\alias{plot}
\alias{plot,timeSeries-method}
\title{ Plot an Object of Class timeSeries }
\description{
  Produces a plot of a timeSeries object from library fCalendar
}
\usage{
\S4method{plot}{timeSeries}(x)
}
\arguments{
  \item{x}{ A timeSeries object, possibly multivariate }
}
\details{
Over write the plot method of the library fCalendar which is
limited to univariate timeSeries objects.
}
\value{
Used for its side effect, producing a plot of the timeSeries object.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\examples{
library(timeDate)
POS <- timeSequence(from="1959-01-01", to="1997-12-01", by="month")
co2.ts <- timeSeries(positions=POS,data=as.vector(co2), units="co2")
plot(co2.ts)

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
\keyword{ts}
\keyword{hplot}
