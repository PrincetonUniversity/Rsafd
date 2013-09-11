\name{plot}
\alias{plot}

\title{ Plot an Object of Class timeSeries }
\description{
  Produces a plot of a timeSeries object from library fCalendar
}
\usage{
plot(x)
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
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\examples{
POS <- timeSequence(from="1959-01-01", to="1997-12-01", by="month")
co2.ts <- timeSeries(positions=POS,data=as.vector(co2), units="co2")
plot(co2.ts)

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

\keyword{ts}
\keyword{hplot}
