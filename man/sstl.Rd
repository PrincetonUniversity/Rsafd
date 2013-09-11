\name{sstl}
\alias{sstl}

\title{ Decomposition of a timeSeries Object into Trend, Seasonal and Remainder Components }

\description{
  Decomposition of a univariate timeSeries object based on the base function stl
}
\usage{
sstl(SERIES, FREQ = 365, TWIND = 0.75, TDEGREE = 0)
}
\arguments{
  \item{SERIES}{ Univariate timeSeries object to be decomposed}
  \item{FREQ}{Frequency of the seasonal component to be identified}
  \item{TWIND}{Factor multiplying the length of the series to obtain the length 
of the smoothing window used in LOESS}
  \item{TDEGREE}{Degree of locally-fitted polynomial in trend extraction. Should be 0 or 1.}
}

\details{
  The function sstl is a simple wrapper around the R function stl so that it 
would apply to objects of class timeSeries as defined in the library fCalendar
}
\value{
  Return a list with three components
  \item{trend}{a timeSeries of the deterministic trend}
  \item{sea}{the seasonal component identified as a periodic timeSeries object with period FREQ}
  \item{rem}{a timeSeries of the remainder}
}

\references{ R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
POS <- timeSequence(from="1959-01-01", to="1997-12-01", by="month")
co2.ts <- timeSeries(positions=POS,data=as.vector(co2),units="co2")

co2.stl <- sstl(co2.ts,FREQ=12,TDEGREE=1)
plot(co2.stl$trend)
plot(co2.stl$sea)
plot(co2.stl$rem)

Charlotte.stl <- sstl(Charlotte.ts[,3])
plot(Charlotte.stl$trend)
plot(Charlotte.stl$sea)
plot(Charlotte.stl$rem)
}
\keyword{ts}
\keyword{manip}
