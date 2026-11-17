\name{DUKE.ts}
\alias{DUKE.ts}
\alias{DYNEGY.ts}
\alias{TXU.ts}
\alias{SO.ts}
\docType{data}
\usage{
data(DUKE.ts)
data(DYNEGY.ts)
data(TXU.ts)
data(SO.ts)
}
\format{Objects of class \code{timeSeries} with one column each.}
\title{
Time series of the daily closing prices of Energy Companies 
}
\description{
timeSeries objects giving the daily closing prices of DUKE Energy, Dynegy, TXU Energy and SOUTHERN Energy}
\details{
Daily values of the closing prices of DUKE Energy between April 6, 1983 and June 17, 2003

Daily values of the closing prices of Dynegy between January 3, 1994 and December 5, 2001

Daily values of the closing prices of SOUTHERN between December 31, 1981 and June 17, 2003

Daily values of the closing prices of TXU Energy between June 10, 1983 and June 17, 2003
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
plot(DUKE.ts) 
}
\keyword{datasets}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
