\name{CPN.ts}
\alias{CPN.ts}
\title{timeSeries of Calpine daily close}
\description{timeSeries object}
\usage{data(CPN.ts)}
\arguments{    
\item{\code{positions}}{dates in timeDate format, 
  (e.g. from 1960-01-01 to 199-03-31)} 
\item{\code{data}}{a numeric vector or matrix containing the values observed at these dates.}
}
\examples{
plot(CPN.ts) 
}
\keyword{datasets}
\format{
Objects of class timeSeries with a slot positions and a slot data. 
}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
