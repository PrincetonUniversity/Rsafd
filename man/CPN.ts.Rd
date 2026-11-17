\name{CPN.ts}
\alias{CPN.ts}
\alias{CPN}
\alias{CPNLRet}
\docType{data}
\title{timeSeries of Calpine daily close}
\description{timeSeries object}
\usage{
data(CPN.ts)
data(CPN)
data(CPNLRet)
}
\examples{
plot(CPN.ts) 
}
\keyword{datasets}
\format{
Objects of class timeSeries with a slot positions and a slot data. 
\describe{
\item{\code{positions}}{dates in timeDate format, 
  (e.g. from 1960-01-01 to 199-03-31)} 
\item{\code{data}}{a numeric vector or matrix containing the values observed at these dates.}
}
\code{CPN.ts} has 2155 observations from 1996-09-20 to 2005-04-15.
\code{CPN} is a 2155 x 1 numeric matrix containing the data of \code{CPN.ts}, and
\code{CPNLRet} is the 2154 x 1 numeric matrix \code{diff(log(CPN))} of the log-returns.
}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
