\name{DukeLRet}
\alias{DukeLRet}
\alias{EnronLRet}
\alias{UtilLRet}
\docType{data}
\title{Duke daily log-returns}
\description{
  Numeric vector giving the daily log-returns for the DUKE.index, ENRON.index, Brazilian and Columbian coffee contracts used
  in the textbook "R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag"
}
\usage{
data(DukeLRet)
data(EnronLRet)
data(UtilLRet)
}
\format{
Numeric vectors of length 259. \code{UtilLRet} is equal to
\code{diff(log(UTIL.index$UTILITY.index))}. The coffee log-returns
\code{BCofLRet} and \code{CCofLRet} are documented in \code{\link{BCofLRet}}.
}
\examples{
plot(DukeLRet,EnronLRet)
plot(BCofLRet,CCofLRet)
}
\keyword{datasets}
\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 2 & Chapter 4}
