\name{DukeLRet}
\alias{BCofLRet}
\alias{CCofLRet}
\alias{EnronLRet}
\docType{data}
\title{Duke daily log-returns}
\description{
  Numeric vector giving the daily log-returns for the DUKE.index, ENRON.index, Brazilian and Columbian coffee contracts used
  in the textbook "R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag"
}
\usage{head(DukeLRet)}
\format{
}
\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 2 & Chapter 4}
\examples{
plot(DukeLRet,EnronLRet)
plot(BCofLRet,CCofLRet)
}
\keyword{datasets}
