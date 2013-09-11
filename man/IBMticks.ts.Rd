\name{IBMticks.ts}
\alias{IBMticks.ts}
\docType{data}
\title{
Tick by tick (transaction) data for IBM on May 10 - 14, 2004
}
\description{
High frequency timeSeries object of 27,343 transactions. The positions (time stamps of the transactions) are given in the format Year/month/day Hour:Minute:Second.

The first entry is the price at which the transaction took place and the second is the volume (number of shares exchanged) of the transaction.}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
plot(IBMticks.ts) 
}
\keyword{datasets}
