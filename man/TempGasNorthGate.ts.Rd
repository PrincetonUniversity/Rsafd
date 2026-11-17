\name{TempGasNorthGate.ts}
\alias{TempGasNorthGate.ts}
\docType{data}
\title{Temperature and Natural Gas Daily Data}
\description{
Object of class \code{timeSeries} with two columns, \code{GateTemp} and \code{NGas}, used in the reference text.
}
\usage{
data(TempGasNorthGate.ts)
}
\format{
An object of class \code{timeSeries} with 10593 daily observations of 2 variables (\code{GateTemp} and \code{NGas}) from 1980-01-01 to 2008-12-31.
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\examples{
plot(TempGasNorthGate.ts[1:1000, ])
}
\keyword{datasets}
