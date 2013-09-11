\name{SPsep98.ts}
\alias{SPsep98.ts}
\docType{data}
\title{ High Frequency Data Time Series of S&P 500 Transaction Data: the September 1998 Futures Contract}
\description{
List giving transaction times and prices for the September 1998 S&P Futures contract
throughout the period the contract was actually traded.
}
\usage{data(SPsep98.ts)}
\format{
  The object of class timeSeries is saved as a list with two components
  \describe{
    \item{\code{pos}}{a vector of character strings giving the positions of the original timeSeries object}
    \item{\code{data}}{a numeric vector or matrix giving the data of the original timeSeries object}
  }
}
\references{
 R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
data(SPsep98.ts)
plot(SPsep98.ts)
}
\keyword{datasets}
