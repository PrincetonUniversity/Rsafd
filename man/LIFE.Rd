\name{LIFE}
\alias{LIFE}
\alias{INT.YIELD}
\docType{data}
\title{Time to Maturity and Yield of US Bonds and Notes on April 17, 2000}
\description{
Numeric vectors extracted from the data set \code{\link{USBN041700}}: \code{LIFE} is the column \code{LIFE} (time to maturity) and \code{INT.YIELD} is the column \code{YIELD}.
}
\usage{
data(LIFE)
data(INT.YIELD)
}
\format{
Numeric vectors of length 99.
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\examples{
plot(LIFE, INT.YIELD)
}
\keyword{datasets}
