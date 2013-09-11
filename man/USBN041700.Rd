\name{USBN041700}
\alias{USBN041700}
\docType{data}
\title{US Bonds and Notes Data on April 17, 2000}
\description{
Quotes used in the reference text as illustration of the construction of forward curves.
}
\usage{data(USBN041700)}
\format{
  A data frame with 99 observations on the following 9 variables.
  \describe{
    \item{\code{NAME}}{a factor with levels describing the type of Treasury, typically \code{US.TRSY.BONDS} or \code{US.TRSY.NOTES}}
    \item{\code{YEAR}}{a numeric vector giving the year of issue}
    \item{\code{EXPIRY}}{a numeric vector giving the expiration year}
    \item{\code{PRICE}}{a numeric vector giving the price}
    \item{\code{YIELD}}{a numeric vector giving the yield}
    \item{\code{RED.YIELD}}{a numeric vector}
    \item{\code{ACCR.INT.}}{a numeric vector giving the accrued interest}
    \item{\code{LIFE}}{a numeric vector the life (time to maturity)}
    \item{\code{VARIATION}}{a numeric vector}
  }
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
head(USBN041700
}
\keyword{datasets}
