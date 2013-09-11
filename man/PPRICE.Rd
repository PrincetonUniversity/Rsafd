\name{PPRICE}
\alias{PPRICE}
\docType{data}
\title{Daily Energy Prices and Temperatures}
\description{
Data set of daily gas spot prices, temperatures, spot and forward electricity prices. 
Data set used in one of the problems of Chapter 5 of  
"R.A. Carmona, Statistical Analysis of Financial Data in R, 2013, Springer Verlag".
}
\usage{PPRICE}
\format{
Data frame with 680 rows and 4 columns.
The columns are named \code{GasSpot}, \code{SDTemp}, \code{PPower}, and  \code{FPower}. 
  \describe{
    \item{\code{GasSpot}}{ gives the spot price of natural gas on that date}
    \item{\code{SDTemp}}{gives the average temperature in San Diego
over the 31 days preceding the date in question}
    \item{\code{PPower}}{gives the average over the $5$ days preceding the date in question,
of the spot price of firm on peak electric power at the Palo Verde station}
    \item{\code{FPower}}{gives the average over the two weeks following the date in question,
of the spot price of firm on peak electric power at the Palo Verde station}
}
}
\references{
 R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
head(PPRICE)}
\keyword{datasets}
