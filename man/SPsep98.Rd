\name{SPsep98}
\alias{SPsep98}
\docType{data}
\title{ High Frequency Data: S&P 500 Data, September 1998 Futures Contract}
\description{
Data frame giving transaction prices of the September 1998 S&P Futures contract
throughout the period the contract was traded.
}
\usage{data(SPsep98)}
\format{
  A data frame with 248155 observations on the following 3 variables.
  \describe{
    \item{\code{date}}{a factor with levels \code{19971008}, .... ,code{19980917} giving the year, month and day of the transaction}
    \item{\code{time}}{a factor with levels \code{08:30:00}, .... ,\code{15:45:00} giving the time of the day the transaction took place}
    \item{\code{close}}{a numeric vector giving the price at which the transaction took place}
  }
}
\references{
 R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
data(SPsep98)
}
\keyword{datasets}
