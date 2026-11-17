\name{stocks}
\alias{stocks}
\title{Daily stock closing prices}
\description{
Daily closing prices for a set of 36 stocks chosen for an application of a graphical Gaussian model with the function \code{glasso}. The dates start on May 18, 2012 to include Facebook.}
\usage{data(stocks)}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
head(stocks)
}
\keyword{datasets}
\format{
  A numerical matrix with 2170 observations (name after the trading days) and 36 columns (one for each stock). The column names are the tickers of the stocks and the row names are the dates
}
