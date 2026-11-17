\name{bitcoin}
\alias{bitcoin}
\title{Daily Bitcoin price between 8/25/2014 and 1/6/2022}
\description{
Example of daily prices use for the purpose of illustration of heavy tail distribution fitting.
}
\usage{
bitcoin}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
plot(1:2684,bitcoin, type="l")
}
\keyword{financial}
\format{
Numeric vector of length 2684 giving the prices of Bitcoin over the period ranging from 8/25/2014 to 1/6/2022. 8 days were removed because of the absence of trading on those days.}
