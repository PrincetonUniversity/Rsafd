\name{FRWRD}
\alias{FRWRD}
\docType{data}
\title{Natural Gas Forward Contract Prices}
\description{
Example of a Natural Gas forward curve used in 
R.A. Carmona, Statistical Analysis of Financial Data in R, 2013, Springer Verlag.
}
\usage{FRWRD}
\format{
Numeric vector of length 36 giving the prices of forward contracts of Henry Hub  Natural Gas (NG) with times to maturity 1,2,3, ..... , 36 months.}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag.
}
\examples{
plot(1:36,FRWRD,main="Polynomial Gas Forward Curves")
lines(1:36,fitted(lm(FRWRD~poly(1:36,3))))
lines(1:36,fitted(lm(FRWRD~poly(1:36,6))),lty=3)
lines(1:36,fitted(lm(FRWRD~poly(1:36,8))),lty=6)
}
\keyword{datasets}
