\name{y4nls1}
\alias{y4nls1}
\alias{y4nls2}
\docType{data}
\title{ Crude Oil Forward Data}
\description{
  Crude Oil Forward Prices on April 16, 1998 (y4nls1) and MArch 23, 1988 (y4nls2) used as examples of data treated by non-linear parametric regression in the reference text.
}
\format{
Numeric vectors}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
par(mfrow=c(1,2))
x <- 1:length(y4nls1)
plot(x,y4nls1,main="4/16/1998")
x <- 1:length(y4nls2)
plot(x,y4nls2,main="3/23/1998")
par(mfrow=c(1,1))
}
\keyword{datasets}
