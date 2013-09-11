\name{DSPLRet}
\alias{DSPLRet}
\docType{data}
\title{
Numeric vector of daily log returns of the S&P500 index
}
\description{
A numerical vector used in chapters 1 & 2 of the textbook "R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag"
for illustration purposes}

\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 1 & Chapter 2}

\examples{
hist(DSPLRet,col="blue",freq=F,ylim=c(0,30))DENS <- density(DSPLRet, adjust=12)
}
\keyword{datasets}
