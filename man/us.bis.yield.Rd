\name{us.bis.yield}
\alias{us.bis.yield}
\docType{data}
\title{ US Treasury Yields}
\description{
Fixed (rolling) maturity yields for 1352 consecutive trading days following January 3, 1995.}
\usage{data(us.bis.yield)}
\format{
A data frame with 1352 observations on the following 11 variables.
The first variable gives (a proxy of) the short instantaneous rate, and the following variables 
give the yields, as computed by the Bank of International Settlements (BIS), for the maturities
1,2,3,4,5,5.5,6.5,7.5,8.5 and 9.5 months.}
\source{
  Bank of International Settlements
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 4.
}
\examples{
data(us.bis.yield)

dim(us.bis.yield)

us.bis.yield.pca <- princomp(us.bis.yield)
plot(us.bis.yield.pca)
title(sub="Proportions of the Variance Explained")
}
\keyword{datasets}
