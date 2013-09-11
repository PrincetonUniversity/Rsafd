\name{swap}
\alias{swap}
\docType{data}
\title{Swap Rates}
\description{
Swap rates collected from May 1998 to March 2000.}
\usage{data(swap)}
\format{
  A data frame with 496 observations on 15 variables
  giving the swap rates for the constant (rolling) maturities 1,2,3,4,5,6,7,8,9,10,12,15,20,25, and 30 years.
}
\source{
Data Stream
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
data(swap)
swap.pca <- princomp(swap)
plot(swap.pca)
title("Proportions of the Variance Explained by the Components")

YEARS <- c(1,2,3,4,5,6,7,8,9,10,12,15,20,25,30)
par(mfrow=c(2,2))
plot(YEARS,swap.pca$loadings[,1],ylim=c(-.6,.6))
lines(YEARS,swap.pca$loadings[,1])
plot(YEARS,swap.pca$loadings[,2],ylim=c(-.6,.6))
lines(YEARS,swap.pca$loadings[,2])
plot(YEARS,swap.pca$loadings[,3],ylim=c(-.6,.6))
lines(YEARS,swap.pca$loadings[,3])
plot(YEARS,swap.pca$loadings[,4],ylim=c(-.6,.6))
lines(YEARS,swap.pca$loadings[,4])
par(mfrow=c(1,1))
title("First Four Loadings of the Swap Rates")
}
\keyword{datasets}
