\name{MSP}
\alias{MSP}
\docType{data}
\title{ Minute by Minute S&P 500 Quotes}
\description{
This numeric vector contains the minute by minute quotes of the S&P 500 index on September 10, 1998.}
\usage{data(MSP)}
\format{
  The format is:
 num [1:404] 1080 1078 1078 1080 1080 ...
}
\details{
Data needed for a Problem in Chapter 1 of the book below.
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
head(MSP)
plot(MSP,type="l") 
}
\keyword{datasets}
