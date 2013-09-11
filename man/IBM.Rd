\name{IBM}
\alias{IBM}
\docType{data}
\title{ High Frequency Data of IBM trades in June 1999}
\description{
The original data were processed to produce regularly spacedquotes, averaging to reduce the bid-ask spread so that we have only one single price quoteevery 30 seconds. There were 22 trading days on that month, so the data ends upbeing a 22 x 779 matrix. Used in Chapter 8 of the textbook "R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag".
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 8
}
\examples{
I <- 1plot(IBM[I,],type = "l", main="IBM High Frequency (30 s) on First Day of Trading in June 1999")
}
\keyword{datasets}
