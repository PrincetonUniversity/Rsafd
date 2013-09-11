\name{TRGSS}
\alias{TRGSS}
\alias{TSTSS}
\docType{data}
\title{(Simulated) Data Sets used in the Regression Problem on the Pricing of Spread Options}
\description{
TYRGSS and TSTSS are numeric matrices with 8 columns and 831 rows for TRGSS and 300 rows for TSTSS. The first two columns give the daily prices of forward contracts of electricity and natural gas respectively, the next six columns giving prices of zero-strike spread options with time to maturity 2,4,6,8,10 and 12 months.
TRGSS is the training sample, TSTSS the testing sample. 
}
\usage{data(MID1)}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
dim(TRGSS)
# [1] 831 8
}
\keyword{datasets}
