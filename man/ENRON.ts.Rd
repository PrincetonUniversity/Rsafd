\name{ENRON.ts}
\alias{ENRON.ts}
\alias{EXXON.ts}
\docType{data}
\title{
DOW JONES Indexes for Energy Companies
}
\description{
Indexes computed from the values of the shares and the capitalizations of several energy companies, ENRON, EXXON, and whose behaviors over time are very similar to the behavior of the corresponding stock prices. The indexes are normalized to be equal to 100 on January 1, 1992, the first day for which they were computed.
}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
plot(ENRON.ts)	
plot(EXXON.ts)
}
\keyword{datasets}
