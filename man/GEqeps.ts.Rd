\name{GEqeps.ts}
\alias{GEqeps.ts}
\alias{IBMqeps.ts}
\alias{PEPqeps.ts}
\docType{data}
\title{ Quarterly Earnings}
\description{Quarterly earnings of a company for the purpose of the filtering analysis done in Chapter 6 of R.Carmona, Statistical Analysis of Financial Data in R.
}
\arguments{
Univariate timeSeries objects with a slot pos and a slot data 
\item{\code{pos}}{Date of the last day of the quarter for which the earnings are reported}
\item{\code{data}}{Quarterly earnings}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\examples{
plot(GEqeps.ts)
}
\keyword{datasets}
