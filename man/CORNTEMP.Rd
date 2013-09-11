\name{CORNTEMP}
\alias{CORNTEMP}
\docType{data}
\title{Corn and Temperature Data}
\description{
Monthly price of corn in Iowa and temperature data in Des Moines to illustrate the dependence of the two via regression, and justify the use of temperature options to hedge commodity price risk. These data are used in a problem at the end of Chapter 4 of the textbook. 
}
\format{
  A data frame with 468 observations on the following 3 variables.
  \describe{
    \item{\code{MCorn}}{Monthly price of a bushel of corn in Iowa}
    \item{\code{MDMHDD}}{Monthly average of the number of heating degree days in Des Moines as measured at the meteorological station of the airport.} 
    \item{\code{MDMHDD}}{Monthly average of the number of cooling degree days in Des Moines as measured at the meteorological station of the airport.} 
  }
}
\details{
    MCorn     MDMHDD    MDMCDD
175  0.90 43.5781250 0.0000000
176  0.88 46.1551724 0.0000000
177  0.90 42.5161290 0.0000000
178  0.96 14.0000000 0.9666667
179  1.00  4.6129032 1.3387097
180  1.03  0.4833333 4.2500000
}
\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 5
}
\keyword{datasets}
