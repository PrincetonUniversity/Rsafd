\name{MONTHLY}
\alias{MONTHLY}
\docType{data}
\title{
Monthly values of the short and long interest rates and S&P 500 index
}
\description{
Data frame with 163 rows and 3 columns of monthly values used in Section 8.5 and in two problems in Chapter 8 of the textbook "Statistical Analysis of Financial Data in R".
The row names are the dates of the measurements. They range between May 1, 1986 and November 1, 1999.
}
\usage{data(MONTHLY)}
\format{
  A data frame with 163 observations on the following 3 variables.
  \describe{
    \item{\code{Short}}{a numeric vector of the monthly values of one year Treasury Bills}
    \item{\code{Long}}{a numeric vector of the monthly values of the 30 year Treasury Bonds}
    \item{\code{SP500}}{a numeric vector of the monthly values of the S&P 500 index}
  }
}

\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 5
}
\examples{
head(MONTHLY)
}
\keyword{datasets}
