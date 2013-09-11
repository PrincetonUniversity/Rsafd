\name{hills}
\alias{hills}
\docType{data}
\title{Scottish Hill Races Record Times}
\description{
Standard data set (see source below) -- already used in the package MASS -- used in the in a problem of the reference text. 
}
\usage{data(hills)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{dist}}{length of the race}
    \item{\code{climb}}{incline in feet}
    \item{\code{time}}{record time in minutes}
  }
}
\source{
A.C. Atkinson (1986) Comment: Aspects of diagnostic regression analysis. Statistical Science 1, 397–402.
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\examples{
head(hills)
}
\keyword{datasets}
