\name{CO2TRG}
\alias{CO2TRG}
\alias{CO2TST}
\docType{data}
\title{CO2 Option Data}
\description{
Data frames of option data used in the reference text. \code{CO2TRG} and \code{CO2TST} have the same structure.
}
\usage{
data(CO2TRG)
data(CO2TST)
}
\format{
Data frames with 880 (\code{CO2TRG}) and 687 (\code{CO2TST}) observations on the following 5 variables.
\describe{
  \item{\code{volume}}{a numeric vector}
  \item{\code{eua}}{a numeric vector}
  \item{\code{tau}}{a numeric vector}
  \item{\code{strike}}{a numeric vector}
  \item{\code{call}}{a numeric vector}
}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\examples{
head(CO2TRG)
summary(CO2TST)
}
\keyword{datasets}
