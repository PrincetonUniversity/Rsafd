\name{PCS}
\alias{PCS}
\docType{data}
\title{ PCS Index Data}
\description{
Year to date aggregate amount of total damage reported in the United States
to the insurance industry. PCS stands for Property Claim Services.
}
\usage{data(PCS)}
\format{
  A data frame with 380 observations on the following 2 variables.
  \describe{
    \item{\code{Col1}}{a numeric vector coding time stamps for the catastrophic events}
    \item{\code{Col2}}{a numeric vector giving the total amount of the claims payed because of this event}
  }
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 2
}
\examples{
data(PCS)
}
\keyword{datasets}
