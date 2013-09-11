\name{TRGSP}
\alias{TRGSP}
\alias{TRGSP2000}
\alias{TRGSP2}
\alias{TRGSP3}
\alias{TSTSP}
\alias{TSTSP2000}
\alias{TSTSP2}
\alias{TSTSP3}
\docType{data}
\title{ S&P 500 Option Data}
\description{
  Each TRGSPxx data set is the Training Sample for tests of the predictive powers 
of the various non-parametric regression methods discussed in R.A.Carmona, Statistical Analysis of
Financial Data in S-Plus. The performance of these regression methods are compared 
on corresponding test samples TSTSPxx.
}
\usage{data(TRGSP),data(TRGSP2000),data(TRGSP2),data(TRGSP3)
data(TSTSP),data(TSTSP2000),data(TSTSP2),data(TSTSP3)
}
\format{
  A data frame with the following 6 variables:
  \describe{
    \item{\code{SP}}{a numeric vector}
    \item{\code{KK}}{a numeric vector}
    \item{\code{TAU}}{a numeric vector}
    \item{\code{RR}}{a numeric vector}
    \item{\code{ISIGMA}}{a numeric vector}
    \item{\code{CALL}}{a numeric vector}
  }
}

\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag.
}
\examples{
data(TRGSP)
}
\keyword{datasets}
