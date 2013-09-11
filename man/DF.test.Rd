\name{DF.test}
\alias{DF.test}

\title{
Dickey-Fuller Unit Root Test
}
\description{Simple wrapper around the Augmented Dickey-Fuller test proposed in the package \code{tseries}, written	to accommodate timeSeries objects
}
\usage{
DF.test(x, alternative = c("stationary", "explosive"), k)
}
\arguments{
  \item{x}{Numeric vector or object of class timeSeries}
  \item{alternative}{Character string "stationary" or "explosive"}
  \item{k}{Integer, for compatibility with \code{tseries}}
}
\value{
Returns a list comprising
  \item{statistic}{The numerical value of the statistics}
  \item{p.value}{p-value of the test}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\keyword{time series}

