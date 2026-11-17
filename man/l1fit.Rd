\name{l1fit}
\alias{l1fit}

\title{Least Absolute Deviations Regression}
\description{
Share the format of lsfit but perform a L1 regression instead of a least squares regression.
}
\usage{
l1fit(x, y, intercept = TRUE)
}
\arguments{
  \item{x}{Vector (simple regression) or matrix (multiple regression) of explanatory variables}
  \item{y}{Vector of observations of response variable. Length of \code{y} should be first dimension of \code{x}}
  \item{intercept}{Boolean (TRUE by default) which should be set to FALSE if regression through the origin is required. }
}
\details{
As Li regression is a particular case of quantile regression,
this function is just a wrapper to rq() from the library robustbase. Use that instead!
}
\value{
See the help file of rq
}
\references{library robustbase}

\keyword{regression}
