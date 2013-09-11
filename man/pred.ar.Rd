\name{pred.ar}
\alias{pred.ar}
\title{ Compute Predictions for an AR Model }
\description{
 
}
\usage{
pred.ar(series, ar.est, ahead = 1)
}
\arguments{
  \item{series}{ Time Series to which the AR model was fitted }
  \item{ar.est}{ Estimated AR model }
  \item{ahead}{ Integer giving the number of lags ahead to be predicted }
}
\details{
 Wrapper for the function \code{ar.pred} given in the ar help file of S-Plus
}
\value{
  Numerical vector of length \code{ahead} giving the desired predictions.
}
\references{  R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag }
\author{ Rene Carmona }
\seealso{ See Also \code{\link{ar}} }
\keyword{ models }
