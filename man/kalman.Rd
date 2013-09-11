\name{kalman}
\alias{kalman}

\title{ Kalman Filter for a Linear Partially Observed System }
\description{
  Use the Kalman filter to compute the one-day ahead prediction
and its standard error in a linear state space system with partial
observations.
}
\usage{
kalman(FF, SigV, GG, SigW, Xhat, Omega, Y)
}

\arguments{
  \item{FF}{ d[X] by d[X] "state matrix" giving the time evolution of the partially observed state }
  \item{SigV}{ d[X] by d[X] variance/covariance matrix of the noise term of the state evolution equation }
  \item{GG}{ d[Y] by d[X] "observation matrix" entering the observation equation}
  \item{SigW}{ d[Y] by d[Y] variance/covariance matrix of the noise term entering the observation equation }
  \item{Xhat}{ d[X] dimensional vector representing the current value of the one step ahead state vector prediction }
  \item{Omega}{ d[X] by d[X] matrix representing the current value of the variance/covariance matrix of Xhat}
  \item{Y}{ d[Y] dimensional vector representing the new observation }
}

\details{
  Update the one step-ahead prediction and the estimate of its variance-covariance matrix using the
Kalman recursive filter.
}

\value{
  Returns a list with two components
  \item{xpred}{Vector of the updates of the components of Xhat}
  \item{error}{Updated estimate of the variance/covariance matrix of Xhat}
}

\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}

\examples{}

\keyword{multivariate}
