\name{fit.copula}
\alias{fit.copula}

\title{Maximum Likelihood Fit of a Copula }
\description{
Given a parametric family of copulae, compute the maximum likelihood estimate of the parameter}

\usage{
fit.copula(data, family = "normal", plotPicture = T, init.est = NA, epsilon = 1e-05, ...)
}
\arguments{
  \item{data}{ n x 2 matrix giving the n observations }
  \item{family}{ ~~Describe \code{family} here~~ }
  \item{plotPicture}{ Produces a contour plot illustrating the fit }
  \item{init.est}{ Vector of the initial values of the parameters of the copula. Used to start the maximization of the likelihood }
  \item{epsilon}{ Tolerance to determine if the optimization procedure converged }
}
\details{
The copula families which are implemented can be can be given specifying one of the strings "normal","frank","kimeldorf.sampson","gumbel","galambos","husler.reiss","bb1","tawn","bb2","bb3","bb5","bb6","bb7","normal.mix","bb4", "joe"}

\value{
Copula (from the chosen family) fitted by maximum likelihood.
}
\references{ 
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}

\keyword{distribution}
\keyword{multivariate}
