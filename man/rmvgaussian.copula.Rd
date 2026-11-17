\name{rmvgaussian.copula}
\alias{rmvgaussian.copula}
\title{
Generate Random Samples from the One-Factor Multivariate Gaussian Copula
}
\description{
Generate Random Samples from the One-Factor Multivariate Gaussian Copula
}
\usage{
rmvgaussian.copula(N, d = 2, rho, Sigma)
}
\arguments{
  \item{N}{
\code{N} integer giving the number of scenarios generated (sample size)
}
  \item{d}{
\code{d} integer giving the dimension of each individual sample
}
  \item{rho}{
\code{rho} numeric in (0,1) giving the correlation parameter of the one-factor correlation matrix when \code{Sigma} is not provided
}
  \item{Sigma}{
\code{Sigma} d x d matrix giving the covariance or correlation matrix of the Gaussian copula
}
}
\details{
The matrix Sigma is used only when the parameter  \code{rho} is not provided. It is automatically normalized to a correlation matrix.
The function may fail if the matrix \code{Sigma} is neither a covariance matrix not a correlation matrix.
}
\value{
The function returns a N x d matrix of N samples of dimension d from the Gaussian copula with covariance or correlation matrix \code{Sigma} when the latter is provided, or with one-factor correlation matrix with \code{rho} off the diagonal.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3
}
\seealso{ \code{\link{rcopula}} }
\examples{
Nsim <- 1000
SD <- rmvgaussian.copula(Nsim, d=5, rho=0.7)
}
\keyword{distribution}
