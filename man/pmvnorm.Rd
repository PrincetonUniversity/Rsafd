\name{pmvnorm}
\alias{pmvnorm}
\alias{qmvnorm}
\title{
Functions giving the values of the cdf and quantile and producing random samples of a multivariate distribution function
}
\description{
Compute the values of the cumulative distribution function of a multivariate normal distribution.
}
\usage{
pmvnorm(x, mean = rep(0, d), cov = diag(d), sd, rho, d = 2, sigma = cov,
        log = FALSE)
}
\arguments{
  \item{x}{
numeric vector of the values where the function is to be evaluated
}
  \item{mean}{
numeric vector giving the mean of the distribution
}
  \item{cov}{
matrix giving the variance - covariance of the distribution
}
  \item{sd}{
vector giving the standard deviations of the marginal distributions}
  \item{rho}{
Pearson correlation coefficient in the case of dimension 2
}
  \item{d}{
dimension of the distribution
}
  \item{sigma}{
same as \code{cov}
}
  \item{log}{
logical, not used in the computation
}
}
\author{
Rene Carmona, \email{rcarmona@princeton.edu}
}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
