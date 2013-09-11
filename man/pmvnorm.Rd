\name{pmvnorm}
\alias{pmvnorm}
\alias{qmvnorm}
\alias{rmvnorm}
\title{
Functions giving the values of the cdf and quantile and producing random samples of a multivariate distribution function
}
\usage{
pmvnorm(x, mean = rep(0, d), cov = diag(d), sd, rho, d = 2, sigma = cov, log = FALSE)
pmvnorm(n, mean = rep(0, d), cov = diag(d), sd, rho, d = 2, sigma = cov, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector of the values where the function is to be evaluated
}
  \item{n}{
integer giving the number of desired random samples
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
}
\references{
 R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\author{
Rene Carmona, \email{rcarmona@princeton.edu}
}
