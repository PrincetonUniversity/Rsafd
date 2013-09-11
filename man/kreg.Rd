\name{kreg}
\alias{kreg}
\title{ Multivariate Kernel Regression }
\description{
 Kernel regression for a given bandwidth and a given kernel function
}
\usage{
kreg(x, y, xpred = x, kernel = "gaussian", bandwidth)
}

\arguments{
  \item{x}{\code{n x p} matrix giving the values of the n observations of the \code{p} explanatory variables}
  \item{y}{\code{n} - dimensional vector of the observations of the response variable}
  \item{xpred}{\code{m x p} matrix giving the values of the m values of the \code{p} explanatory variables for which we want to estimate the response}
  \item{kernel}{Character string giving the name of the kernel function (see Details below)}
  \item{bandwidth}{Scalar giving the value of the bandwidth used in the kernel estimation}
}
\details{
The kernel functions implemented are \code{box}, \code{triangle},
and \code{parzen}. Their name should be passed to the parameter
\code{kernel}.
}
\value{
Returns a list with components
  \item{xpred}{matrix of the explanatory variables for which response variables are computed}
  \item{ypred}{Kernel estimates of the response variables for the values of the explanatory variables contained in \code{xpred}}
  \item{bandwidth}{The value of the bandwidth used in the kernel regression}

  Returns the fitted values for the values of the explanatory variables in the 
sample are returned if \code{xpred} is not given.
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
X <- matrix(data=rnorm(2048),ncol=2)
Y <- sin(X[,1]^2 + X[,2]^2) +rnorm(1024,sd=.2)
SEQU <- seq(from=0, to=3,length=30)
XPRED <- cbind(SEQU,SEQU)
REG <- kreg(X,Y,XPRED,kernel="gaussian", bandwidth=.3)
}
\keyword{regression}
