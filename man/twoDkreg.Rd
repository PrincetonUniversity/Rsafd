\name{twoDkreg}
\alias{twoDkreg}
\title{ Two Dimensional Gaussian Kernel Regression }
\description{
  Implementation (in R) of a two-dimensional kernel regression using the Gaussian kernel
intended for plotting purposes.
}
\usage{
twoDkreg(X, Y, B, XPRED = X, FITTED = F, N = 256, lims = apply(X, 2, range), X1NAME = "X1", X2NAME = "X2", YNAME = "Y", PLOT = T)
}

\arguments{
  \item{X}{ a n x 2 numeric matrix giving the n observations of the two explanatory variables }
  \item{Y}{ a numeric vector of length n giving the corresponding values of the response variable }
  \item{B}{ a positive number or a vector of two positive numbers used in the computation of the 
 bandwith(s). See details below.}
  \item{XPRED}{ a p x 2 numeric matrix giving the values where the regression is to be computed }
  \item{FITTED}{ Boolean, TRUE if the regression is to be computed at the values of the observa }
  \item{N}{ Regression surface is plotted on a regular N x N grid if PLOT is TRUE }
  \item{lims}{ Limits of the regular grid over which the regression surface is plotted }
  \item{X1NAME}{ Name of the first explanatory variable }
  \item{X2NAME}{ Name of the second explanatory variable }
  \item{YNAME}{ Name of the response variable }
  \item{PLOT}{ Boolean, TRUE if the regression surface is to be plotted }
}
\details{
  If B is missing, the one dimensional "asymptotically optimal" bandwidth for a Gaussian 
kernel is used for each explanatory variable. If a single number is provided for B, 
each explanatory variable has for bandwidth B times the univariate optimal bandwidth 
described earlier. Finally, if B is a vector of length two, its entries are used as bandwidths
for the two explanatory variables. So using a vector with two identical entries is the only
way to have the same bandwidth for the two explanatory variables.
A plot of the regression surface is produced as side effect if XPRED is not provided and 
if PLOT is TRUE. 
}
\value{
  Return a LIST comprising 
  \item{xpred }{ matrix of values of explanatory variables where the regression is computed}
  \item{ypred }{ regression estimate for the values of xpred}
  if XPRED is given
  \item{fitted }{ values of the regression estimate for the observations}
  \item{ssr }{ sum of squared errors}
  if FITTED is TRUE
  \item{gridx1 }{ values of the first explanatory variable in the grid}
  \item{gridx2 }{ values of the second explanatory variable in the grid}
  \item{y }{ values of the regression estimate over the grid}
  otherwise	
}
\seealso{ \code{\link{kreg}}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
X <- matrix(data=runif(2048,-4,+4),ncol=2)
Y <- sin(X[,1]^2 + X[,2]^2) +rnorm(1024,sd=.1)
twoDkreg(X,Y)
}
\keyword{regression}
