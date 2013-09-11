\name{lm.diag}
\alias{lm.diag}
\title{ Diagnostics of a linear model }
\description{
  Wrapper on the diagnostic functions of lm objects.
}
\usage{
lm.diag(LM)
}
\arguments{
  \item{LM}{ should  be an object of class lm, as the output of the function lm. }
}
\details{
  This function is just a wrapper for the diagnostic functions of lm objects.
}
\value{
  Return a list comprising
  \item{hat}{vector containing the diagonal of the hat matrix (see the hat function).}
  \item{std.res}{vector or matrix containing the standardized residuals. This uses hat and std.dev to standardize the residuals.}
  \item{stud.res}{vector or matrix containing the studentized residuals. This uses hat and a different estimate of the standard deviation.}
  \item{cooks}{vector or matrix containing Cook's distance for each observation. The element in the i-th row and j-th column is the measure of the distance between the parameter estimates for the j-th regression with and without the i-th observation}
  \item{covratio}{covariance ratio for the parameter estimates} 
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{ Rene Carmona }
\seealso{ See Also as \code{\link{ls.diag}} }
\examples{
X <- seq(from=1,to=10,length=1024)+rnorm(1024,mean=0,sd=0.2)
Y <- 2*X+1 +rnorm(1024,mean=0,sd=.2)
TST.lm <- lm(Y~X)
plot(lm.diag(TST.lm)$hat)
}
\keyword{ regression }
