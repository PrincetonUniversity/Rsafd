\name{fns}
\alias{fns}
\alias{yns}
\title{
Function computing the values of forward rate curve from the Nelson-Siegel parametric family
}
\description{
Compute a forward (interest rate) curve from the Nelson-Siegel family}
\usage{
fns(x, THETA)
}
\arguments{
  \item{x}{
time to maturity of the rate
}
  \item{THETA}{
four dimensional vector of the parameters of the curve
}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\author{
Rene A. Carmona
}
\examples{
THETA <- c(.07,-.03,.1,2.5)
XX <- seq(from=0,to=30,by=1/12)
FORW <- fns(XX,THETA)
YIELD <- yns(XX,THETA)
par(mfrow=c(1,2))
plot(XX,FORW,type="l",ylim=c(0,.25))
title("Example of a Nelson Siegel Forward Curve")
plot(XX,YIELD,type="l",ylim=c(0,.25))
title("Example of a Nelson Siegel Yield Curve")
par(mfrow=c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
