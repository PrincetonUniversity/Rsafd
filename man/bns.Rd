\name{bns}
\alias{bns}
\alias{fns}
\alias{yns}

\title{Coupon Bond Price, Forward Rate, and Yield Curves in the Nelson - Siegel Family}
\description{
Compute the term structure of interest rate from the parametric  family proposed by Nelson and Siegel.
}
\usage{
fns(x, THETA)
yns(x, THETA)
bns(COUPON, AI, LIFE, X = 100, THETA = c(0.06, 0, 0, 1))
}
\arguments{
  \item{x}{Time to maturity (positive numeric)}
  \item{THETA}{Vector of the four parameters of a forward curve from the Nelson - Siegel family}  
  \item{COUPON}{Rate of the semi-annual coupon payment}
  \item{AI}{Accrued Interest}
  \item{LIFE}{Time to maturity of the bond}
  \item{X}{Nominal of the bond}
}

\value{
\item{fns}{Numerical value of the Nelson - Siegel forward curve determined by the parameter THETA for the time to maturity x}
\item{yns}{Numerical value of the Nelson - Siegel yield curve determined by the parameter THETA for the time to maturity x}
\item{bns}{Numerical value giving the price of a zero coupon bond with these parameters}
}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 3}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
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
\keyword{math}
