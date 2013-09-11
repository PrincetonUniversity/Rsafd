\name{qqexp}
\alias{qqexp}
\alias{qqnorm}
\title{Empirical Q-Q PLOT WITH THE EXPONENTIAL and NORMAL DISTRIBUTIONS}
\description{
Produce a Q-Q plot by comparing the empirical quantiles of a data 
sample to the theoretical quantiles of the exponential or normal
distribution}

\usage{
qqexp(x, nq=50)
qqnorm(x)
}

\arguments{
  \item{x}{ Data numeric vector}
  \item{nq}{ Number of quantiles of the exponential distribution actually computed (equal to 50 by default}
}
\value{
  Used for its (graphical) side effects, one of them being the plot 
of a robust estimate of the line on which the points of the Q-Q plot should be
if the sample distribution was indeed from the exponential or normal family.
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{ R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\seealso{     
\code{\link{qqplot}}
 }

\examples{
TST <- rnorm(1024)
qqnorm(TST)

TST <- rgpd(1024)
qqnorm(TST)
}

\keyword{distribution}
\keyword{univar}
