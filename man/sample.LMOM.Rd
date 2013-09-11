\name{sample.LMOM}
\alias{sample.LMOM}

\title{SAMPLE L-MOMENTS}

\description{
      Computes unbiased estimates of mean, second L-moment,
        L-skewness, and L-kurtosis.
}
      
\usage{
         sample.LMOM(x)
}

\arguments{ 
         \item{x}{numeric vector of data values. Missing values (NA's) are not allowed}
}

\value{    
    Vector of length four with the estimates for l_1, l_2, tau_3,
    and tau_4.
}

\seealso{
  \code{\link{gev.lmom}}
}

\examples{
x <- rgpd(1000)
sample.LMOM(x) #Compute L-moments of the sample
Mean (l_1) L-mom 2 (l_2)    L-skewness    L-kurtosis 
 1.0144528     0.5045187     0.3278689     0.1626310 

# Theoretical L-Moments of GPD with parameters
# m=0, lambda = 1, and k = 0 are:
# Mean (l_1) L-mom 2 (l_2) L-skewness L-kurtosis
#      1           1/2         1/3      1/6
}
\references{ R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{distribution}
\keyword{univar}
