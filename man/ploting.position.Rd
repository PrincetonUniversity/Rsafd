\name{plotting.position}
\alias{plotting.position}
\alias{plotting.positions}
\alias{ploting.position.estimator.LMOM}
\title{POTTING POSITION ESTIMATOR}
\description{
        Computes plotting position estimate of the sample L-Moments.
}
\usage{
plotting.positions(sample, gamma = -0.35, delta = 0)
}
\arguments{ 
    \item{sample}{Numeric vector of data values. Missing values (NA's) are not allowed.}
     \item{gamma}{A scalar (see below).}
     \item{delta}{A scalar (see below).}
}
\value{    
    Vector of length four with the estimates l_1, l_2, tau_3,
    and tau_4, based of the plotting position.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\seealso{
  \code{gev.lmom}
}
\examples{
x <- rpareto(50,xi = 0.4)
PPLM <- plotting.positions(x)
PPLM
#    ell_1     ell_2     tau_3     tau_4 
#2.2074045 1.7393349 0.5459595 0.3498482 
SLM <- sample.LMOM(x)
SLM
#    ell_1     ell_2     tau_3     tau_4 
#2.6754741 1.7584511 0.5554916 0.3645081 
}
\keyword{distribution}
