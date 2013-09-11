\name{fit.gpd}
\alias{fit.gpd}
\title{POT ESTIMATION OF A GENERALIZED PARETO DISTRIBUTION}

\description{
      Semi-parametric estimation based on the Peaks over Threshold (POT)
      method.of a distribution with upper and
      possibly lower tail have the same size as the tail 
      of a generalized Pareto distribution.
}
      
\usage{
         fit.gpd(data, tail="two", upper = NA, lower = NA, upper.method ="ml", lower.method="ml", plot = T, ...)
}

\arguments{ 
	\item{data}{A numerical vector of sample data values}
	\item{\strong{Optional Arguments}}{}
        	\item{tail}{Character string which can be "two" if two tails are estimated, "upper" if the upper tail is estimated, and "lower" if the lower tail is estimated} 
       	\item{upper}{The value of the upper threshold. A Generalized Pareto 
           distribution is fitted to the data exceeding this threshold.
           If its value is missing, \code{upper} is estimated so that either
           15\% of the data samples, or exactly 150 values, whichever number is greater,
           are greater than the upper threshold.}
       \item{lower}{The value of the lower threshold. A Generalized Pareto 
           distribution is fitted to the data below this threshold.
          If its value is missing, it is estimated so that either
          15\% of the data samples, or exactly 150 values, whichever number is smaller,
           is smaller than the lower threshold.}
       \item{upper.method}{Character string which can be either "ml" or "lmom" depending upon
            the choice of the method used to estimate the parameters of the generalized
            Pareto distribution fitted to the part of the data exceeding \code{upper}}
       \item{lower.method}{Character string which can be either "ml" or "lmom" depending upon
            the choice of the method used to estimate the parameters of the generalized
            Pareto distribution fitted to the part of the data below \code{lower}}
       \item{plot}{If \code{TRUE} (default), for each tail, a plot of the excesses over threshold
           against the quantiles of the GPD with the same shape parameter xi
           is produced. The fact that the points line up supports the assumption about the
           (generalized) Pareto nature of the tail in question.}
}

\value{    
       The value of \code{fit.gpd} is an object of class \code{gpd}
       containing slots for n -- length of the data, data --
       sorted vector of data points, upper.threshold -- the value of upper threshold, and
       upper.par.ests -- estimates of GPD shape and scale parameter for the tail.
       and lower.threshold and lower.par.ests representing
       threshold and parameter estimates for the lower tail.
}


\note{     If \code{plot = T} a plot of the exceedances over the threshold against the
     quantiles of the estimated GPD is produced for each of the tails studied.

The semiparametric estimation performed by this function was done by the functions
\code{pot.1tail.est} and \code{pot.2tails.est} in the original \code{EVANESCE} library.
These functions were included in the \code{FinMetrics} library of \code{S-Plus} under the
name of \code{fit.gpd} to avoid confusion with the similar function \code{gpd} contributed
by Mc Neil. We named it \code{fit.gpd} for backward compatibility
with \code{FinMetrics}.
}

\seealso{
\code{\link{gpd}},\code{\link{gpd.ml}}, \code{\link{gpd.lmom}}}

\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\examples{
SHAPE.XI <- TRUE
data(BCofLRet)
data(CCofLRet)
NZ <- (BCofLRet != 0 & CCofLRet  != 0)
BLRet <- BCofLRet[NZ]
CLRet <- CCofLRet[NZ]

# One tail
X <- BLRet[BLRet > 0]
x.est <- fit.gpd(X,tail="upper",upper=0.04)

# Two tails
X <- BLRet[BLRet != 0]
x.est <- fit.gpd(X, upper = 0.015, lower = -0.015, method = "lmom")
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{univar}
\keyword{distribution}
