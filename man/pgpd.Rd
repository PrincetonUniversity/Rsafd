\name{pgpd}
\alias{pgpd}
\alias{dgpd}
\alias{qgpd}
\alias{rgpd}
\alias{pgpd,gpd.upper.tail-method}
\alias{pgpd,gpd.lower.tail-method}
\alias{pgpd,gpd.two.tails-method}
\alias{dgpd,gpd.upper.tail-method}
\alias{dgpd,gpd.lower.tail-method}
\alias{dgpd,gpd.two.tails-method}
\alias{qgpd,gpd.upper.tail-method}
\alias{qgpd,gpd.lower.tail-method}
\alias{qgpd,gpd.two.tails-method}
\title{RANDOM GENERATOR, DENSITY, CUMULATIVE DISTRIBUTION AND QUANTILE FUNCTIONS OF A GPD OBJECT}
\description{
        Random number generator, density, cumulative distribution function and quantiles for the distribution
        of a GPD object (as produced for example by the semi-parametric estimation 
        procedure \code{fit.gpd}).
}
\usage{
rgpd(dist, n = 50)
pgpd(dist, x)
dgpd(dist, x)
qgpd(dist, x)
}
\arguments{
       \item{n}{Integer giving the size of the desired random sample}
       \item{x}{A numeric vector of values at which the density and the cdf are computed
         (\code{pgpd}, \code{dgpd}), or a numeric vector of probabilities at which
         the quantiles are computed (\code{qgpd})}
       \item{dist}{An object of class \code{gpd} like the output of \code{fit.gpd}.}
}
\value{
       Function \code{rgpd}  returns a vector of 
       length \code{n} comprising independent samples from the distribution 
       determined by \code{dist} 
       
       Functions \code{dgpd} and \code{pgpd} return a vector of the same 
       length as \code{x} comprising the values of the density and the cumulative distribution 
       function of the distribution determined by \code{dist} computed at the
        points \code{x}.

       Function \code{gqpd}  returns a vector of the same 
       length as \code{x} comprising the values of the quantiles of the distribution 
       determined by \code{dist} computed at the points \code{x}. It returns NA for p = 0
		and Inf for p = 1. So it is recommended to avoid these values.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\seealso{
 \code{\link{fit.gpd}}
}
\section{Methods}{
\code{pgpd}, \code{dgpd} and \code{qgpd} have methods for the classes
\code{gpd.upper.tail}, \code{gpd.lower.tail} and \code{gpd.two.tails}
(see \linkS4class{gpd}). The methods for \code{rgpd} are listed in
\code{\link{rgpd-methods}}.
}
\examples{
# One tail
data(BCofLRet)
data(CCofLRet)
NZ <- (BCofLRet != 0 & CCofLRet  != 0)
BLRet <- BCofLRet[NZ]
CLRet <- CCofLRet[NZ]
X <- BLRet[BLRet > 0]
x.est <- fit.gpd(X,tail="upper")

y <- c(10:500)/1000
plot(y, pgpd(x.est,y), log = "x", type = "l")
big.X <- sort(X) > 0.01
points((sort(X))[big.X], (ppoints(sort(X)))[big.X])

# random generation from this distribution:
X <- BLRet[BLRet != 0]
x.est <- fit.gpd(X, upper = 0.015, lower = -0.015, method = "lmom")
n <- length(X)
Y <- qgpd(x.est,runif(n))
# statistically equivalent to 
Y <- rgpd(x.est,n)
plot(X, ylim = c(-0.3,0.3))
plot(Y, col = 4, , ylim = c(-0.3,0.3))
}
\keyword{univar}
\keyword{distribution}
