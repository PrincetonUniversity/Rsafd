\name{plot.gpd}
\alias{plot.gpd}

\title{PLOT OF THE TAIL OF A GENERALIZED PARETO DISTRIBUTION}

\description{
        Cumulative distribution function and quantiles for the distribution
        of a GPD object (as produced for example by the semi-parametric estimation procedure \code{gpd.tail}).
}

\usage{
plot.gpd(gpd.obj, tail="upper",  optlog = NA, extend = 1.5, labels = T, ...)
}

\arguments{
       \item{gpd.obj}{An object of class \code{gpd} as created by \code{gpd.tail}.}
       \item{tail}{Character string which can be either "upper" or "lower" depending upon
            the choice of the upper or lower part of the data}
       \item{optlog}{Character string which can be "x", "y", "xy", "yx" or "NA" indicating 
            which axis if any, should be labelled an a logarithmic scale.}
}

\value{
       Plot the upper or lower tail of the density inferred from the empirical data
       together with the tail of the theoretical GPD with the estimated shape parameter. 
}

\seealso{
    \code{\link{gpd.tail}}, \code{\link{shape.plot}} 
}
    
\examples{
# One tail
data(PCS)
PCS.index <- PCS[,2]
PCS.est <- gpd.tail(PCS.index)
}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{hplot}
