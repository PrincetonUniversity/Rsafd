\name{tailplot}
\alias{tailplot}
\alias{tailplot,gpd.upper.tail-method}
\alias{tailplot,gpd.lower.tail-method}
\alias{tailplot,gpd.two.tails-method}
\title{GENERALIZED EXTREME VALUE DISTRIBUTION}
\description{
        Cumulative distribution function and quantiles for the distribution
        of a GPD object (as produced for example by the semi-parametric estimation 
        procedure \code{gpd.tail}).
}
\usage{
tailplot(gpd.obj, optlog = NA, extend = 1.5, labels = TRUE, ...)
}
\arguments{ 
       \item{gpd.obj}{ An object of class \code{gpd} as created by \code{fit.gpd}.}
 \strong{Optional Arguments}
       \item{optlog}{ Character string which can be "x", "y", "xy", "yx", "" or "NA" indicating 
            which axis if any, should be labelled an a logarithmic scale.}
       \item{extend}{ Factor by which the x-axis should extend beyond the data} 
       \item{labels}{ Title and axis labels produced if TRUE} 
       \item{\dots}{ additional arguments, currently not used}
}
\value{
        Plot(s) of the estimate(s) of the tail(s) of the cumulative distribution function.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\section{Methods}{
\describe{
\item{\code{signature(gpd.obj = "gpd.upper.tail")}}{plot of the upper tail}
\item{\code{signature(gpd.obj = "gpd.lower.tail")}}{plot of the lower tail}
\item{\code{signature(gpd.obj = "gpd.two.tails")}}{plots of both tails}
}}
\seealso{\code{\link{fit.gpd}}, \linkS4class{gpd}}
\examples{
# One tail
X <- rpareto(1024)
x.est <- fit.gpd(X, tail = "upper")
tailplot(x.est)
}
\keyword{hplot}
\keyword{distribution}
\keyword{univar}
