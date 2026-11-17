\name{timeSeries}
\alias{timeSeries}
\alias{timeSeries,ANY-method}
\title{ Constructor for Objects of Class timeSeries }
\description{
Creates timeSeries objects}
\usage{
timeSeries(positions, data, units = NULL, zone = "GMT", FinCenter = "GMT", ...)
}
\arguments{
  	\item{data}{Numeric vector (univariate case) or matrix}
 	 \item{positions}{Vector of length dim(data)[1] of timeDate objects giving the time stamps}
  	\item{units}{character vector giving the names (if any) of the columns of the data matrix}
	\item{zone}{ "GMT" by default}
	\item{FinCenter}{"GMT" by default}
	\item{\dots}{additional arguments, currently not used}
}
\value{
Return an object of class timeSeries.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag}
\seealso{
    \code{\link{plot}},\code{\link[timeDate]{align}},  \code{\link{merge}}, \linkS4class{timeSeries}
}
\examples{
require(timeDate)
POS <- timeSequence(from="1959-01-01", to="1997-12-01", by="month")
co2.ts <- timeSeries(positions=POS,data=as.vector(co2), units="co2")
plot(co2.ts)
}
\keyword{ts}
