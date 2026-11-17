\name{seriesData}
\alias{seriesData}
\alias{seriesPositions}
\alias{seriesUnits}
\alias{is.timeSeries}
\alias{seriesData,timeSeries-method}
\alias{seriesPositions,timeSeries-method}
\alias{seriesUnits,timeSeries-method}
\alias{is.timeSeries,timeSeries-method}
\title{Accessor Functions for timeSeries Objects}
\description{
Extract the data, the positions (time stamps) or the units of a \code{timeSeries}
object, and test whether an object is a \code{timeSeries}.
}
\usage{
seriesData(object)
seriesPositions(object)
seriesUnits(object)
is.timeSeries(object)
}
\arguments{
  \item{object}{object of class \code{timeSeries}}
}
\value{
\code{seriesData} returns the numeric vector or matrix of the data,
\code{seriesPositions} returns the \code{timeDate} object of the positions,
\code{seriesUnits} returns the character vector of the units, and
\code{is.timeSeries} returns \code{TRUE}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{\code{\link{timeSeries}}, \linkS4class{timeSeries}}
\examples{
head(seriesData(CPN.ts))
seriesPositions(CPN.ts)[1:5]
seriesUnits(CPN.ts)
is.timeSeries(CPN.ts)
}
\keyword{ts}
