\name{dump.timeSeries}
\alias{dump.timeSeries}
\alias{timeSeriesFromList}
\title{Conversion between timeSeries Objects and Lists}
\description{
\code{dump.timeSeries} writes a \code{timeSeries} object to a file, in the format of
\code{\link{dump}}, as a list with components \code{pos} (the positions as character
strings) and \code{data}. \code{timeSeriesFromList} creates a \code{timeSeries}
object from such a list.
}
\usage{
dump.timeSeries(TS, FILENAME)
timeSeriesFromList(LIST)
}
\arguments{
  \item{TS}{object of class \code{timeSeries}}
  \item{FILENAME}{character string giving the name of the file}
  \item{LIST}{list with a component \code{pos} giving the positions and a component
    \code{data} giving the data}
}
\value{
\code{dump.timeSeries} returns invisibly the list with components \code{pos} and
\code{data}. \code{timeSeriesFromList} returns an object of class \code{timeSeries}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\seealso{\code{\link{timeSeries}}}
\examples{
DM.ts <- DesMoines.ts[1:100, ]
FILE <- tempfile(fileext = ".R")
DM.list <- dump.timeSeries(DM.ts, FILE)
unlink(FILE)
DM2.ts <- timeSeriesFromList(DM.list)
}
\keyword{ts}
