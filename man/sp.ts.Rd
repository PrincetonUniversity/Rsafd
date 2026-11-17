\name{sp.ts}
\alias{sp_data}
\alias{sp.ts}
\docType{data}
\title{Daily Closing values of the SP500 index}
\description{
                \code{sp_data} is a numeric vector of the daily closing values of the SP500 index from July 2, 1962 to December 31, 2020. \code{sp.ts} is a \code{timeSeries} object constructed from the data in the numeric vector \code{sp_data}.
}
\usage{
	data(sp_data)
	data(sp.ts)
}
\format{
\code{sp_data} is a numeric vector of length 14727 and \code{sp.ts} is an object of class
\code{timeSeries} with the same 14727 values.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\examples{
plot(sp.ts)
}
\keyword{datasets}
