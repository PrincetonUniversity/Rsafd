\name{WSP}
\alias{WSP}
\alias{WSPLRet}
\docType{data}
\title{ Weekly quotes of the S&P 500 Index}
\description{ Weekly quotes of the S&P 500 index
from January 4, 1960 to June 2, 2003.}
\usage{data(WSP),data(WSPLR), data(WSP.lts),  }
\format{
WSP is the numeric vector of the raw values of the index, WSPLR is the vector of
(log) returns obtained by WSPLR <- diff(log(WSP)), WSPLRet is the short version (subset) 
of WSPLR used in Subsection 1.7.3 of the textbook as example, and WSP.lts is a list with a component
WSP.ts$pos containing the dates, and a component WSP.ts$data containing the entries of WSP.
This format is suitable for the function timeSeriesFromList which will create an object of class
timeSeries. See example below.
}
\details{
}
\source{
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
data(WSPLRet)
par(mfrow=c(1,2))
shape.plot(WSPLRet,tail="lower")
shape.plot(WSPLRet,tail="upper")
par(mfrow=c(1,1))
}
\keyword{datasets}
