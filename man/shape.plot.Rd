\name{shape.plot}
\alias{shape.plot}

\title{ESTIMATE OF XI AS FUNCTION OF THE THRESHOLD}

\description{
        Plot of estimates of the shape parameter of one tail of the data.
}

\usage{
         shape.plot(data, tail="upper", method = "ml", from = 0.5, to = 0.98, nint = 30)
}

\arguments{
       \item{data}{A numeric vector of values at which the cdf is computed}
       \item{tail}{Character string which can be either "upper", "lower" or "two" depending upon
            the choice of the upper tail only,  the lower tail only, or both tails}
       \item{method}{Character string which can be either "ml" or "lmom" depending upon
            the choice of the method used to estimate the shape parameter}
       \item{from}{Percentage of the data values left as not belonging to the tail}
       \item{to}{Maximum percentage of the data values excluded from the tail}
       \item{nint}{Integer giving the number of thresholds and hence the number of estimates
            of the shape parameter}
}

\value{
        Plot of the estimates of the shape parameter xi as function of the threshold
        (negative of the threshold for the lower tail). The percentages of points
        exceeding the threshold are given on the upper part of each plot.
}

\note{
        Based on Mc Neil's original \code{shape} function.
}

\seealso{
    \code{\link{gpd.tail}} 
}

\examples{
PCS.index <- PCS[,2]
shape.plot(PCS.index)

shape.plot(WSPLRet,tail="two")
}
\references{ R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{hplot}
