\name{timeSeries-class}
\docType{class}
\alias{timeSeries-class}
\alias{[,timeSeries-method}
\alias{[,timeSeries,ANY,ANY,ANY-method}
\alias{Arith,timeSeries,numeric-method}
\alias{Arith,numeric,timeSeries-method}
\alias{Arith,timeSeries,timeSeries-method}
\alias{Ops,timeSeries,timeSeries-method}
\alias{Math,timeSeries-method}
\alias{head,timeSeries-method}
\alias{tail,timeSeries-method}
\title{Class timeSeries}
\description{
Objects of class \code{timeSeries} hold the values of one or several time series
together with the time stamps (positions) at which they are observed.
}
\section{Objects from the Class}{
Objects are created by the function \code{\link{timeSeries}}.
}
\section{Slots}{
\describe{
  \item{\code{positions}:}{object of class \code{timeDate} giving the time stamps}
  \item{\code{data}:}{numeric vector or matrix giving the values of the series (one column per variable)}
  \item{\code{units}:}{character vector giving the names of the variables}
  \item{\code{zone}:}{character string giving the time zone}
  \item{\code{FinCenter}:}{character string giving the financial center}
}
}
\section{Methods}{
\describe{
  \item{[}{\code{signature(x = "timeSeries")}: \code{x[i, j]} returns the \code{timeSeries}
    made of the rows \code{i} and the columns \code{j} of \code{x}}
  \item{Arith}{\code{signature(e1 = "timeSeries", e2 = "numeric")},
    \code{signature(e1 = "numeric", e2 = "timeSeries")},
    \code{signature(e1 = "timeSeries", e2 = "timeSeries")}: arithmetic operations
    applied to the data; the positions of the first \code{timeSeries} argument are kept
    and the units are set to \code{"TBA"}}
  \item{Ops}{\code{signature(e1 = "timeSeries", e2 = "timeSeries")}: same for the other operators}
  \item{Math}{\code{signature(x = "timeSeries")}: mathematical functions (\code{log},
    \code{exp}, \code{abs}, ...) applied to the data; the units are set to \code{"TBA"}}
  \item{head}{\code{signature(x = "timeSeries")}: first rows of the data}
  \item{tail}{\code{signature(x = "timeSeries")}: last rows of the data}
}
See also \code{\link{seriesData}}, \code{\link{plot}}, \code{\link{diff}},
\code{\link{merge}}, \code{\link{concat}} and \code{\link{sstl}}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag}
\examples{
showClass("timeSeries")
head(CPN.ts)
LRET <- diff(log(CPN.ts))
tail(LRET)
}
\keyword{classes}
\keyword{ts}
