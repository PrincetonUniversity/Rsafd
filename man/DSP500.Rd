\name{DSP500}
\alias{DSP500}
\docType{data}
\title{
Vector of the daily close and log-returns of the S&P 500 index
}
\description{
Daily quotes of S&P 500 index from January 4, 1960 to June 6, 2003
}
\usage{data(DSP500)}
\format{A data frame with 10929 observations on the following 6 variables.}
\arguments{
	\item{\code{Date}}{A character string giving the date}
    \item{\code{Open}}{a numeric vector}
    \item{\code{High}}{a numeric vector}
    \item{\code{Low}}{a numeric vector}
    \item{\code{Close}}{a numeric vector}
    \item{\code{Volume}}{a numeric vector}
}
\details{
This is a large file. Use it at your own risk}
\source{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\examples{
head(DSP500) 
}
\keyword{datasets}
