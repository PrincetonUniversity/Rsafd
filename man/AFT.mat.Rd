\name{AFT.mat}
\alias{AFT.mat}
\alias{MORN.mat}
\alias{SAFT.mat}
\alias{SMORN.mat}
\docType{data}
\title{Matrices of Afternoon High Frequency S&P Indicators}
\description{
AFT.mat is a 59 x 6 matrix of the afternoon values of the 6 indicators computed described in the reference text, for 59 trading days 
The matrices MORN.mat, SAFT.mat and SMORN.mat have similar structure}
\details{
The columns are named after the indicators \code{range}, \code{nbticks}, \code{vol},   \code{volratio}, \code{l2slope} and  \code{ticksep} described in the reference text.
}
\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\examples{
head(AFT.mat)
}
\keyword{datasets}
