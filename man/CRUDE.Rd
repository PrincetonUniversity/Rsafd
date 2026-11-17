\name{CRUDE}
\alias{CRUDE}
\title{Term Structure of Crude Oil Forward Prices}
\description{
\code{CRUDE} is a 335 x 12 numeric matrix, each row corresponds to a date, the first one being 4/18/1989, and the last one being 8/12/2002.  Each row of the matrix \code{CRUDE} gives the prices of the 12 futures contracts of crude oil as traded the day before. Two rows are used in the text as illustration of the nonlinear least squares function \code{nls} The whole matrix is used in a problem at the end of Chapter 5 on Local and NonParametric Regression.
}
\docType{data}
\usage{data(CRUDE)}
\format{A numeric matrix with 335 rows and 12 columns.}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 5
}
\keyword{datasets}
