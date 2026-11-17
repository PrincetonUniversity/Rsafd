\name{Kend.table}
\alias{Kend.table}
\docType{data}
\title{Table of Copula Parameters Indexed by Kendall's Tau}
\description{
Data frame used by \code{\link{fit.copula}} to choose the initial values of the
parameters of the copula families from the empirical value of Kendall's tau.
}
\usage{Kend.table}
\format{
A data frame with 11 rows and 23 columns. The first column, \code{tau}, gives the
values 0, 0.1, ..., 1 of Kendall's tau. The other columns (\code{B1}, \code{B3},
\code{B4}, \code{B6}, \code{B7}, \code{B8}, \code{BB1.delta}, \code{bb1.theta},
\code{Twan.copula}, \code{TC.a.b}, \code{BB2.theta}, \code{BB2.delta},
\code{BB3.theta}, \code{BB3.delta}, \code{BB5.theta}, \code{BB5.delta},
\code{BB6.theta}, \code{BB6.delta}, \code{BB7.theta}, \code{BB7.delta},
\code{BB4.theta}, \code{BB4.delta}) give the corresponding values of the parameters
of the copula families.
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\seealso{\code{\link{fit.copula}}, \code{\link{Kendalls.tau}}}
\examples{
Kend.table[, 1:4]
}
\keyword{datasets}
