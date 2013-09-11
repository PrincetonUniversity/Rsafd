\name{wCAPM.ts}
\alias{wCAPM.ts}
\alias{wCAPM1.ts}
\alias{wCAPM2.ts}
\docType{data}
\title{Weekly Excess Returns for Tests of CAPM}
\description{
\code{timeSeries} object with a slot \code{pos} and a slot \code{data} 
giving the weekly excess returns of the energy companies AEP, DUKE, PCG, SO, TXU and TJU, 
together with the weekly returns of the S&P500 chosen as a proxy for the market 
portfolio and the yield on
a 3 month T-Bill given for discounting purposes.
}
\arguments{
\item{\code{pos}}{Last day of the week on which the returns are computed}
\item{\code{data}}{Matrix with one column for each of the variables described above}
For the purpose of the static analysis done in the textbook R.A. Carmona, Statistical 
Analysis of Financial Data in R, and for the dynamic analysis done later in the book based on Kalman filtering theory, there are 3 timeSeries with 8 variables each 
\item{\code{wCAPM}}{Data for the period starting on 22/10/1993 and ending 1/7/2003}
\item{\code{wCAPM1}}{Data for the period starting on 1/1/1995 and ending 1/1/1999}
\item{\code{wCAPM2}}{Data for the period starting on 1/1/1999 and ending 1/1/2003}
}

\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag
}
\examples{
plot(wCAPM.ts[,1])
}
\keyword{datasets}
