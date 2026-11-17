\name{AEP.wer.ts}
\alias{AEP.wer.ts}
\alias{MARKET.wer.ts}
\docType{data}
\title{Weekly Excess Returns for AEP and for a Market Portfolio}
\description{
Objects of class \code{timeSeries} giving weekly excess returns: \code{AEP.wer.ts} (units \code{AEP.wer}) and \code{MARKET.wer.ts} (units \code{Market.wer}).
}
\usage{
data(AEP.wer.ts)
data(MARKET.wer.ts)
}
\format{
Objects of class \code{timeSeries} with 831 weekly observations from 1995-01-09 to 2010-12-13.
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\examples{
plot(AEP.wer.ts)
plot(MARKET.wer.ts)
}
\keyword{datasets}
