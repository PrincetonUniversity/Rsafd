\name{DUKE.wer}
\alias{DUKE.wer}
\alias{DUKE1.wer}
\docType{data}
\alias{DUKE2.wer}
\title{Weekly Excess Returns for Duke Energy}
\description{
 Numeric vector giving the weekly excess returns for Duke Energy. The vectors DUKE1.wer and DUKE2.wer are subsets of DUKE.wer covering disjoint periods discussed in the reference text. 
}
\usage{
data(DUKE.wer)
data(DUKE1.wer)
data(DUKE2.wer)
}
\format{Numeric vectors of length 568 (\code{DUKE.wer}) and 236 (\code{DUKE1.wer}, \code{DUKE2.wer}).}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag
}
\examples{
plot(DUKE.wer,type="l")
}
\keyword{datasets}
