\name{bscall}
\alias{bscall}
\alias{isig}
\title{Black-Scholes Price of an European Call Option}
\description{
Evaluate the Black-Scholes formula for the price of an European call option
}
\usage{
bscall(TAU = 0.03, K = 1, S, R = 0.1, SIG = 0.15)
isig(TAU, K, S, R, CALL, TOL = 1e-005, JMAX = 1000)
}

\arguments{
  \item{TAU}{Time to maturity (in years)}
  \item{K}{Strike of the option}
  \item{S}{Current value of the underlier of the option }
  \item{R}{Yearly short interest rate}
  \item{SIG}{Volatility of the underlier}
  \item{CALL}{Price of the call option}
}

\value{
  \code{bscall} gives the value of the call option as given by the Black-Scholes formula, while \code{isig} gives the implied volatility.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 5 & Appendix}
\examples{
BSTRG <- bscall(TAU=TRGSP[,3]/252,K=TRGSP[,2],S=TRGSP[,1],R=TRGSP[,4],SIG=TRGSP[,5])
range(TRGSP[,6] - BSTRG)
}
\keyword{math}