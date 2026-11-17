\name{isig}
\alias{isig}
\title{
Function computing the Implied Volatility of an Option
}
\description{
Compute the implied volatility of a call option by numerically inverting the Black-Scholes formula.
}
\usage{
isig(TAU, K, S, R, CALL, TOL = 1e-05, JMAX = 1000)
}
\arguments{
  \item{TAU}{Time to maturity (in years) of the option}
  \item{K}{Strike of the option}
  \item{S}{Value of the underlying at the time the option is priced}
  \item{R}{(Annualized) short interest rate}
  \item{CALL}{type of the option}
  \item{TOL}{Tolerance for the numeric approximations}
  \item{JMAX}{Maximum number of iterations}
}
\details{
Inverse the function Call numerically}
\value{An approximation of the value of the implied volatility}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
	R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag
}
\seealso{\code{\link{bscall}}}
\keyword{financial}
