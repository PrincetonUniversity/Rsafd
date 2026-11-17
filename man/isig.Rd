\name{isig}
\alias{isig}
\title{
Function computing the Implied Volatility of an Option
}
\usage{
isig(TAU, K, S, R, CALL, TOL = 1e-05, JMAX = 1000)
}
\arguments{
  \item{\code{TAU}}{Time to maturity (in years) of the option}
  \item{\code{K}}{Strike of the option}
  \item{\code{S}}{Value of the underlying at the time the option is priced}
  \item{\code{R}}{(Annualized) short interest rate}
  \item{\code{CALL}}{type of the option}
  \item{\code{TOL}}{Tolerance for the numeric approximations}
  \item{\code{JMAX}}{Maximum number of iterations}
}
\details{
Inverse the function Call numerically}
\value{An approximation of the value of the implied volatility}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
	R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag
}
\seealso{\code{bscall}}
\keyword{financial}
