\name{sim.garch}
\alias{sim.garch}
\title{Simulation of a GARCH Model
}
\description{
Generates sample scenarios for a GARCH model}
\usage{
sim.garch(model, n.ahead = 1024, n.start = 0, n.sim = 1, presigmat = NA, preXt = NA, innovations = NA)
}
\arguments{
  \item{model}{A list comprising a number sig2 for the mean volatility, a vector ar containing the AR coefficients and a vector ma containing the MA coefficients}
  \item{n.ahead}{The length of the desired scenarios}
  \item{n.start}{The number of steps generation is performed and ignored for the purpose of reaching a stationary regime}
  \item{n.sim}{Number of desired scenarios}
  \item{presigmat}{Values of the conditional volatility to use as initialization}
  \item{preXt}{Values of the time series to use as initialization}
  \item{innovations}{random samples of the noise for the simulations}
}
\details{
n.start is set to 0 when presigmat and preXt are provided. The innovations are independent N(0,1) when not provided as argument}
\value{
Returns a list comprising
  \item{X.t }{A matrix with n.sim row (one row per Monte Carlo scenario generated) and n.ahead columns containing the values of the simulated time series}
  \item{sigma.t }{A matrix with n.sim row (one row per Monte Carlo scenario generated) and n.ahead columns containing the values of conditional volatility}
}
\seealso{
\code{\link{garch}}}

\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\examples{
LRET <- diff(log(EnronIndex.ts))
LRET.ar <- ar(LRET) 
LRES.g <- garch(LRET.ar$resid)
plot(LRES.g$sigma.t)

MODEL <- list(sig2=LRES.g$coef[1],ma=LRES.g$coef[2],ar=LRES.g$coef[3])
PRESIG <- seriesData(LRES.g$sigma.t)[LL-1]
PREX <- FITTED[LL-1]
SIM <- sim.garch(model=MODEL,n.ahead = 1024, n.start = 0, n.sim = 1, presigmat=PRESIG, preXt=PREX )

LRES.sim <- sim.garch(model=MODEL,n.ahead = 1024, n.start = 1000, n.sim = 1)
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{time series}
