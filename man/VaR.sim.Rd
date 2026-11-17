\name{VaR.sim}
\alias{VaR.sim}
\title{Monte Carlo Computation of the Value at Risk of a Portfolio of Two Assets}
\description{
Estimate by Monte Carlo simulation the value at risk and the expected shortfall of a
portfolio of two assets whose log-returns have marginal distributions estimated by
\code{\link{fit.gpd}} and whose dependence is given by a copula.
}
\usage{
VaR.sim(n, p, copula, x.est, y.est, lambda1, lambda2)
}
\arguments{
  \item{n}{integer giving the number of Monte Carlo samples}
  \item{p}{numeric vector of levels (probabilities) at which the value at risk and the
    expected shortfall are computed}
  \item{copula}{object of class \code{copula}}
  \item{x.est}{object of class \code{gpd} giving the distribution of the log-returns of the first asset}
  \item{y.est}{object of class \code{gpd} giving the distribution of the log-returns of the second asset}
  \item{lambda1}{weight of the first asset in the portfolio}
  \item{lambda2}{weight of the second asset in the portfolio}
}
\details{
\code{n} samples are generated from the copula, the log-returns are obtained with
\code{\link{qgpd}}, and the log-return of the portfolio is computed as
\code{log(lambda1*exp(X) + lambda2*exp(Y))}.
}
\value{
Named numeric vector containing the simulation size, the values at risk and the
expected shortfalls for the levels \code{p}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\code{\link{fit.gpd}}, \code{\link{fit.copula}}, \code{\link{rcopula}}}
\examples{
\donttest{
NZ <- (BCofLRet != 0 & CCofLRet != 0)
B.est <- fit.gpd(BCofLRet[NZ], upper = 0.04, lower = -0.045, plot = FALSE)
C.est <- fit.gpd(CCofLRet[NZ], upper = 0.04, lower = -0.045, plot = FALSE)
VaR.sim(n = 10000, p = c(0.01, 0.05), copula = gumbel.copula(1.4),
        x.est = B.est, y.est = C.est, lambda1 = 0.5, lambda2 = 0.5)
}
}
\keyword{distribution}
\keyword{multivariate}
