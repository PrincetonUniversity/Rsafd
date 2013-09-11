\name{Kendalls.tau}
\alias{Kendalls.tau}
\alias{Spearmans.rho}

\title{ Computation of the Kendall's TAU and the Spearman's RHO Correlation Coefficients }
\description{
Compute the Kendall's TAU and the Spearman's RHO of a copula. 
Give the empirical estimates when the copula is the empirical 
copula of a bivariate sample.  
}

\usage{
Kendalls.tau(copula, tol = 1e-05)
Spearmans.rho(copula, tol= 1e-05)
}
\arguments{
  \item{copula}{ Object of class copula }
  \item{tol}{ Small number giving the tolerance of the numerical approximation }
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 3
}
\examples{
SHAPE.XI <- TRUE
X <- rmvnorm(1024, mean=c(0,0),sd=c(1,1),rho=0.7)
X1 <- exp(X)[,1]
X2 <- exp(X)[,2]
X1.est <- fit.gpd(X1)
X2.est <- fit.gpd(X2)

U <- pgpd(X1.est,X1)
V <- pgpd(X2.est, X2)
EMPCOP <- empirical.copula(U,V)
ESTC <- fit.copula(EMPCOP,"gumbel")

Kendalls.tau(EMPCOP)
Kendalls.tau(ESTC)
Spearmans.rho(EMPCOP)
Spearmans.rho(ESTC)
}
\keyword{ distribution }
\keyword{ multivariate }
