\name{Kendalls.tau}
\alias{Kendalls.tau}
\alias{Spearmans.rho}
\alias{Kendalls.tau,bb1.copula-method}
\alias{Kendalls.tau,bb2.copula-method}
\alias{Kendalls.tau,bb3.copula-method}
\alias{Kendalls.tau,bb4.copula-method}
\alias{Kendalls.tau,bb5.copula-method}
\alias{Kendalls.tau,bb6.copula-method}
\alias{Kendalls.tau,bb7.copula-method}
\alias{Kendalls.tau,empirical.copula-method}
\alias{Kendalls.tau,frank.copula-method}
\alias{Kendalls.tau,galambos.copula-method}
\alias{Kendalls.tau,gumbel.copula-method}
\alias{Kendalls.tau,husler.reiss.copula-method}
\alias{Kendalls.tau,joe.copula-method}
\alias{Kendalls.tau,kimeldorf.sampson.copula-method}
\alias{Kendalls.tau,normal.copula-method}
\alias{Kendalls.tau,normal.mix.copula-method}
\alias{Kendalls.tau,tawn.copula-method}
\alias{Spearmans.rho,bb1.copula-method}
\alias{Spearmans.rho,bb2.copula-method}
\alias{Spearmans.rho,bb3.copula-method}
\alias{Spearmans.rho,bb4.copula-method}
\alias{Spearmans.rho,bb5.copula-method}
\alias{Spearmans.rho,bb6.copula-method}
\alias{Spearmans.rho,bb7.copula-method}
\alias{Spearmans.rho,empirical.copula-method}
\alias{Spearmans.rho,frank.copula-method}
\alias{Spearmans.rho,galambos.copula-method}
\alias{Spearmans.rho,gumbel.copula-method}
\alias{Spearmans.rho,husler.reiss.copula-method}
\alias{Spearmans.rho,joe.copula-method}
\alias{Spearmans.rho,kimeldorf.sampson.copula-method}
\alias{Spearmans.rho,normal.copula-method}
\alias{Spearmans.rho,normal.mix.copula-method}
\alias{Spearmans.rho,tawn.copula-method}
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
\section{Methods}{
Both functions have methods for all the parametric copula classes and for
the class \code{empirical.copula} (see \linkS4class{copula}).
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2014), Springer Verlag, Chapter 3
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
