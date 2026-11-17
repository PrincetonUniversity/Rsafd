\name{pcopula}
\alias{pcopula}
\alias{dcopula}
\alias{pcopula,bb1.copula-method}
\alias{pcopula,bb2.copula-method}
\alias{pcopula,bb3.copula-method}
\alias{pcopula,bb4.copula-method}
\alias{pcopula,bb5.copula-method}
\alias{pcopula,bb6.copula-method}
\alias{pcopula,bb7.copula-method}
\alias{pcopula,empirical.copula-method}
\alias{pcopula,frank.copula-method}
\alias{pcopula,galambos.copula-method}
\alias{pcopula,gumbel.copula-method}
\alias{pcopula,husler.reiss.copula-method}
\alias{pcopula,joe.copula-method}
\alias{pcopula,kimeldorf.sampson.copula-method}
\alias{pcopula,normal.copula-method}
\alias{pcopula,normal.mix.copula-method}
\alias{pcopula,tawn.copula-method}
\alias{dcopula,bb1.copula-method}
\alias{dcopula,bb2.copula-method}
\alias{dcopula,bb3.copula-method}
\alias{dcopula,bb4.copula-method}
\alias{dcopula,bb5.copula-method}
\alias{dcopula,bb6.copula-method}
\alias{dcopula,bb7.copula-method}
\alias{dcopula,frank.copula-method}
\alias{dcopula,galambos.copula-method}
\alias{dcopula,gumbel.copula-method}
\alias{dcopula,husler.reiss.copula-method}
\alias{dcopula,joe.copula-method}
\alias{dcopula,kimeldorf.sampson.copula-method}
\alias{dcopula,normal.copula-method}
\alias{dcopula,normal.mix.copula-method}
\alias{dcopula,tawn.copula-method}
\title{Cumulative Distribution Function and Density of a Copula}
\description{
Compute the values of a copula (\code{pcopula}) and of its density (\code{dcopula})
at the points \code{(u, v)} of the unit square.
}
\usage{
pcopula(copula, u, v)
dcopula(copula, u, v)
}
\arguments{
  \item{copula}{object of class \code{copula} (or one of its subclasses), or, for \code{pcopula}, of class \code{empirical.copula}}
  \item{u}{numeric vector of values in [0,1]}
  \item{v}{numeric vector of values in [0,1], of the same length as \code{u}}
}
\value{
Numeric vector of the values of the copula (\code{pcopula}) or of its density (\code{dcopula}).
}
\section{Methods}{
Both functions have methods for the parametric copula classes
\code{normal.copula}, \code{normal.mix.copula}, \code{frank.copula},
\code{kimeldorf.sampson.copula}, \code{gumbel.copula}, \code{galambos.copula},
\code{husler.reiss.copula}, \code{tawn.copula}, \code{joe.copula} and
\code{bb1.copula}, ..., \code{bb7.copula}. \code{pcopula} also has a method
for the class \code{empirical.copula}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\linkS4class{copula}, \code{\link{normal.copula}}, \code{\link{rcopula}}, \code{\link{contour.pcopula}}, \code{\link{persp.pcopula}}}
\examples{
COP <- gumbel.copula(1.4)
pcopula(COP, c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75))
dcopula(COP, c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75))
}
\keyword{distribution}
\keyword{multivariate}
