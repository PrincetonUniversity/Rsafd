\name{tail.index}
\alias{tail.index}
\alias{lambda}
\alias{tail.index,bb1.copula-method}
\alias{tail.index,bb2.copula-method}
\alias{tail.index,bb3.copula-method}
\alias{tail.index,bb4.copula-method}
\alias{tail.index,bb5.copula-method}
\alias{tail.index,bb6.copula-method}
\alias{tail.index,bb7.copula-method}
\alias{tail.index,empirical.copula-method}
\alias{tail.index,galambos.copula-method}
\alias{tail.index,gumbel.copula-method}
\alias{tail.index,husler.reiss.copula-method}
\alias{tail.index,tawn.copula-method}
\alias{lambda,bb1.copula-method}
\alias{lambda,bb2.copula-method}
\alias{lambda,bb3.copula-method}
\alias{lambda,bb6.copula-method}
\alias{lambda,bb7.copula-method}
\alias{lambda,empirical.copula-method}
\alias{lambda,frank.copula-method}
\alias{lambda,joe.copula-method}
\alias{lambda,kimeldorf.sampson.copula-method}
\title{Tail Indexes and Lambda Function of a Copula}
\description{
\code{tail.index} computes the lower and upper tail dependence indexes of a copula.
\code{lambda} computes, for an Archimedean copula with generator \eqn{\phi}{Phi},
the function \eqn{\lambda(t) = \phi(t)/\phi'(t)}{lambda(t) = Phi(t)/Phi'(t)}.
}
\usage{
\S4method{tail.index}{bb1.copula}(copula, ...)
\S4method{tail.index}{bb2.copula}(copula, ...)
\S4method{tail.index}{bb3.copula}(copula, ...)
\S4method{tail.index}{bb4.copula}(copula, ...)
\S4method{tail.index}{bb5.copula}(copula, ...)
\S4method{tail.index}{bb6.copula}(copula, ...)
\S4method{tail.index}{bb7.copula}(copula, ...)
\S4method{tail.index}{empirical.copula}(copula, ...)
\S4method{tail.index}{galambos.copula}(copula, ...)
\S4method{tail.index}{gumbel.copula}(copula, ...)
\S4method{tail.index}{husler.reiss.copula}(copula, ...)
\S4method{tail.index}{tawn.copula}(copula, ...)
lambda(copula, t)
}
\arguments{
  \item{copula}{object of class \code{copula} (or one of its subclasses) or \code{empirical.copula}}
  \item{t}{numeric vector of values in (0,1)}
  \item{\dots}{additional arguments, currently not used}
}
\value{
For the parametric copula families, \code{tail.index} returns a numeric vector with
components \code{lower.tail} and \code{upper.tail}. For an empirical copula, it
plots estimates of the lower and upper indexes as functions of the quantile level.

\code{lambda} returns a numeric vector of the same length as \code{t}.
}
\section{Methods}{
\code{tail.index} has methods for the classes \code{empirical.copula},
\code{galambos.copula}, \code{gumbel.copula}, \code{husler.reiss.copula},
\code{tawn.copula} and \code{bb1.copula}, ..., \code{bb7.copula}.
\code{lambda} has methods for the Archimedean families \code{frank.copula},
\code{kimeldorf.sampson.copula}, \code{joe.copula}, \code{bb1.copula},
\code{bb2.copula}, \code{bb3.copula}, \code{bb6.copula}, \code{bb7.copula},
and for \code{empirical.copula}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\linkS4class{copula}, \code{\link{Kendalls.tau}}}
\examples{
tail.index(gumbel.copula(1.5))
lambda(frank.copula(2), c(0.25, 0.5, 0.75))
}
\keyword{distribution}
\keyword{multivariate}
