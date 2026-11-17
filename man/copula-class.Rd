\name{copula-class}
\docType{class}
\alias{copula-class}
\alias{empirical.copula-class}
\alias{normal.copula-class}
\alias{normal.mix.copula-class}
\alias{frank.copula-class}
\alias{kimeldorf.sampson.copula-class}
\alias{gumbel.copula-class}
\alias{galambos.copula-class}
\alias{husler.reiss.copula-class}
\alias{tawn.copula-class}
\alias{joe.copula-class}
\alias{bb1.copula-class}
\alias{bb2.copula-class}
\alias{bb3.copula-class}
\alias{bb4.copula-class}
\alias{bb5.copula-class}
\alias{bb6.copula-class}
\alias{bb7.copula-class}
\alias{show,copula-method}
\title{Classes of Copulas}
\description{
The class \code{copula} holds the parameters of a parametric copula. Each parametric
family is represented by a subclass of \code{copula}: \code{normal.copula},
\code{normal.mix.copula}, \code{frank.copula}, \code{kimeldorf.sampson.copula},
\code{gumbel.copula}, \code{galambos.copula}, \code{husler.reiss.copula},
\code{tawn.copula}, \code{joe.copula}, and \code{bb1.copula}, ..., \code{bb7.copula}.
The class \code{empirical.copula} holds a bivariate sample used to define an empirical copula.
}
\section{Objects from the Class}{
Objects of the parametric classes are created by the constructors documented in
\code{\link{normal.copula}} or returned by \code{\link{fit.copula}}. Objects of class
\code{empirical.copula} are created by \code{\link{empirical.copula}}.
}
\section{Slots}{
Class \code{copula} and its subclasses:
\describe{
  \item{\code{parameters}:}{numeric vector of the parameters of the copula}
  \item{\code{param.names}:}{character vector of the names of the parameters}
  \item{\code{param.lowbnd}:}{numeric vector of lower bounds for the parameters}
  \item{\code{param.upbnd}:}{numeric vector of upper bounds for the parameters}
  \item{\code{message}:}{character string describing the copula family}
}
Class \code{empirical.copula}:
\describe{
  \item{\code{x}:}{numeric vector of the first coordinates of the sample}
  \item{\code{y}:}{numeric vector of the second coordinates of the sample}
}
}
\section{Extends}{
Each parametric family class extends \code{"copula"} directly.
}
\section{Methods}{
\describe{
  \item{show}{\code{signature(object = "copula")}: prints the description of the copula family (slot \code{message}) and the names and values of the parameters}
}
See \code{\link{pcopula}}, \code{\link{dcopula}}, \code{\link{rcopula}},
\code{\link{Kendalls.tau}}, \code{\link{Spearmans.rho}}, \code{\link{tail.index}},
\code{\link{lambda}}, \code{\link{contour.pcopula}} and \code{\link{persp.pcopula}}
for the other methods.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\examples{
showClass("copula")
COP <- gumbel.copula(1.5)
COP
}
\keyword{classes}
