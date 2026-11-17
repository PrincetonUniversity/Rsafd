\name{bivd-class}
\docType{class}
\alias{bivd-class}
\alias{gpd.bivd-class}
\alias{empirical.bivd-class}
\title{Classes of Bivariate Distributions}
\description{
The class \code{bivd} represents a bivariate distribution given by a copula and two
marginal distributions identified by their names and parameters. The class
\code{gpd.bivd} represents a bivariate distribution given by a copula and two
marginal distributions of class \code{gpd}. The class \code{empirical.bivd} holds
a bivariate sample.
}
\section{Objects from the Class}{
Objects of class \code{bivd} are created by the function \code{\link{bivd}}.
Objects of the other classes can be created by calls of the form \code{new("gpd.bivd", ...)}
and \code{new("empirical.bivd", ...)}.
}
\section{Slots}{
Class \code{bivd}:
\describe{
  \item{\code{copula}:}{object of class \code{copula}}
  \item{\code{Xmarg}, \code{Ymarg}:}{character strings giving the names of the marginal distributions}
  \item{\code{param.Xmarg}, \code{param.Ymarg}:}{numeric vectors giving the parameters of the marginal distributions}
}
Class \code{gpd.bivd}:
\describe{
  \item{\code{copula}:}{object of class \code{copula}}
  \item{\code{Xmarg}, \code{Ymarg}:}{objects of class \code{gpd} giving the marginal distributions}
}
Class \code{empirical.bivd}:
\describe{
  \item{\code{x}, \code{y}:}{numeric vectors giving the sample}
}
}
\section{Methods}{
See \code{\link{pbivd}}, \code{\link{dbivd}}, \code{\link{rbivd}},
\code{\link{contour.dbivd}} and \code{\link{persp.dbivd}}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 3}
\seealso{\linkS4class{copula}, \linkS4class{gpd}}
\examples{
showClass("bivd")
}
\keyword{classes}
