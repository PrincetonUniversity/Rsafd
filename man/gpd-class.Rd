\name{gpd-class}
\docType{class}
\alias{gpd-class}
\alias{gpd.upper.tail-class}
\alias{gpd.lower.tail-class}
\alias{gpd.two.tails-class}
\title{Classes of Distributions with Generalized Pareto Tails}
\description{
\code{gpd} is a virtual class (a class union) whose members \code{gpd.upper.tail},
\code{gpd.lower.tail} and \code{gpd.two.tails} hold the results of the
semi-parametric estimation performed by \code{\link{fit.gpd}} when the upper tail,
the lower tail, or both tails are fitted by generalized Pareto distributions.
}
\section{Objects from the Class}{
Objects are created by the function \code{\link{fit.gpd}}.
}
\section{Slots}{
\describe{
  \item{\code{tail}:}{character string, \code{"upper"}, \code{"lower"} or \code{"two"}}
  \item{\code{n}:}{length of the data}
  \item{\code{data}:}{numeric vector of the data}
  \item{\code{upper.exceedances}, \code{lower.exceedances}:}{data values beyond the thresholds}
  \item{\code{upper.threshold}, \code{lower.threshold}:}{thresholds}
  \item{\code{p.less.upper.threshold}, \code{p.larger.lower.threshold}:}{proportions of the data below the upper threshold and above the lower threshold}
  \item{\code{n.upper.exceedances}, \code{n.lower.exceedances}:}{numbers of exceedances}
  \item{\code{upper.method}, \code{lower.method}:}{estimation methods (\code{"ml"} or \code{"lmom"})}
  \item{\code{upper.par.ests}, \code{lower.par.ests}:}{parameter estimates of the generalized Pareto distributions}
  \item{\code{upper.converged}, \code{lower.converged}:}{logicals indicating whether the estimation procedures converged}
}
\code{gpd.upper.tail} only has the slots \code{tail}, \code{n}, \code{data} and the
\code{upper.*} slots, \code{gpd.lower.tail} only has the slots \code{tail}, \code{n},
\code{data} and the \code{lower.*} slots, and \code{gpd.two.tails} has all of them.
}
\section{Extends}{
\code{gpd.upper.tail}, \code{gpd.lower.tail} and \code{gpd.two.tails} are the members
of the class union \code{"gpd"}.
}
\section{Methods}{
See \code{\link{pgpd}}, \code{\link{dgpd}}, \code{\link{qgpd}}, \code{\link{rgpd}}
(and \code{\link{rgpd-methods}}) and \code{\link{tailplot}}.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 2}
\examples{
showClass("gpd.two.tails")
}
\keyword{classes}
