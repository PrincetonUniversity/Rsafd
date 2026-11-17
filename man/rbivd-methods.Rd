\name{rbivd-methods}
\alias{rbivd-methods}
\alias{rbivd}
\alias{rbivd,bivd-method}
\alias{rbivd,gpd.bivd-method}
\title{Methods for Function \code{rbivd} in Package \pkg{Rsafd}}
\description{
 Methods for function \code{rbivd} in package \pkg{Rsafd}
}
\usage{
rbivd(dist, n)
}
\arguments{
  \item{dist}{object of class \code{bivd} or \code{gpd.bivd}}
  \item{n}{integer giving the size of the desired random sample}
}
\value{
A data frame with columns \code{x} and \code{y} containing the \code{n} samples.
}
\section{Methods}{
\describe{

\item{\code{signature(dist = "bivd")}}{
generate a random sample of size \code{n} from the bivariate distribution \code{dist}
}

\item{\code{signature(dist = "gpd.bivd")}}{
generate a random sample of size \code{n} from the bivariate distribution \code{dist}
}
}}
\seealso{\code{\link{bivd}}, \code{\link{pbivd}}}
\keyword{methods}
