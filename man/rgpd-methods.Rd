\name{rgpd-methods}
\alias{rgpd-methods}
\alias{rgpd,gpd-method}
\alias{rgpd,gpd.lower.tail-method}
\alias{rgpd,gpd.two.tails-method}
\alias{rgpd,gpd.upper.tail-method}
\title{Methods for Function \code{rgpd} in Package \pkg{Rsafd}}
\description{
 Methods for function \code{rgpd} in package \pkg{Rsafd}
}
\section{Methods}{
\describe{

\item{\code{signature(dist = "gpd")}}{
generate a random sample of size \code{n} from the distribution fitted by \code{\link{fit.gpd}}
}

\item{\code{signature(dist = "gpd.lower.tail")}}{
generate a random sample of size \code{n} from the distribution fitted by \code{\link{fit.gpd}}
}

\item{\code{signature(dist = "gpd.two.tails")}}{
generate a random sample of size \code{n} from the distribution fitted by \code{\link{fit.gpd}}
}

\item{\code{signature(dist = "gpd.upper.tail")}}{
generate a random sample of size \code{n} from the distribution fitted by \code{\link{fit.gpd}}
}
}}
\keyword{methods}
