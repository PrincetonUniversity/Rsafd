\name{SHAPE.XI}
\alias{SHAPE.XI}
\docType{data}

\title{PARAMETRIZATION OF GENERALIZED ETREME VALUES AND PARETO DISTRIBUTIONS}

\description{For financial applications, it is customary to
   define the GEV and GPD distributions as we did in the text with the parameter xi. However, for other applications, e.g. hydrology, it is customary to use a parametrization based on the shape parameter k = -xi. The user  may choose to work with the parametrization he or she is more comfortable with
 by setting the global variable \code{SHAPE.XI} to \code{TRUE} if parameterization by xi is preferred, 
 and to \code{FALSE} if the alternative parameterization by k is desired.
}
      
\note{
    The values of most of the heavy tail distribution fitting functions of the library \code{Rsafd} are affected by the value of the boolean variable {\code{SHAPE.XI}}. Its default value is {\code{TRUE}}
}
\references{
 R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag}
\examples{
m <- c(0,0,0)
lambda <- c(1.3, 1.0, .7)
xi <- c(-0.2,0.0,0.2)
x <- rep(1.0,3)
pgev(x,m,lambda,xi)
SHAPE.XI <- FALSE
pgev(x,m,lambda,xi)
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{misc}
