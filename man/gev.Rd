\name{dgev}
\alias{gev}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}

\title{GENERALIZED EXTREME VALUE DISTRIBUTION}

\description{
Density, cumulative distribution function, quantiles and random sample
generation for the Generalized Extreme Value (GEV) distribution
}
\usage{
         dgev(x, m=0, lambda = 1, xi = 0)    
         pgev(q, m=0, lambda = 1, xi = 0)
         qgev(p, m=0, lambda = 1, xi = 0)
         rgev(n, m=0, lambda = 1, xi = 0)
}
\arguments{
  \item{x}{numeric vector. Missing values (NA's) are allowed}
  \item{q}{vector of quantiles. Missing values (NA's) are allowed}
  \item{p}{vector of probabilities. Missing values (NA's) are allowed}
  \item{n}{sample size. If \code{length(n)} is larger than 1, then \code{length(n)} random values are returned}
	\item{\strong{Optional Arguments}}{}      
        \item{m}{Location parameter. Can be a vector (see details)}
        \item{lambda}{Scale parameter. Can be a vector (see details)}
        \item{xi}{Shape parameter (note the sign convention described in \code{SHAPE.XI}). Can be a vector (see details)}
}
\details{
    The length of vectors \code{m}, \code{lambda}, and \code{xi} should either be 1 or equal to the length of the vectors \code{x}, \code{q}, \code{p}, or \code{n}. Missing values and NA's are not allowed.
}
\value{
   Density (\code{dgev}), cumulative distribution function (\code{pgev}),
   quantile (\code{qgev}), or random sample (\code{rgev})
   for the Generalized Pareto Distribution with parameters
   \code{m}, \code{lambda}, and \code{xi}.

\bold{SIDE EFFECTS} The function \code{rgpd} causes the creation of \code{.Random.seed} if it does not already exist, otherwise its value is updated.
}
\seealso{    \code{\link{rgev}}, \code{\link{gpd.lmom}}, \code{\link{gpd.ml}}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag,
}
\examples{
    x <- rgev(200) # Generates a sample of size 200 from
                               # the Gumbel distribution
    
# Evaluation at x = 1.0 of GEV density functions
# with different parameters
    m <- c(0,0,0)
    lambda <- c(1.3, 1.0, .7)
    xi <- c(-0.2,0.0,0.2)
    x <- rep(1.0,3)
    pgev(x,m,lambda,xi)
#    [1] 0.6480696 0.6922006 0.7522940
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{distribution}
\keyword{univar}