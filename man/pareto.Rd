\name{pareto}
\alias{pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}

\title{PARETO DISTRIBUTION}

\description{
     Density, cumulative distribution function, quantiles and random samples for the Pareto Distribution 
}
      
\usage{
         dpareto(x, m=0, lambda = 1, xi = 0)
         ppareto(q, m=0, lambda = 1, xi = 0)
         qpareto(p, m=0, lambda = 1, xi = 0)
         rpareto(n, m=0, lambda = 1, xi = 0)
}

\arguments{ 
         \item{x}{numeric vector. Missing values (NA's) are allowed}
         \item{q}{vector of quantiles. Missing values (NA's) are allowed}
         \item{p}{vector of probabilities. Missing values (NA's) are allowed}
         \item{n}{sample size. If \code{length(n)} is larger than 1, then \code{length(n)} random values are returned}

\item{\strong{Optional Arguments}}{}      
        \item{m}{Location parameter}
        \item{lambda}{Scale parameter}
        \item{xi}{Shape parameter (note the sign convention described in \code{SHAPE.XI})}
}

\details{  
    The length of vectors \code{m}, \code{lambda}, and \code{xi} should either be 1 or equal to the length of the vectors \code{x}, \code{q}, \code{p}, or \code{n}. Missing values and NA's are not allowed.
}

\value{    
    Density (\code{dpareto}), cumulative distribution function (\code{ppareto}),
    quantile (\code{qpareto}), or random sample (\code{rpareto})
     for the Generalized Extreme Value distribution with parameters
     \code{m}, \code{lambda}, and \code{xi}.

\bold{SIDE EFFECTS}
     The function \code{rpareto} causes the creation of \code{.Random.seed} if it does not already exist, otherwise its value is updated.
}

\examples{
x <- rpareto(200) # Generates a sample of size 200 from
               # Exponential distribution

# Evaluation at x = 1.0 of Pareto probability
# density functions with different parameters 
m <- c(0,0,0)
lambda <- c(1.3, 1.0, .7)
xi <- c(-0.2,0.0,0.2)
x <- rep(1.0,3)
ppareto(x,m,lambda,xi)
#[1] 0.5662428 0.6321206 0.7153720
}
\seealso{
    \code{\link{gev}}, \code{\link{gev.lmom}}, \code{\link{gev.ml}}
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 2
}
\keyword{distribution}
\keyword{datagen}
\keyword{univar}