\name{estimate.mix.gev}
\alias{estimate.mix.gev}

\title{MIX PARAMETER ESTIMATION FOR A GEV DISTRIBUTION}

\description{
      Computes MIX parameter estimates of a GEV
      distribution
}
      
\usage{
         estimate.mix.gev(sample, init.est = NA,  sampLmom = NA, epsilon = 1e-005)
}

\arguments{ 
       \item{sample}{Numeric vector of data values. Missing values (NA's) are not allowed.}
	\item{\strong{OPTIONAL ARGUMENTS}}{}
         \item{init.est}{If parameter estimates
         are available, they may be passed to the
         function. They will be used as initialization
         for the optimization routine in the MIX method. If no
         argument is specified, L-moment estimates will be
         computed and used to initialize the MIX routine.}
         \item{epsilon}{The relative tolerance of the optimization
         routine for MIX. The default value is 10^{-5}.}
}

\value{    
       A list whose first component is the vector \code{$param.est} containing the MIX estimates of the location parameter \code{m}, the scale parameter \code{lambda}, and the shape parameter \code{xi}, and whose second component is the logical \code{$converged} which is equal to \code{TRUE} or \code{FALSE} depending on whether or not the optimization routine converged.)
}

\details{  
      \code{estimate.mix.gev} is an implementation of the first of the
      two mixed methods introduced in Morrison's PhD and discussed in the text.
      The minimization routine \code{nlm} from the core \code{R} 
      distribution is used to find \code{nlm}
      the conditional maximum of the log-likelihood function.
      If \code{nlm} does not converge, a warning is issued.
}  

\seealso{
  \code{\link{gev.lmom}}, \code{\link{gev.ml}}
}
  
  \references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}

\examples{
x <- rgev(500, lambda = 3.5, xi = 0.4)

estimate.mix.gev(x)
#$param.est
#       m   lambda       xi 
#5.409156 6.855572 4.814889 
#
#$convereged
#[1] TRUE
}

\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{univar}
