\name{gev.ml}
\alias{gev.ml}
\alias{gpd.ml}

\title{Maximum Likelihood Parameter Estimates}

\description{
      Computes maximum likelihood estimates of the parameters of GPD
      and GEV distributions.
}
      
\usage{
         gev.ml(sample, init.est = <<see below>>, epsilon = 1e-005)
         gpd.lm(sample, location = NA, init.est = <<see below>>, epsilon=1e-005)
}

\arguments{ 
         \item{sample}{Numerical vector of data values. Missing values (NA's) are not allowed}
         \item{init.est}{A vector of parameter estimates
         obtained by a different method may be passed to the
         function. It will be used as starting point
         for the optimization routine attempting to maximize the 
         likelihood. If no \code{init.est}
         argument is specified, L-moment estimates will be
         computed and used for initialization purposes.}
         \item{epsilon}{The relative tolerance of the optimization
         routine in the maximization of the likelihood. 
         The default value is 10^{-5}.}
         \item{location}{(for \code{gpd.lm} only) \code{epsilon}}
         If the location parameter \code{m} of a GPD
         is known, it may be specified here. The remaining two parameters
         will then be estimated by maximum likelihood.}


\value{    
       A list with fields \code{param.est} containing the \code{param.est}
       maximum likelihood estimates of the location parameter \code{m}, the scale 
       parameter \code{lambda}, and the shape parameter \code{xi}, and 
       \code{converged} which is equal to \code{TRUE}
       or \code{FALSE} depending on whether or not the optimization routine
       converged.
}

\details{  
      These functions compute the log-likelihood from the sample data, and 
      then find the parameters which maximize
      this log-likelihood. The nonlinear minimization routine \code{nlm} 
      from the core \code{R} distribution is used to find
      the maximum of the log-likelihood function. If \code{nlm} does not converge,
      a warning is issued.
}  

\seealso{
  \code{\link{gev.lmom}}, \code{\link{estimate.mix.gev}}
}
 
\examples{
x <- rgev(500, lambda = 3.5, xi = 0.4)
gev.ml(x)
#$param.est
#         m     lambda         xi 
#-0.1632381  3.2012076  0.3976820 

x <- rpareto(500, lambda = 3.5, xi = 0.4)
gpd.ml(x)$param.est
#$param.est
#          m      lambda          xi 
#0.003792366 3.079291618 0.433570363 
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 2}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{univar}
\keyword{distribution}