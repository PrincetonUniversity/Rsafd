\name{gev.lmom}
\alias{gev.lmom}
\alias{gpd.lmom}

\title{L-MOMENT ESTIMATORS}

\description{
      Computes L-moment parameter estimates for GEV and GPD
      distributions.
}
      
\usage{
         gev.lmom(lmom)
         gpd.lmom(lmom, location = NA, sample = NA) 
}

\arguments{ 
       \item{lmom}{Vector of estimates of sample L-moments in this order: l_1, l_2, and tau_3. Missing values (NA's) are not allowed. These estimates can be obtained using the function \code{sample.LMOM} or the function \code{plotting.position}. For \code{gpd.lmom}, at least one of the two optional arguments \code{location} or \code{sample} should be specified.
}

For \code{gpd.lmom} function:
         \item{\code{location}}{If the value of the GPD location parameter $m$
         is known, it may be specified here.}
         \item{\code{sampl}}{Numeric vector of data values}.

\strong{OPTIONAL ARGUMENTS}       
        \item{\code{m}}{Location parameter}
        \item{\code{lambda}}{Scale parameter}
        \item{\code{xi}}{Shape parameter (note the sign convention
       described in \code{SHAPE.XI})}
}


\value{    
       A list containing the field \code{param.est} for the
       L-moment estimates of the location parameters \code{m},
       the scale parameter \code{lambda}, and the shape parameter \code{xi}.
}

\details{  
       The argument of the function should be a vector of L-moments. The function
       \code{sample.LMOM} computes unbiased estimates of sample
       L-moments, but other estimates can also be used for that
       purpose (for example, plotting position estimates). If the
       argument of the function is a vector with length greater
       than four, the function will assume that a sample is passed
       as a first argument instead of a vector of L-moment estimates, and a warning
       message will be given.
}

\seealso{
      \code{sample.LMOM}, \code{plotting.position}.
       \code{gev.ml}
}

\examples{
x <- rgev(500, lambda = 3.5, xi = 0.4)
lmomx <- sample.LMOM(x)
#lmomx
#    ell_1     ell_2     tau_3     tau_4 
#4.1185870 4.1703614 0.4757740 0.3459382 

gev.lmom(lmomx)
#$param.est
#         m     lambda         xi 
#-0.2008167  3.3246279  0.4271664 

x <- rpareto(500,xi = 0.4)
SLM <- sample.LMOM(x)

gpd.lmom(SLM,location=0)
#$param.est
#               m        lambda                 xi 
#0.0000000 1.0368496 0.4809592 

gpd.lmom(SLM,sample=x)
#$param.est
#           m       lambda           xi 
#-0.001486515  1.039880886  0.479828864 
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\keyword{univar}
