\name{normal.copula}
\alias{normal.mix.copula}
\alias{kimeldorf.sampson.copula}
\alias{gumbel.copula}
\alias{galambos.copula}
\alias{husler.reiss.copula}
\alias{tawn.copula}
\alias{frank.copula}
\alias{archm.copula}
\alias{bb1.copula}
\alias{bb2.copula}
\alias{bb3.copula}
\alias{bb4.copula}
\alias{bb5.copula}
\alias{bb6.copula}
\alias{joe.copula}

\title{ Copula Constructor }
\description{
  Create an object of class \code{copula} from a specific family
}
\usage{
normal.copula(delta),normal.mix.copula(p,delta1,delta2),
kimeldorf.sampson.copula(delta),gumbel.copula(delta),
galambos.copula(delta),husler.reiss.copula(delta),
tawn.copula(delta),frank.copula(delta),archm.copula(delta),joe.copula(delta),
bb1.copula(theta,delta),bb2.copula(theta,delta),bb3.copula(theta,delta),
bb4.copula(theta,delta),bb5.copula(theta,delta),bb6.copula(theta,delta)
}
\arguments{
  \item{delta}{ Positive number giving the parameter of the normal (<1), 
kimeldorf.sampson , gumbel (>=1),galambos, husler.reiss, tawn, 
frank, and bbx (>=1 for bb1 and bb6), copula families }
  \item{p}{ Number (between 0 and 1)giving the proportion in the normal mixture }
  \item{delta1}{ Positive number giving the parameter of the first normal entering the mixture }
  \item{delta2}{ Positive number giving the parameter of the second normal entering the mixture }
  \item{a}{ Positive number (<=1) giving the 1st parameter of the tawn copula family}
  \item{b}{ Positive number (<=1) giving the 2nd parameter of the tawn copula family}
  \item{r}{ Number (>=1) giving the 3rd parameter of the tawn copula family }
  \item{theta}{ 1st parameter of the bbx copula families (>=1 for bb3, bb5, bb6, bb7)}
}
\value{
  Return an object of class \code{copula} from the specified family
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in R, (2013), Springer Verlag, Chapter 3
}
\examples{
contour.copula(normal.copula(0.7))
persp.copula(gumbel.copula(1.7))
}
\keyword{distribution}

