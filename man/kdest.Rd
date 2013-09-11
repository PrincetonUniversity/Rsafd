\name{kdest}
\alias{kdest}
\title{ Gaussian Kernel 2-D Density Estimation }
\description{
Perform a two-dimensional kernel density estimation with a Gaussian
kernel.
}
\usage{
kdest(X, Y, H, N = 256, lims = c(range(X), range(Y)), COND = 0, PLOT = T, XNAME = "X", YNAME = "Y")
}
\arguments{
  \item{X}{ Numeric vector giving a sample from the first variable }
  \item{Y}{ Numeric vector giving a sample from the second variable }
  \item{H}{ Bandwidth }
  \item{N}{ The density estimate is computed at NxN regularly spaced points }
  \item{lims}{ Limits of the ranges over which the density is to be computed }
  \item{COND}{ Integer, 0 for joint density, 1 for conditional density of Y given X, 2 for conditional density of X given Y }
  \item{PLOT}{ Boolean, a surface plot of the density estimate is produced when TRUE }
  \item{XNAME}{ Character string to appear as label on the x-axis }
  \item{YNAME}{ Character string to appear as label on the y-axis }
}
\value{
List comprising
  \item{deltax }{Mesh size in the x-direction}
  \item{deltay }{Mesh size in the y-direction}
  \item{gridx }{Numerical vector giving the subdivision in the x-direction}
  \item{gridy }{Numerical vector giving the subdivision in the y-direction}
  \item{ z }{Matrix giving the values of the density estimate over the grid}
}
\seealso{ \code{\link{kreg}},\code{\link{twoDkreg}} }
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 3}
\examples{
X <- rnorm(1024,sd=1)
Y <- rnorm(1024,sd=4)
kdest(X,Y)
}

\keyword{ distribution }
