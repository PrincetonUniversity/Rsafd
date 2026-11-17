\name{eda.shape}
\alias{eda.shape}
\title{Exploratory Data Analysis of the Shape of a Distribution}
\description{
Produce four plots of a numeric sample on the same page: a histogram, a box plot,
a kernel density estimate, and a normal Q-Q plot with a robust regression line.
}
\usage{
eda.shape(x, xname = deparse(substitute(x)))
}
\arguments{
  \item{x}{numeric vector}
  \item{xname}{character string used in the title of the plots}
}
\value{
Used for its side effect of producing the plots.
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{R. A. Carmona: Statistical Analysis of Financial Data in R, (2014) Springer Verlag, Chapter 1}
\seealso{\code{\link{shape.plot}}}
\examples{
eda.shape(DukeLRet)
}
\keyword{hplot}
\keyword{univar}
