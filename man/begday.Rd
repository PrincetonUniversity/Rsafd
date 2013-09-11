\name{begday}
\alias{begday}
\alias{noon}

\title{ Manipulation of timeDate Objects}

\description{
  Given a timeDate, compute the timeDate corresponding to the beginning 
of the day, or noon, on the day identified by timeDate
}
\usage{
begday(DAY)
noon(DAY)
}
\arguments{
  \item{DAY}{ a timeDate object }
}
\value{  
  \item{begday}{A timeDate with the same year, month and day and hour:minute:second equal to 00:00:00}
  \item{noon}{A timeDate with the same year, month and day and hour:minute:second equal to 12:00:00}
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag}
\examples{
DAYCUT <- timeDate("02/15/1990")
noon(DAYCUT)
}
\keyword{ts}
\keyword{chron}
