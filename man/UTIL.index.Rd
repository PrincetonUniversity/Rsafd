\name{UTIL.index}
\alias{UTIL.index}
\alias{ENRON.index}
\alias{DUKE.index}
\alias{UTILITY.index}

\docType{data}
\title{Utilities Indexes}
\description{
Indexes computed daily in 1993 for DUKE, ENRON and the Utility sector. In the case of Enron and Duke, the index is merely the value of the stock multiplied by a factor appropriately chosen so that the resulting value is exactly 100 on a day chosen to initialize the index. These data are used in the textbook "R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag" to illustrate the problems with time series regression.}


\usage{
head(UTIL.index)
head(DUKE.index)
head(ENRON.index)
}
\format{
  \code{UTIL.index} is a data frame with 260 observations on the following 3 variables. 
The row names are the dates for which the indexes are computed. 
They start from 01/04/1993 and end 12/31/1993
  \describe{
    \item{\code{ENRON.index}}{a numeric vector}
    \item{\code{DUKE.index}}{a numeric vector}
    \item{\code{UTILITY.index}}{a numeric vector}
  }
}
\details{
\code{
           ENRON.index DUKE.index UTILITY.index
01/04/1993    135.0000   104.2857      99.94170
01/05/1993    135.3714   103.5714      99.49463
01/06/1993    132.8571   104.2857      99.86034
01/07/1993    130.7143   103.5714      98.70023
01/08/1993    126.8000   101.8000      97.93630
01/11/1993    127.5143   101.8000      98.69736
}
}
\source{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 4.
}
\examples{
plot(DUKE.index,UTIL.index)
}
\keyword{datasets}
