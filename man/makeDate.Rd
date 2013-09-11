\name{makeDate}
\alias{makeDate}
\title{
Function returning an object of class timeDate from a character string
}
\usage{
makeDate(x, in.format = "d-m-y")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
character string giving a date
}
  \item{in.format}{
standard date format
}
}
\details{
Assume that x a vector of character strings representing dates. The argument in.format gives the format of the entries of x
 in.format = c("Ymd","Y-m-d","y-m-d","Y/m/d","y/m/d","d-m-Y","d-m-y","d/m/Y","d/m/y")
 where Y is a string with four digits, y and d are one or two digit strings 
 representing integers, and m is either a one or two digit string 
 (like 7 or 12) or a three character string like Apr
 Returns a vector of character strings of lengths 10, xxxx-xx-xx for Y-m-d 
 like 1997-01-07  
}
\author{Rene Carmona, \email{rcarmona@princeton.edu}}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag, Chapter 6}
\examples{
}
