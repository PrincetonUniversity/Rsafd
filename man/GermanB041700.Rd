\name{GermanB041700}
\alias{GermanB041700}
\docType{data}
\title{German Treasury Bonds}
\description{
German Treasury bonds outstanding on April 17, 2000}
\usage{data(GermanB041700)}
\format{
  A data frame with 47 observations on the following 8 variables.
  \describe{
    \item{\code{Issue}}{Label of the bond (often the year of issuance}
    \item{\code{Coupon}}{Rate of (yearly) coupon payments}
    \item{\code{Maturity}}{Year the bond matures}
    \item{\code{Price}}{Price of the bond (per nominal of 100}
    \item{\code{Intrst.Yield}}{Yield of the bond on that day}
    \item{\code{Redemp.Yield}}{Redemption yield}
    \item{\code{Accrud.Intrst}}{Accrued interests accrued on 4/17/2000 from the date of the last coupon payment}
    \item{\code{Life}}{Time (in years) remaining on 4/17/2000 until the maturity of the bond}
  }
}
\details{
   Issue Coupon Maturity  Price Intrst.Yield Redemp.Yield Accrud.Intrst  Life
1   1992   8.00     2002 105.28         7.60        5.202          7.67  2.04
3   1993   6.75     2003 104.37         6.47        5.158          6.47  3.04
5   1999   3.75     2009  90.43         4.15        5.135          1.07  8.72
7     G3   3.00     2010  77.00         3.90        6.080          0.12 10.46
9     G4   3.00     2010  77.00         3.90        6.080          0.12 10.46
11    G5   3.00     2010  77.00         3.90        6.080          0.12 10.46

}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag,, Chapter 4
}
\examples{
data(GermanB041700)
}
\keyword{datasets}
