\name{Rsafd-package}
\alias{Rsafd-package}
\alias{Rsafd}
\docType{package}
\title{
Method Functions and Data Sets for the Book Statistical Analysis of Financial Data
}
\description{
ALl the illustrations, examples and problems of the book can be run with this package}
\details{
\tabular{ll}{
Package: \tab Rsafd\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-09-15\cr
License: \tab GPL version 2 or newer\cr
}
Provides R versions of all the functions and the data sets used in the book\cr 
	'Statistical Analysis of Financial Data in S-Plus. 
}
\author{
Rene Carmona, \email{rcarmona@princeton.edu}
}
\references{
R.A. Carmona, Statistical Analysis of Financial Data in S-Plus, 2004, Springer Verlag.
}
\examples{
data(BCofLRet)
data(CCofLRet)
NZ <- (BCofLRet != 0 & CCofLRet  != 0)
BLRet <- BCofLRet[NZ]
CLRet <- CCofLRet[NZ]

par(mfrow=c(1,2))
shape.plot(BLRet,tail="lower")
shape.plot(BLRet,tail="upper")
par(mfrow=c(1,1))

B.est <- fit.gpd(BLRet, upper = 0.04, lower = -0.045)

tailplot(B.est)

C.est <- fit.gpd(CLRet, upper = 0.04, lower = -0.045)

U <- pgpd(B.est, BLRet)
V <- pgpd(C.est,CLRet)
plot(U,V)
EMPCOP <- empirical.copula(U,V)

FAM <- "gumbel"
ESTC <- fit.copula(EMPCOP,FAM)

print(ESTC)
Kendalls.tau(EMPCOP)
Kendalls.tau(ESTC, tol = 1e-2)
Spearmans.rho(EMPCOP)
Spearmans.rho(ESTC)
}
\keyword{univar}
\keyword{multivariate}
\keyword{ts}
\keyword{distribution}
\keyword{ package }
