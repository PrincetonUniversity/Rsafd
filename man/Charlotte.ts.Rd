\name{Charlotte.ts}
\alias{TEMPS.ts}
\alias{TEMPS1.ts}
\alias{TEMPS2.ts}
\alias{TEMPS3.ts}
\alias{Covington.ts}
\alias{DesMoines.ts}
\alias{Portland.ts}
\alias{Reno.ts}
\alias{LasVegas.ts}
\alias{Sacramento.ts}
\alias{Newark.ts}
\alias{LaGuardia.ts}
\alias{Philadelphia.ts}
\docType{data}
\title{Daily Temperature Data}
\description{{\code{timeSeries}} objects giving the daily temperatures at various meteorological stations.}
\format{
Objects with a slot positions and a slot data. \code{positions} gives the dates in timeDate format, 
  (e.g. from 1960-01-01 to 199-03-31 in the case of Charlotte)
  and data is a numeric vector or matrix, e.g. matrix with 13969 rows and 3 
  columns in the case of Charlotte
}
\arguments{
	\item{Charlotte.ts}{The data slot has 3 variables  "MinT" "MaxT" "AvgT". On any given day, they give the minimum temperature, the maximum temperature, and the average temperature on that day. Recall that "AvgT" is merely computed as ("MinT"+"MaxT")/2}
	\item{TEMPS.ts}{The 3 variables forming the data slot are "CVG" for Covington, "DM" for Des Moines and "PDX" for Portland}
	\item{TEMPS1.ts}{The 3 variables forming the data slot are "RENO" for Reno, "LV" for Las Vegas and "SACRA" for Sacramento}
	\item{TEMPS2.ts}{The 3 variables forming the data slot are "NEW" for Newark, "LG" for La Guardia and "PHIL" for Philadelphia}
	\item{TEMPS3.ts}{The 2 variables forming the data slot are "AIX" for Aix en Provence, "STETIENNE" for Saint Etienne are daily temperature reading simulated for the purpose of a problem set}
	\item{Covington.ts}{The data slot is a vector giving the daily average temperature at the Covington meteorological station. First component of the multivariate timeSeries TEMPS.ts}
	\item{DesMoines.ts}{The data slot is a vector giving the daily average temperature at the Des Moines meteorological station. Second component of the multivariate timeSeries TEMPS.ts}
	\item{PDX.ts}{The data slot is a vector giving the daily average temperature at the Portland meteorological station. Third component of the multivariate timeSeries TEMPS.ts}
	\item{Reno.ts}{The data slot is a vector giving the daily average temperature at the Reno meteorological station. First component of the multivariate timeSeries TEMPS1.ts}
	\item{LasVegas.ts}{The data slot is a vector giving the daily average temperature at the Las Vegas meteorological station. Second component of the multivariate timeSeries TEMPS1.ts}
	\item{Sacramento.ts}{The data slot is a vector giving the daily average temperature at the Sacramento meteorological station. Third component of the multivariate timeSeries TEMPS1.ts}
	\item{Newark.ts}{The data slot is a vector giving the daily average temperature at the Newark meteorological station. First component of the multivariate timeSeries TEMPS2.ts}
	\item{La Guardia.ts}{The data slot is a vector giving the daily average temperature at the La Guardia meteorological station. Second component of the multivariate timeSeries TEMPS2.ts}
	\item{Philadelphia.ts}{The data slot is a vector giving the daily average temperature at the Philadelphia meteorological station. Third component of the multivariate timeSeries TEMPS2.ts}
}
\references{
R. A. Carmona: Statistical Analysis of Financial Data in R, (2013) Springer Verlag
}
\examples{
plot(DesMoines.ts[1:1000,])
}
\keyword{datasets}
