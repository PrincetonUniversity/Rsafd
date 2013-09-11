########################################
### Defining Classes & Methods for time series ###
########################################

# source("/Users/rcarmona/Desktop/Rstuff/timeSeries.tex")

##########
# Defining class timeSeries

setClassUnion("numericORmatrix",c("numeric","matrix"))
setClass("timeSeries", representation(positions  = "timeDate", data = "numericORmatrix", units = "character", zone = "character", FinCenter = "character"))
   

# Methods to extract the different slots of a timeSeries

if (!isGeneric("seriesPositions")) 
{
	if (is.function("seriesPositions"))
		 fun <- seriesPositions
	else fun <- function(object) standardGeneric("seriesPositions")
	setGeneric("seriesPositions", fun)
}
setMethod("seriesPositions", "timeSeries", function(object) object@positions)

if (!isGeneric("seriesData")) 
{
	if (is.function("seriesData"))
		 fun <- seriesData
	else fun <- function(object) standardGeneric("seriesData")
	setGeneric("seriesData", fun)
}
setMethod("seriesData", "timeSeries", function(object) object@data)

if (!isGeneric("seriesUnits")) 
{
	if (is.function("seriesUnits"))
		 fun <- seriesUnits
	else fun <- function(object) standardGeneric("seriesUnits")
	setGeneric("seriesUnits", fun)
}
setMethod("seriesUnits", "timeSeries", function(object) object@units)


if(isGeneric("timeSeries")) removeGeneric("timeSeries")

timeSeries <- function (positions, data, units = NULL,  zone = "GMT", FinCenter = "GMT", ...) 
{
    if (missing(positions)) 
    {
        N <- dim(as.matrix(data))[1]
        positions <- timeSequence(from = "1960-01-01", length.out = N, zone = "GMT", FinCenter = "GMT")
    }
    if (is.character(positions)) 
    {
          format <- .whichFormat(positions)
       	 positions <- timeDate(charvec=positions, format=format, zone=zone, FinCenter=FinCenter)
    }
    timeDates <- positions
    data = as.matrix(data)
    rownames(data) = as.character(timeDates)
    if (is.null(units)) 
    {
            if (is.null(colnames(data))) 
                units = paste("TS.", 1:dim(data)[2], sep = "")
            else 
                units = colnames(data)
     }
     colnames(data) = units
    OUT = new("timeSeries", data = data, positions = timeDates, FinCenter = timeDates@FinCenter, zone = timeDates@FinCenter, units = as.character(units))
    return(OUT)
}
setGeneric("timeSeries")

if (!isGeneric("is.timeSeries")) 
{
	if (is.function("is.timeSeries"))
		 fun <- is.timeSeries
	else fun <- function(object) standardGeneric("is.timeSeries")
	setGeneric("is.timeSeries", fun)
}
setMethod("is.timeSeries", "timeSeries", function(object) is(object,"timeSeries"))

diff.timeSeries <- function (x, lag = 1, diff = 1, ...) 
{
# If diff=1, diff.timeSeries(x)[j]= x[j]-x[j-k]  if lag=k>=0, and x[j+|k|]-x[j] of lag=k<0
    y = as.matrix(seriesData(x))
    L <- length(seriesPositions(x))
    if (lag>=0)
    {
    	NEWDATA <- diff(y, lag=lag, differences=diff)
       Y <- timeSeries(positions=timeDate(rownames(NEWDATA)),data=NEWDATA,
              units = seriesUnits(x), zone = x@FinCenter, FinCenter = x@FinCenter)
    }
    else
    {
      k <- as.integer(abs(lag))
      NEWDATA <- diff(y, lag=k, differences=diff)
      Y <- timeSeries(positions=seriesPositions(x)[1:(L-k)],data=NEWDATA,
      	units = seriesUnits(x), zone = x@FinCenter, FinCenter = x@FinCenter)
    }
    Y
}
setMethod("diff","timeSeries", function(x, lag = 1, diff = 1, ...)  diff.timeSeries(x, lag = 1, diff = 1, ...) )

if (!isGeneric("plot"))
{
    if (is.function("plot"))
    fun <- plot
    else fun <- function(object) standardGeneric("plot")
    setGeneric("plot", fun)
}
setMethod("plot","timeSeries", function(x) {
    stack <- sys.calls( )
    stack.fun <- Filter( function(.) .[[1]] == as.name("plot"), stack )
    NAME <- paste("Time Series Plot of ",deparse( stack.fun[[1]][[2]] ),sep="")
    MAT <- x@data
    ROWNAMES <- as.character(rownames(MAT))
    DATES <- as.POSIXct(ROWNAMES,tz="GMT")
    DAT <- as.numeric(MAT)
    MINX <- min(DAT[!is.na(DAT)])
    MAXX <- max(DAT[!is.na(DAT)])
    p <- dim(MAT)[2]
    plot(DATES,MAT[ ,1],type="l",col=1,xaxs="i",lty=1,ylab="",ylim=c(MINX,MAXX))
    title(main=NAME,cex.main=0.8)
    if (p>1)
    {
        for (K in 2:p)        
        {
            points(DATES, MAT[ ,K], type="l", lty=2*K, col=K)
            legend("topright",colnames(MAT),col=(1:p),lty=1+2*(0:(p-1)),cex=.6,bty="n")
        }
    }    
    invisible(x)
}
)


#if (!isGeneric("ar"))
#{
#    if (is.function("ar")) fun <- ar
#    else fun <- function(object) standardGeneric("ar")
#    setGeneric("ar", fun)
#}

setMethod("ar",signature(x = "timeSeries"),  function(x,aic=TRUE,order.max=NULL,method=c("yule-walker", "burg", "ols", "mle", "yw"), na.action = na.fail, series = deparse(substitute(x)), ...) 
 {
#	if (class(x) != "timeSeries")
#        		stop("The function ar.timeSeries is restricted to timeSeries objects")
	stack <- sys.calls( )
    	stack.fun <- Filter( function(.) .[[1]] == as.name("ar.timeSeries"), stack )
    	NAME <- paste("Time Series Plot of ",deparse( stack.fun[[1]][[2]] ),sep="")
	TMP <- ar(x@data)
    	POS <- x@positions
	LL <- length(POS)
	POS <- POS[(TMP$order +1):LL]
	DATA <- as.matrix(TMP$resid)
	DATA <- DATA[(TMP$order +1):LL,]
    	TS <- timeSeries(positions = POS,data=DATA, units="Residuals", FinCenter=x@FinCenter, zone=x@zone)
	OUT <- list(order=TMP$order, ar=TMP$ar, var.pred=TMP$var.pred, x.mean=TMP$x.mean, aic=TMP$aic,partialacf=TMP$partialacf, resid =TS, method=TMP$method, series = NAME, call=TMP$call)
	return(OUT)
} 
)


if (!isGeneric("begday")) 
{
	if (is.function("begday"))
		 fun <- begday
	else fun <- function(object) standardGeneric("begday")
	setGeneric("begday", fun)
}
setMethod("begday","timeDate", function(object) {trunc(object,units="days")})

head.timeSeries <- function(x) {utils::head(x@data)}
setMethod("head","timeSeries", function(x) head.timeSeries(x))
tail.timeSeries <- function(x) {utils::tail(x@data)}
setMethod("tail","timeSeries", function(x) tail.timeSeries(x))

#setGeneric("head", function(x) standardGeneric("head"))
#head.timeSeries <- function(x) {head(x@data)}
#setMethod("head","timeSeries", function(x) head.timeSeries(x))
#setGeneric("head", function(x) standardGeneric("head"))
#head.timeSeries <- function(x) {UseMethod("head",x@data)}
#tail.timeSeries <- function(x) {UseMethod("tail",x@data)}
#tail.timeSeries <- function(x) {tail(x@data)}


setMethod("Math", "timeSeries",
            function(x) {
			x@data <- callGeneric(x@data)
			slot(x,"units") <- "TBA"
			x
})

setMethod("Arith", c("timeSeries","numeric"),
 	function(e1, e2) {
  		e1@data <- callGeneric(e1@data , e2)
		slot(e1,"units") <- "TBA"
  		e1
})

setMethod("Arith", c("numeric","timeSeries"),
 	function(e1, e2) {
  		e2@data <- callGeneric(e1, e2@data)
		slot(e2,"units") <- "TBA"
  		e2
})

setMethod("Ops", c("timeSeries","timeSeries"),
            function(e1,e2) {
            	e1@data <- callGeneric(e1@data,e2@data)
			slot(e1,"units") <- "TBA"
            	e1
})

#If this works, it will have to be fixed for discrepancies in the positions
setMethod("Arith", c("timeSeries","timeSeries"),
            function(e1,e2) {
            	e1@data <- callGeneric(e1@data,e2@data)
			slot(e1,"units") <- "TBA"
            	e1
})

setMethod("[", "timeSeries",
function(x, i,j ) {
	POS <- as.character(x@positions)[i]
	TMP <- as.matrix(x@data[i,j])
#	rownames(TMP) <- POS
	TMP <- timeSeries(positions=timeDate(POS,zone=x@zone,FinCenter=x@FinCenter),data=TMP,units=x@units[j], FinCenter=x@FinCenter,zone=x@zone)
	TMP
})

setGeneric("concat", function(x,y ) standardGeneric("concat"))
setMethod("concat", c("timeSeries","timeSeries"),
            function(x,y) {
# Concatenate two timeSeries
#	if (!is.timeSeries(x) | !is.timeSeries(x))
#		stop("Arguments of c.timeSeries should be objects of class timeSeries")
	if ( dim(seriesData(x))[2] != dim(seriesData(y))[2])
		stop("arguments x and y should be timeSeries with the same numbers of variables")
	if ( seriesUnits(x) != seriesUnits(y))
		stop("arguments x and y should have the same units")
	NewData <- rbind(seriesData(x),seriesData(y))
	return(timeSeries(positions=timeDate(rownames(NewData)),data=NewData))
})

#c.timeSeries <- function(x,y)
#{
# Concatenate two timeSeries
#	if (!is.timeSeries(x) | !is.timeSeries(x))
#		stop("Arguments of c.timeSeries should be objects of class timeSeries")
#	if ( dim(seriesData(x))[2] != dim(seriesData(y))[2])
#		stop("arguments x and y should be timeSeries with the same numbers of variables")
#	if ( seriesUnits(x) != seriesUnits(y))
#		stop("arguments x and y should have the same units")
#	NewData <- rbind(seriesData(x),seriesData(y))
#	return(timeSeries(positions=timeDate(rownames(NewData)),data=NewData))
#}
#setMethod("c",c("timeSeries","timeSeries"), function(x,y) {c.timeSeries(x,y)})

#if (!isGeneric("concat")) 
#{
#	if (is.function("concat"))
#		 fun <- concat
#	else fun <- function(x,y) standardGeneric("concat")
#	setGeneric("concat", fun)
#}
#setMethod("concat",signature=c("timeSeries","timeSeries"), function(x,y) {c.timeSeries(x,y)})

#######################################
# Example from class copula (to be eventually removed)
# if(isGeneric("contour.pcopula")) removeGeneric("contour.pcopula")
# setGeneric("contour.pcopula", function(cop,n=50, ... ) standardGeneric("contour.pcopula"))
# setMethod("contour.pcopula","copula",
# 	function(cop, n=50, ... )
# 	{
#     		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
#     		divis <- divis[2:(n+1)]
#     		xmat <- rep(divis, n )
#   		ymat <- rep(divis, each = n)
#    		zmat <- pcopula(cop,  xmat,  ymat)
#    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
#    		contour(val)
#    		title("Contour Plot of the Copula", xlab="u",ylab="v")
#    		invisible(val)
#	})
#########################################
if(isGeneric("merge.timeSeries")) removeGeneric("merge.timeSeries")
setGeneric("merge.timeSeries", function(x,y) standardGeneric("merge.timeSeries"))
merge.timeSeries <- function (x, y, INTERSECT = FALSE, units = NULL, ...)
{
    if (INTERSECT == TRUE)
    {
        x.data <- seriesData(x)
        x.names <- rownames(x.data)
        y.data <- seriesData(y)
        y.names <- rownames(y.data)
        xy <- intersect(x.names,y.names)
        xy <- sort(xy)
        xy.data <- cbind(x.data[xy,],y.data[xy,])
        rownames(xy.data) <- xy
        ans <- timeSeries(positions=timeDate(xy), data=xy.data, units=c(x@units,y@units))
    }
    else
    {
        POS <- as.character(c(as.character(x@positions), as.character(y@positions)))
        LENGTH = length(as.character(seriesPositions(x)))
        DUP = duplicated(POS)[1:LENGTH]
        DUP2 = duplicated(POS)[-(1:LENGTH)]
        M1 = seriesData(x)
        M2 = seriesData(y)
        dim1 = dim(M1)
        dim2 = dim(M2)
        X1 = matrix(rep(NA, times = dim1[1] * dim2[2]), ncol = dim2[2])
        X2 = matrix(rep(NA, times = dim2[1] * dim1[2]), ncol = dim1[2])
        colnames(X1) = colnames(M2)
        NC = (dim1 + dim2)[2] + 1
        Z = rbind(cbind(M1, X1, DUP), cbind(X2, M2, DUP2))
        Z = Z[order(rownames(Z)), ]
        NC1 = dim1[2] + 1
        IDX = (1:(dim1 + dim2)[1])[Z[, NC] == 1]
        Z[IDX - 1, NC1:(NC - 1)] = Z[IDX, NC1:(NC - 1)]
        Z = Z[!Z[, NC], -NC]
        timeDates <- timeDate(charvec=rownames(Z), zone=x@zone, FinCenter=x@FinCenter)
        ans <- timeSeries(positions = timeDates, data = Z, zone = x@zone,
        FinCenter = x@FinCenter, units = c(seriesUnits(x), seriesUnits(y)))
        if (!is.null(units)) {
            seriesUnits(ans) <-  units
            colnames(ans@Data) <- units
        }
        ans
    }
}
setMethod("merge",c("timeSeries","timeSeries"), function(x,y) merge.timeSeries(x,y))


if(isGeneric("sstl")) removeGeneric("sstl")
setGeneric("sstl", function(SERIES, FREQ=365, TWIND = 0.05,TDEGREE=0) standardGeneric("sstl"))
setMethod("sstl","timeSeries", function (SERIES, FREQ=365, TWIND = 0.05,TDEGREE=0) 
{
    if (class(SERIES) != "timeSeries") 
        stop("The function sstl is restricted to timeSeries objects")
    TMP.ts <- ts(as.numeric(seriesData(SERIES)), start=1, frequency=FREQ)
    TMP.stl <- stl(TMP.ts, s.window="periodic",t.degree=TDEGREE,t.window=TWIND*length(TMP.ts))
    TREND <- timeSeries(positions = seriesPositions(SERIES), data = TMP.stl$time.series[,2], units = paste(SERIES@units,"_trend",sep=""))
    SEA <- timeSeries(positions = seriesPositions(SERIES), data = TMP.stl$time.series[,1], units = paste(SERIES@units,"_sea",sep=""))
    REM <- timeSeries(positions = seriesPositions(SERIES), data = TMP.stl$time.series[,3], units = paste(SERIES@units,"_rem",sep=""))
    OUT <- list(trend = TREND, sea = SEA, rem = REM) 
    return(OUT)
})

#setGeneric("as.ts", function(x) standardGeneric("as.ts"))
#setMethod("as.ts","timeSeries", function(x) {return(ts(data=x@data,start=1,end=length(x@positions)))})

#setOldClass("ts", S4Class = "timeSeries")
#selectMethod("ar","ts")
#setAs(from="timeSeries", to="ts", function(from) {ts(data=from@data,start=1,end=length(from@positions))})

#setGeneric("arima", function(x,order) standardGeneric("arima"))
#setMethod("arima","timeSeries", function(x,order) arima(as(x,"ts"),order))
#setGeneric("ar", function(x,aic,order.max,method) standardGeneric("ar"))
#setMethod("ar","timeSeries", function(x,aic,order.max,method) ar(as(x,"ts"),aic=FALSE,order.max,method))
#setGeneric("acf", function(x,y) standardGeneric("acf"))
#setMethod("acf","timeSeries", function(x,y) acf.timeSeries(x,type=y))

as.ts.timeSeries <- function(x) {ts(data=x@data,start=1,end=length(x@positions))}




arima.timeSeries <- function(x,order) {arima(x@data,order)}
acf.timeSeries <- function(x,y)
{
	TST <- acf(x@data,type=y, plot=FALSE)
	plot(TST, main=paste(y," - acf of ",deparse(substitute(x)),sep=""),cex.main=0.8)
	return(TST)
}

timeSeriesFromList <- function(LIST)
{
    return(timeSeries(positions=LIST$pos,data=LIST$data))
}


################################
# Two methods for the timeDate class

if (!isGeneric("begday")) 
{
	if (is.function("begday"))
		 fun <- begday
	else fun <- function(object) standardGeneric("begday")
	setGeneric("begday", fun)
}
setMethod("begday","timeDate", function(object) {trunc(object,units="days")})

if (!isGeneric("noon")) 
{
	if (is.function("noon"))
		 fun <- noon
	else fun <- function(object) standardGeneric("noon")
	setGeneric("noon", fun)
}
setMethod("noon","timeDate", function(object) timeDate(paste(trunc(object,units="days"), "12:00:00"), format = "%Y-%m-%d %H:%M:%S"))

#####################################
### Utilities

makeDate <- function(x,in.format="d-m-y")
{
# I assume that x a vector of character strings representing dates
# The argument in.format gives the format of the entries of x
# in.format = c("Ymd","Y-m-d","y-m-d","Y/m/d","y/m/d","d-m-Y","d-m-y","d/m/Y","d/m/y","m/d/y","m/d/Y","m-d-y","m-d-Y")
# where Y is a string with four digits, y and d are one or two digit strings 
# representing integers, and m is either a one or two digit string 
# (like 7 or 12) or a three character string like Apr
# Returns a vector of character strings of lengths 10, xxxx-xx-xx for Y-m-d 
# like 1997-01-07  
# Usage: makeDate(DD,format="d-m-y")
	if (in.format == "Ymd")
	{
		y <- substr(x,1,4)
		m <- substr(x,5,6)
		Nd <- substr(x,7,8)
	}
	else
	{	
		x <- gsub("/","-",x)
		x <- strsplit(x,"-")
		n <- length(x)
		x <-unlist(x)
		if (substr(in.format,1,1)=="m") 
		{
			m <- as.integer(x[1+3*(0:(n-1))])
			if (substr(in.format,3,3)=="d") 
			{
				d <- as.integer(x[2+3*(0:(n-1))])
				tmp <- as.integer(x[3*(1:n)])
				if (substr(in.format,5,5)=="y") 
				{
					Ntmp <- rep(0,length(tmp))
					if (sum(tmp<40)>0) Ntmp[tmp<40] <- 2000+tmp[tmp<40]
					if (sum(tmp>=40)>0) Ntmp[tmp>=40] <- 1900+tmp[tmp>=40]
					y <- Ntmp
				}
				if (substr(in.format,5,5)=="Y") 
					y <- tmp
			}
			if (substr(in.format,3,3)=="y") 
			{
				Ntmp <- rep(0,length(tmp))
				if (sum(tmp<40)>0) Ntmp[tmp<40] <- 2000+tmp[tmp<40]
				if (sum(tmp>=40)>0) Ntmp[tmp>=40] <- 1900+tmp[tmp>=40]
				y <- Ntmp
				d <- as.integer(x[3*(1:n)])
			}
			if (substr(in.format,3,3)=="Y") 
			{
				y <- as.integer(x[2+3*(0:(n-1))])
				d <- as.integer(x[3*(1:n)])
			}
		}
		else
		{
			m <- x[2+3*(0:(n-1))]
			if (sum(is.na(as.integer(m)))!=0)
			{
				m[m=="Jan"] <- "01"
				m[m=="Feb"] <- "02"
				m[m=="Mar"] <- "03"
				m[m=="Apr"] <- "04"
				m[m=="May"] <- "05"
				m[m=="Jun"] <- "06"
				m[m=="Jul"] <- "07"
				m[m=="Aug"] <- "08"
				m[m=="Sep"] <- "09"
				m[m=="Oct"] <- "10"
				m[m=="Nov"] <- "11"
				m[m=="Dec"] <- "12"
			}
			tmp <- as.integer(x[1+3*(0:(n-1))])
			if (substr(in.format,1,1)=="Y") 
			{
				y <- tmp
				d <- x[3*(1:n)]
			}
			if (substr(in.format,1,1)=="y")
			{
				Ntmp <- rep(0,length(tmp))
				if (sum(tmp<40)>0) Ntmp[tmp<40] <- 2000+tmp[tmp<40]
				if (sum(tmp>=40)>0) Ntmp[tmp>=40] <- 1900+tmp[tmp>=40]
				y <- Ntmp
				d <- as.integer(x[3*(1:n)])
			}
			if (substr(in.format,5,5)=="Y") 
			{
				d <- tmp
				y <- x[3*(1:n)]
			}
			if (substr(in.format,5,5)=="y")
			{
				Ntmp <- rep(0,length(tmp))
				d <- tmp
				tmp <- as.integer(x[3*(1:n)])
				if (sum(tmp<40)>0) Ntmp[tmp<40] <- 2000+tmp[tmp<40]
				if (sum(tmp>=40)>0) Ntmp[tmp>=40] <- 1900+tmp[tmp>=40]
				y <- Ntmp
			}
		}
		Nd <- rep("00",length(d))
		if (sum(d<10)>0) Nd[d<10] <- paste("0",d[d<10],sep="")
		if (sum(d>=10)>0) Nd[d>=10] <- d[d>=10]
	}
	return(paste(y,m,Nd,sep="-"))
}

dump.timeSeries <- function (TS, FILENAME) 
{
    if (!is.timeSeries(TS)) 
        stop("This function is restricted to timeSeries objects")
    DD <- list(pos = as.character(seriesPositions(TS)), data = seriesData(TS))
    OBJECTNAME <- deparse(substitute(TS))
    assign(OBJECTNAME, DD)
    dump(OBJECTNAME, file = FILENAME)
    invisible(DD)
}
	
	
#####################################
# Standard time Series stuff

kalman <- function (FF, SigV, GG, SigW, Xhat, Omega, Y) 
{
    Delta <- GG %*% Omega %*% t(GG) + SigW
    Theta <- FF %*% Omega %*% t(GG)
    X <- FF %*% Xhat + Theta %*% solve(Delta) %*% (Y - GG %*% 
        Xhat)
    Om <- FF %*% Omega %*% t(FF) + SigV - Theta %*% solve(Delta) %*% 
        t(Theta)
    Ret <- list(xpred = X, error = Om)
    Ret
}

pred.ar <- function (series, ar.est, ahead = 1) 
{
    order <- ar.est$order
    series <- as.matrix(series)
    pred.out <- array(NA, dim = c(order + ahead, ncol(series)), 
        dimnames = list(NULL, dimnames(series)[[2]]))
    mean.ser <- apply(series, 2, mean)
    ser.cent <- sweep(series, 2, mean.ser)
    pred.out[seq(order), ] <- ser.cent[rev(nrow(series) - seq(order) + 
        1), ]
    for (i in (order + 1):nrow(pred.out)) {
        pred.out[i, ] <- apply(aperm(ar.est$ar, c(1, 3, 2)) * 
            as.vector(pred.out[i - seq(order), ]), 3, sum)
    }
    sweep(pred.out[-seq(order), , drop = FALSE], 2, mean.ser, "+")
}


garch <- function (x, order = c(1, 1), series = NULL, control = garch.control(...),  ...) 
{
    	stack <- sys.calls( )
    	stack.fun <- Filter( function(.) .[[1]] == as.name("garch"), stack )
    	NAME <- deparse( stack.fun[[1]][[2]] )
	TS <- (class(x) == "timeSeries")
	if (TS)
	{
		POS <- x@positions
		x <- x@data
	}
   	 if (NCOL(x) > 1) 
        		stop("Argument x is neither a vector nor a univariate time series")
    	if (any(is.na(x))) 
        		stop("NAs are present")
	OUT <- tseries::garch(x, order, series, control, ...)
	LL <- length(x)
	P <- order[1]
	if (TS)
		RES.ts <- timeSeries(positions=POS[(P+1):LL], data=OUT$residuals[(P+1):LL], units="Residuals")
	else RES.ts <- OUT$residuals[(P+1):LL]
	if (TS)
		SIGt.ts <- timeSeries(positions=POS[(P+1):LL], data=OUT$fitted.values[(P+1):LL,1], units="sig_t")
	else SIGt.ts <- OUT$fitted.values[(P+1):LL,1]
	FITTED <- OUT$fitted.values[(P+1):LL,1]*OUT$residuals[(P+1):LL]
	if (TS)
		FITTED.ts <- timeSeries(positions=POS[(P+1):LL], data=FITTED, units="hat_x_t")
	else FITTED.ts <- FITTED
	OUTPUT <- list(order=OUT$order, coef=OUT$coef, n.likeli = OUT$n.likeli, n.used = OUT$n.used, 
        resid = RES.ts, fitted.values = FITTED.ts, sigma.t = SIGt.ts, vciv = OUT$vcov, name = NAME)
}        
        
        
DF.test <- function (x, alternative = c("stationary", "explosive"), k ) 
{
    	stack <- sys.calls( )
    	stack.fun <- Filter( function(.) .[[1]] == as.name("DF.test"), stack )
    	NAME <- deparse( stack.fun[[1]][[2]] )
	if (class(x) == "timeSeries") 
		x <- x@data
   	 if (NCOL(x) > 1) 
        		stop("Argument x is neither a vector nor a univariate time series")
    	if (any(is.na(x))) 
        		stop("NAs are present")
    	if (missing(k))
		k = trunc((length(x) - 1)^(1/3))
   	 if (k < 0) 
        		stop("k negative")
    	alternative <- match.arg(alternative)
	k <- k + 1
    	y <- diff(x)
    	n <- length(y)
    	z <- embed(y, k)
    	yt <- z[, 1]
    	xt1 <- x[k:n]
    	tt <- k:n
    	if (k > 1) 
	{
        		yt1 <- z[, 2:k]
        		res <- lm(yt ~ xt1 + 1 + tt + yt1)
    	}
    	else res <- lm(yt ~ xt1 + 1 + tt)
    	res.sum <- summary(res)
    	STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2, 2]
    	table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), c(3.95, 
        3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 3.45, 3.43, 
        3.42, 3.41), c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), c(1.14, 
        1.19, 1.22, 1.23, 1.24, 1.25), c(0.8, 0.87, 0.9, 0.92, 
        0.93, 0.94), c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15, 
        0.24, 0.28, 0.31, 0.32, 0.33))
    	table <- -table
    	tablen <- dim(table)[2]
    	tableT <- c(25, 50, 100, 250, 500, 1e+05)
    	tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    	tableipl <- numeric(tablen)
    	for (i in (1:tablen)) 
		tableipl[i] <- approx(tableT, table[, i], n, rule = 2)$y
    	interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
   	 if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y)) 
        		if (interpol == min(tablep)) 
            		warning("p-value smaller than printed p-value")
        		else warning("p-value greater than printed p-value")
   	 if (alternative == "stationary") 
        		PVAL <- interpol
    	else if (alternative == "explosive") 
        			PVAL <- 1 - interpol
   		else stop("irregular alternative")
    	PARAMETER <- k - 1
    	METHOD <- "Augmented Dickey-Fuller Test"
    	names(STAT) <- "Dickey-Fuller"
    	names(PARAMETER) <- "Lag order"
    	list(statistic = STAT, parameter = PARAMETER, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = NAME)
}

sim.garch <- function (model, n.ahead = 1024, n.start = 0, n.sim = 1, presigmat = NA, preXt = NA , innovations = NA )
{
	# model is assumed to be a list c(sig2,ar,ma) with 
	# sig a constant, ar the vector of length p of ar coefs and ma the vector of length q of ma coefs
	
	if ((!missing(presigmat)) && (!missing(preXt))) 
		n.start <- 0
    	if (!is.list(model)) 
        		stop("'model' should be a list")
	
    	p <- length(model$ar)
	if( !missing(presigmat) && length(presigmat)!=p )
		stop("presigmat should be a vector of length the order of the ar part")
    	if (p) 
	{
        		minroots <- min(Mod(polyroot(c(1, -model$ar))))
        		if (minroots <= 1) 
           		 stop("'ar' part of model is not stationary")
    	}
    	q <- length(model$ma)
	if( !missing(preXt) && length(preXt)!=q )
		stop("preXt should be a vector of length the order of the ma part")
	if (missing(presigmat) || missing(preXt))
	{
		presigmat <- rep(0,p)
		preXt <- rep(0,q)
	}
    	if (!missing(innovations) && NROW(innovations) != n.sim) 
		stop("the matrix of innovations does not have the right number of rows")
    	if (!missing(innovations) && NCOL(innovations) < (n.start+n.ahead)) 
		stop("the matrix of innovations does not have enough columns")
   	if(missing(innovations)) INNOV <- array(rnorm(n.sim*(n.start + n.ahead)),dim=c(n.sim,n.start + n.ahead))
    	else INNOV <- innovations[,1:(n.start + n.ahead)]
    	AR <- array(rev(model$ar), dim=c(p,1))
	MA <- array(rev(model$ma), dim=c(q,1))
	SIG2 <- rep(model$sig2,n.sim)
    	Lsig <- p + n.start+n.ahead
    	Lx <- q + n.start+n.ahead
    	SIG2t <- array(0,dim=c(n.sim,Lsig))
	SIG2t[,1:p] <- outer(rep(1,n.sim),presigmat^2,"*")
    	Xt <- array(0,dim=c(n.sim,Lx))
	Xt[,1:q] <- outer(rep(1,n.sim),preXt,"*")
    	X2t <- array(0,dim=c(n.sim,Lx))
	X2t[,1:q] <- outer(rep(1,n.sim),preXt^2,"*")
    	for (J in 1:(n.start+n.ahead))
	{
		SIG2t[,p+J] <- SIG2 + SIG2t[,J:(p+J-1)]%*%AR + X2t[,J:(q+J-1)]%*%MA
		Xt[,q+J] <- sqrt(SIG2t[,p+J])*INNOV[,J]
		X2t[,q+J] <- Xt[,q+J]^2
	}
    	OUT <- list(X.t = Xt[,(q+n.start+1):Lx], sigma.t=sqrt(SIG2t[,(p+n.start+1):Lsig]))
}
