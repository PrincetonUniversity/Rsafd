##################################
### Defining Classes & Methods ###
##################################

SHAPE.XI <- TRUE

##########
# Defining class gpd

setClassUnion("numericORlogical",c("numeric","logical"))

setClass("gpd.upper.tail", representation(
		tail = "character",
            n = "numeric",   
            data = "numeric",
            upper.exceedances = "numeric", 
            upper.threshold = "numeric",
            p.less.upper.threshold= "numeric", 
            n.upper.exceedances = "numeric", 
            upper.method = "character", 
            upper.par.ests = "numeric", 
            upper.converged = "logical" )
)            
             

setClass("gpd.lower.tail", representation(
		tail = "character",
            n = "numeric",   
            data = "numeric",
            lower.exceedances = "numeric", 
            lower.threshold = "numeric", 
            p.larger.lower.threshold =  "numeric", 
            n.lower.exceedances = "numeric", 
            lower.method = "character",
            lower.par.ests = "numeric", 
            lower.converged = "logical") 
)            
 
setClass("gpd.two.tails", representation(
		tail = "character",
            n = "numeric",   
            data = "numeric",
            upper.exceedances = "numeric", 
            lower.exceedances = "numeric", 
            upper.threshold = "numeric",
            lower.threshold = "numeric", 
            p.less.upper.threshold= "numeric", 
           p.larger.lower.threshold =  "numeric", 
            n.upper.exceedances = "numeric", 
            n.lower.exceedances = "numeric", 
            upper.method = "character", 
            lower.method = "character",
            upper.par.ests = "numeric", 
            lower.par.ests = "numeric", 
            upper.converged = "logical", 
            lower.converged = "logical") 
)            
 
#setClass("gpd", contains=c("gpd.upper.tail", "gpd.lower.tail","gpd.two.tails"))
setClassUnion("gpd", c("gpd.upper.tail", "gpd.lower.tail","gpd.two.tails"))

SHAPE.XI <- TRUE
   
############
## Remove methods if they already exist
  
 if(isGeneric("tailplot")) removeGeneric("tailplot")




####################
####### GEV stuff

dgev <- function(x,  m=0, lambda = 1, xi = 0) 
{    
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(x) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as x")
	val <- rep(Inf, length(x))  
	uu <- exp(-(x-m)/lambda)
	vv <- (1+(x-m)*k/lambda)^(-1/k -1)
	if (n==1)
	{
		if (k==0)
            val <- uu/lambda*exp(-uu) 
		else
			val <-   vv/lambda*exp(-vv*(1+(x-m)*k/lambda))
	}
	else
	{
		K0 <- (k==0)
		KN0 <- !K0
		val[K0] <- uu[K0]/lambda[K0]*exp(-uu[K0])
		val[KN0] <- vv[KN0]/lambda[KN0]*exp(-vv[KN0]*(1+(x[KN0]-m[KN0])*k[KN0]/lambda[KN0]))
	}
    val
}

pgev <- function(q, m=0, lambda = 1, xi = 0)
{   
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(q) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as q")
	val <- rep(Inf, length(q))  
	uu <- exp(-exp(-(q - m)/lambda)) 
	vv <- exp( - (1. + (k * (q - m))/lambda)^(-1./k)) 
	if (n==1)
	{
		if (k==0)
            val <- uu 
		else
			val <- vv
	}
	else
	{
		K0 <- (k==0)
		KN0 <- !K0
		val[K0] <- uu[K0]
		val[KN0] <- vv[KN0]
	}
    val
}

qgev <- function(p,  m=0, lambda = 1, xi = 0) 
{
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(p) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as p")
    if (sum(p <= 0 | p >= 1)>0) 
        warning("Argument of qpareto should be between 0 and 1, NA will be returned in other cases")
	val <- rep(Inf, length(p))  
	if (n==1)
	{
		if (k==0)  
			val <- m - lambda * log(-log(p))
		else 
			val <- m - lambda/k * ( 1 - (-log(p))^(-k))
	}
	else
	{
		K0 <- (k==0)
		val[K0] <- m[K0] - lambda[K0] * log(-log(p[K0]))
		val[!K0] <- m[!K0] - lambda[!K0]/k[!K0] * ( 1 - (-log(p[!K0]))^(-k[!K0]))
	}
    val [p < 0 | p > 1] <- NA
    val
}

rgev <- function(n ,  m=0, lambda = 1, xi = 0) 
{
   qgev(runif(n),m,lambda,xi)
}



#######################
##Estimation from data

gev.lmom <- function(lmom)
{

        if (length(lmom) > 4) {
            message <- "It looks like the argument of gev.lmom is a sample and not sample L-moments."
            message <- paste(message, " Parameter estimation proceeds with this assumption,")
            message <- paste(message, " computing the L-moments from this sample and then,")
            message <- paste(message, " applying the function gev.lmom to this set of L-moments")
            warning(message)
            lmom <- sample.LMOM(lmom)
        }
    
		epsilon <- 1e-6;
       
        c1 <- 7.817740; c2 <- 2.930462; c3 <- 13.641492; c4 <- 17.206675;
        eu <- -0.577215664901532861;
        z0 <- log(2.0)/log(3.0);

        t3 <- lmom[3];
        
        if(is.na(lmom[3]) | lmom[3] > 1)
            { stop(" The L-skewness is not specified or greater than 1 " ); }
        if(is.na(lmom[2]) | lmom[2] < 0 )
            { stop("Negative second L-moment!" ); }
       
        z <- 2.0/(3.0+t3) - z0;
         
        # Initial guess for k
        g <- z*(c1+z*(c2 + z*(c3 + z*c4)));

        # don't solve the equation, since the approximation
        # is good for  -0.1 < k < 0.5
        dosolve <- TRUE        
        if ( t3 >= -0.1 & t3 <= 0.5) {dosolve <- FALSE}
        if (dosolve & t3 < -0.9) {g <- 1.0 - log(1.0+t3)/log(2.0)}
        
        t0 <- 1;
        t00 <- (t3 + 3.0)/2.0
        assign("t0",t00, pos = 1)
        
        equation.error <- function(g) {
              t <- (1.0 - 3.0^(-g))/(1.0 - 2.0^(-g)) - t0
              abs(t)
        }
    
        if (dosolve)
        {   
           fit <- nlm(equation.error, g, iterlim = 100, steptol = epsilon)
           if ((fit$code==3)|(fit$code==4)|(fit$code==5))
                warning("LMOM estimation for GEV: Solution of the equation for k did not converge")
           g <- fit$estimate
        }
        paramest <- c(NA,NA,NA)
        names(paramest) <- c("m","lambda","k")
       
        paramest[3] <- g;
        gam <- gamma(1.0+g);

        # I changed the original order of EVANESCE: location is now first, scale second           
        paramest[2] <- lmom[2] *g / (gam*(1.0-2.0^(-g)));
        paramest[1] <- lmom[1] - paramest[2] *(1.0 - gam)/g;
        
        if (SHAPE.XI) 
        { 
            paramest[3] <- -paramest[3]
            names(paramest) <-  c("m","lambda","xi")
        }
        
        val <- list(param.est = paramest)
        val
}

gev.ml <- function(sample, init.est = NA, epsilon = 1e-5) 
{
        lmomest <- init.est
        assign("tempX",sample,pos= 1)
        assign("tempN",length(sample), pos =1)

        if (is.na(lmomest[1])) 
        { 
            lmomest <- gev.lmom(sample.LMOM(tempX))$param.est
            if (SHAPE.XI)  lmomest[3] <- -lmomest[3]
        }
        negative.log.likelihood <- function(theta)
        {
            k <- theta[3]
            m <- theta[1]
            lambda <- theta[2]

            xsc <- 1 - k*(tempX-m)/lambda
            ll <- NA 
            if (sum(xsc < 0) > 0 | lambda < 0) { ll <- NA}
            else  ll <- -tempN*log(lambda) - sum(xsc^(1/k)) + (1/k - 1)*sum(log(xsc))
            
            -ll   
        }  
        
        fit <- nlm(negative.log.likelihood, lmomest, iterlim = 100, steptol = epsilon)
        if (fit$code>3) 
            warning("Maximum Likelihood Method for GEV did not converge")

        paramest <- fit$estimate
        names(paramest) <- c("m","lambda","k")
        if (SHAPE.XI)
        {
            paramest[3] <- -paramest[3]
            names(paramest) <-  c("m","lambda","xi")
        }

        val <- list(param.est = paramest, converged = as.logical(fit$code<4))
        val
}

sample.LMOM <- function(x) 
{
          x <- sort(x)
          N <- length(x)
          i <- c(1:N)

    
          fn1 <- N - i
          fn2 <- (N - i - 1)*fn1
          fn3 <- (N - i - 2)*fn2
          
          a1 <- sum(fn1/(N-1) * x)/N
          a2 <- sum(fn2/(N-1)/(N-2) * x)/N
          a3 <- sum(fn3/(N-1)/(N-2)/(N-3) * x)/N

          l1 <- mean(x)
          l2 <- l1 - 2*a1
          tau2 <- (l1 - 6.0*a1 + 6.0*a2)/l2
          tau3 <- (l1 - 12.0*a1 + 30.0*a2 - 20.0*a3)/l2   

          val <- c(l1,l2,tau2,tau3)
          names(val) <- c("ell_1", "ell_2", "tau_3", "tau_4")        
          
          val
}

plotting.positions <- function(sample, gamma = -0.35, delta = 0 ) 
{
       x <- sort(sample)
       n <- length(x)

       pp.pos <- function(i,n, gamma, delta) {
           (i + gamma)/(n + delta)
        }
       
       pp <- pp.pos(c(1:n),n,gamma,delta)
       ell1 <- sum(pp*x)/n
       ell2 <- sum((2*pp-1)*x)/n
       ell3 <- sum((6*pp*pp - 6*pp+1)*x)/n
       ell4 <- sum((20*pp^3-30*pp*pp +12*pp-1)*x)/n
      
       tau3 <- ell3/ell2
       tau4 <- ell4/ell2

       val <- c(ell1,ell2,tau3,tau4)
       names(val) <- c("ell_1", "ell_2", "tau_3", "tau_4")       
       val
}



######################################
########## GPD stuff

dpareto <- function(x, m = 0, lambda = 1, xi = 0) 
{
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(x) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as x")
	val <- rep(Inf, length(x))  
	if (n==1)
	{
		if (k==0) 
		{
			val <- exp(-(x - m)/lambda)/lambda
			val[x < m] <- 0
		}
		else 
		{ 
			val <- 1. + k * (x - m)/lambda
			if (k>0) 
			{
				val[x<=m] <- 0
				val[x>m] <- val[x>m]^(-1/k-1)/lambda
			}
			else 
			{
				val[x<=m | x >= m-lambda/k] <- 0
				val[x>m & x < m-lambda/k] <- val[x>m & x < m-lambda/k]^(-1/k-1)/lambda
			}
		}
	}
	else
	{
		K0 <- (k==0)
		KPOS <- (k>0)
		KNEG <- (k<0)
		val[K0] <- exp(-(x[K0] - m[K0])/lambda[K0])/lambda[K0]
		val[(x < m) & K0] <- 0
		val[!K0] <- 1. + k[!K0] * (x[!K0] - m[!K0])/lambda[!K0]
		val[KPOS & (x<=m)] <- 0
		val[KPOS & (x>m)] <- val[(!K0) & (x>m)]^(-1/k[(!K0) & (x>m)]-1)/lambda[(!K0) & (x>m)]
		val[KNEG &(x<=m | x >= m-lambda/k)] <- 0
		val[KNEG & x>m & x < m-lambda/k] <- val[KNEG & x>m & x < m-lambda/k]^(-1/k[KNEG & x>m & x < m-lambda/k]-1)/lambda[KNEG & x>m & x < m-lambda/k]
	}
    val
}

ppareto <- function(q, m = 0, lambda = 1, xi = 0)
{   
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(q) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as q")
	val <- rep(Inf, length(q))  
	if (n==1)
	{
		if (k==0) {
			val <- 1. - exp(-(q - m)/lambda)
			val[q<m] <- 0
		}
		else 
		{ 
			val <- 1. + k * (q - m)/lambda
			if (k>0) {
				val[q<=m] <- 0
				val[q>m] <- 1 - val[q>m]^(-1/k)
			}
			else {
				val[q<=m | q >= m-lambda/k] <- 0
				val[q>m & q < m-lambda/k] <- 1 - val[q>m & q < m-lambda/k]^(-1/k)
			}
		}
	}
	else
	{
		K0 <- (k==0)
		KPOS <- (k>0)
		KNEG <- (k<0)
		val[K0] <- 1. - exp(-(q[K0] - m[K0])/lambda[K0])
		val[(q < m) & K0] <- 0
		val[!K0] <- 1. + k[!K0] * (q[!K0] - m[!K0])/lambda[!K0]
		val[KPOS & (q<=m)] <- 0
		val[KPOS & (q>m)] <- 1 - val[KPOS & (q>m)]^(-1/k[KPOS & (q>m)])
		val[KNEG &(q<=m | q >= m-lambda/k)] <- 0
		val[KNEG & q>m & q < m-lambda/k] <- 1 - val[KNEG & q>m & q < m-lambda/k]^(-1/k[KNEG & q>m & q < m-lambda/k])
	}
    val
}

qpareto <- function(p, m = 0, lambda = 1, xi = 0) 
{
    k <- xi
    if (!SHAPE.XI) k <- -xi
	if ((length(m)!=length(lambda))|(length(m)!=length(xi))|(length(xi)!=length(lambda)))
		stop("m, lambda and xi should have the same lengths")
	n <- length(m)
	if ((n > 1) & (length(p) != n))
		stop("When vectors of lengths > 1, m, lambda and xi should have the same length as p")
    if (sum(p <= 0 | p >= 1)>0) 
        warning("Argument of qpareto should be between 0 and 1, NA will be returned in other cases")
	val <- rep(Inf, length(p))  
	if (n==1)
	{
		if (k==0)  
			val <- m - lambda*log(1 - p)
		else 
			val <- m - lambda*( 1 - (1 - p)^(-k))/k
	}
	else
	{
		K0 <- (k==0)
		val[K0] <-	m[K0] - lambda[K0]*log(1 - p[K0])
		val[!K0] <- m[!K0] - lambda[!K0]*( 1 - (1 - p[!K0])^(-k[!K0]))/k[!K0]
	}
	val[p <= 0 | p >= 1] <- NA
    val
}

rpareto <- function(n , m = 0, lambda = 1, xi = 0) 
{
   qpareto(runif(n),m,lambda,xi)
}



#######################
### Estimation from data

gpd.lmom <- function(lmom, location = NA, sample = NA)
{ 
        paramest <- c(NA,NA,NA)
        names(paramest) <- c("m","lambda","k")
        if (SHAPE.XI)  names(paramest) <- c("m","lambda","xi")

        if (!is.na(location)) 
        {
            if (length(lmom) > 4)
            {
                sample <- lmom
                lmom <- sample.LMOM(sample)
            }
        k <- lmom[1]/lmom[2] -2
        lambda <- (1+k) * lmom[1]
        m <- location
        }
        else 
        {
           if (length(lmom) > 4)
           {
               sample <- lmom
               lmom <- sample.LMOM(lmom)
           }
           if (is.na(sample[1]))
           {
               stop(paste("Problem in function gpd.lmom: either location parameter",
                    " or the sample observations should be specified"))
           }
           xx <- min(sample)
           n <- length(sample)
           k <- (n*(lmom[1] - xx) - 2*(n - 1) * lmom[2])/
                   ((n - 1)*lmom[2] - (lmom[1] - xx))
           lambda <- (1 + k) *(2 + k) * lmom[2]
           m <- xx - lambda/(n + k)
        }
        paramest[3] <- k;
        paramest[2] <- lambda;
        paramest[1] <- m;     
       
        if (SHAPE.XI) 
        { 
            paramest[3] <- -paramest[3]
            names(paramest) <-  c("m","lambda","xi")
        }
        
        val <- list(param.est = paramest)       
        val
}

gpd.ml <- function(sample, location = NA, init.est = NA, epsilon = 1e-6) 
{
    if(is.data.frame(sample))
        sample <- as.matrix(sample)

    n <- length(sample)
    lmomest <- init.est   
    if (!is.na(location)) 
    {   
        tmpX <- sample-location
        tmpX <- tmpX[tmpX>0]
        tmpN <- length(tmpX)
        assign("tempX",tmpX, pos = 1)
        assign("tempN",tmpN, pos =1)
        if (is.na(lmomest[1]))
        {
            lmomest <- gpd.lmom(lmom=sample.LMOM(sample), location = location)$param.est
            if (SHAPE.XI)  lmomest[3] <- -lmomest[3]
        }
        x0 <- c(lmomest[2],lmomest[3])
        negative.log.likelihood <- function(theta)
        {
            # I use ll <- -10^10 for the function "optim" does not
            # seem to like NA's or Inf
            k <- theta[2]
            lambda <- theta[1]
            xsc <- 1 - (k*(tempX))/lambda
			# ll <- NA
			ll <- -10^10
            if (sum(xsc < 0) > 0 | lambda < 0)
              # ll <- NA
              ll <- -10^10
            else
                ll <- -tempN*log(lambda) + (1/k - 1)*sum(log(xsc))
            -ll
        }

        fit <- nlm(negative.log.likelihood,x0, iterlim = 200, steptol = epsilon)
        converged <- FALSE
        if(fit$code<4) converged <- TRUE              
        if (!converged)
        {
            #LMOM estimate might be bad... Try moment etimate as the intial starting point...
            tempMean <- mean(tempX)
            CV <- (tempMean * tempMean)/var(tempX)
            x0[1] <- 0.5 * tempMean * (CV+ 1)
            x0[2] <- 0.5 *(CV - 1)
            fit <- nlm(negative.log.likelihood,x0, iterlim = 200)
            if (fit$code>3)
                warning("Maximum Likelihood Method for the GPD did not converge")              
        }
        paramest <- c(location,fit$estimate[1], fit$estimate[2])
        names(paramest) <- c("m","lambda","k")
        if (SHAPE.XI)  names(paramest) <- c("m","lambda","xi")
        if (SHAPE.XI)  paramest[3] <- -paramest[3]    
    }        
    else
    {           
        assign("tempX",sample,pos= 1)
        assign("tempN",length(sample), pos =1)
        if (is.na(lmomest[1]))
        {
              lmomest <- gpd.lmom(lmom=sample.LMOM(sample), sample = sample)$param.est
              if (SHAPE.XI)  lmomest[3] <- -lmomest[3]
        }
        negative.log.likelihood <- function(theta)
        {
            # I use ll <- -10^10 for the function "optim" does not
            # seem to like NA's or Inf
              k <- theta[3]
              m <- theta[1]
              lambda <- theta[2]
              xsc <- 1 - k*(tempX-m)/lambda
              # ll <- NA
              ll <- -10^10
              if (sum(xsc < 0) > 0)
              #      ll <- NA
                    ll <- -10^10
              else
                    ll <- -tempN*log(lambda) + (1/k - 1)*sum(log(xsc))
              -ll
        }

        fit <- optim(lmomest, negative.log.likelihood, method="L-BFGS-B", lower=c(-Inf,0,-Inf), upper=c(min(tempX),Inf,Inf))
        converged <- FALSE
        if (fit$convergence==0) converged <- TRUE
        if (!converged)
               warning("Maximum Likelihood Method for the GPD did not converge")

        paramest <- fit$par
        names(paramest) <- c("m","lambda","k")
        if (SHAPE.XI)  names(paramest) <- c("m","lambda","xi")
        if (SHAPE.XI)  paramest[3] <- -paramest[3]
    }
    val <- list(n=n, data=sample, param.est = paramest, converged = converged)
    return(val)
}


sample.LMOM <- function(x) 
{
          x <- sort(x)
          N <- length(x)
          i <- c(1:N)

    
          fn1 <- N - i
          fn2 <- (N - i - 1)*fn1
          fn3 <- (N - i - 2)*fn2
          
          a1 <- sum(fn1/(N-1) * x)/N
          a2 <- sum(fn2/(N-1)/(N-2) * x)/N
          a3 <- sum(fn3/(N-1)/(N-2)/(N-3) * x)/N

          l1 <- mean(x)
          l2 <- l1 - 2*a1
          tau2 <- (l1 - 6.0*a1 + 6.0*a2)/l2
          tau3 <- (l1 - 12.0*a1 + 30.0*a2 - 20.0*a3)/l2   

          val <- c(l1,l2,tau2,tau3)
          names(val) <- c("ell_1", "ell_2", "tau_3", "tau_4")        
          
          val
}

fit.gpd <- function(data, tail = "two", upper = NA, lower = NA, upper.method = "ml", lower.method = "ml", plot = TRUE, ...)
{
	# returns an object of class "gpd". Used to be gpd.tail
    	# Called "gpd.tail" to avoid confusion with McNeil's function "gpd"
    	# Was "pot.1tail.est' and "pot.2tails.est" in EVANESCE

	if (!( tail=="upper" | tail=="lower" | tail=="two"))
		stop("The parameter should be one of the character strings 'upper', 'lower' or 'two'") 
    	if(is.data.frame(data)) 
        data <- as.matrix(data)
    	n <- length(data)
    	if(is.na(upper)) 
    	{
        		sorted.data <- sort(data)
        		if(n <= 150/0.15)
        		{
           		 	uu1 <- sorted.data[n - trunc(n * 0.15)]
           		 	uu2 <- sorted.data[n - trunc(n * 0.15) - 1]
            			upper <- (uu1 + uu2)/2
        		}
        		else upper <- sorted.data[n - 150]
    	}
    	if(is.na(lower)) 
    	{
        		sorted.data <- sort(data)
       		  	if(n <= 150/0.15)
        		{
           			 uu1 <- sorted.data[trunc(n * 0.15)]
            			uu2 <- sorted.data[trunc(n * 0.15) + 1]
            			lower <- (uu1 + uu2)/2
        		}
        		else lower <- sorted.data[150]
    	}

    # Analysis of the upper tail!
    	upper.exceed <- data[data > upper]
    	upper.excess <- upper.exceed - upper
    	if(casefold(upper.method, upper = FALSE) == "ml") 
    	{
        		gpd.est.res <- gpd.ml(sample = upper.excess, location = 0)
        		gpd.est <- gpd.est.res$param.est
        		upper.converged <- gpd.est.res$converged
        		if(!upper.converged)
           		 warning(" MLE method for GPD did not converge for the upper tail. \n",
           			 "You can try to set the option upper.method = \"lmom\" for upper tail")
    	}
    	else if(casefold(upper.method, upper = FALSE) == "lmom") 
    		{
       		 	lmom <- sample.LMOM(upper.excess)
       		 	gpd.est <- gpd.lmom(lmom, sample = upper.excess, location = 0)$param.est
        			upper.converged <- NA
    		}
    		else
        			stop(paste("Unknown method for the parameter estimation: ", method))

    	upper.par.ests <- c(gpd.est["lambda"], gpd.est["xi"])
    	n.upper.exceed <- length(upper.excess)
    	p.less.upper.thresh <- 1 - n.upper.exceed/n
    
      # Analysis of the lower tail!
        	lower.exceed <- data[data < lower]
       	 	lower.excess <- lower - lower.exceed
        	if(casefold(lower.method, upper = FALSE) == "ml")
       	 {
            	gpd.est.res <- gpd.ml(sample = lower.excess,location = 0)
            	gpd.est <- gpd.est.res$param.est
            	lower.converged <- gpd.est.res$converged
            	if(!lower.converged)
            	warning(" MLE method for GPD did not converge for the lower tail. \n",
                		"You can try to set the option lower.method = \"lmom\" for the lower tail")
        	}
        	else if(casefold(lower.method, upper = FALSE) == "lmom")
        		{
            			gpd.est <- gpd.lmom(sample = lower.excess, location = 0)$param.est
            			lower.converged <- NA
        		}
        		else
            		stop(paste("Unknown method for the parameter estimation: ", method))

        lower.par.ests <- c(gpd.est["lambda"], gpd.est["xi"])
        n.lower.exceed <- length(lower.excess)
        p.larger.lower.thresh <- 1 - n.lower.exceed/n

	if (tail == "two")
	{
		out <- new("gpd.two.tails", 		
			tail = tail,
            	n = length(data),   
            	data = sort(data),
            	upper.exceedances = upper.exceed, 
            	lower.exceedances = lower.exceed, 
            	upper.threshold = upper,
            	lower.threshold = lower, 
            	p.less.upper.threshold= p.less.upper.thresh, 
            	p.larger.lower.threshold =  p.larger.lower.thresh, 
            	n.upper.exceedances = n.upper.exceed, 
            	n.lower.exceedances = n.lower.exceed, 
            	upper.method = upper.method, 
            	lower.method = lower.method,
            	upper.par.ests = upper.par.ests, 
            	lower.par.ests = lower.par.ests, 
            	upper.converged = upper.converged, 
            	lower.converged = lower.converged) 
		if(plot) 
    		{
        		par.orig <- par()
       		 	par(mfrow = c(2, 1))
       		 	qq <- qpareto(ppoints(upper.excess), xi = upper.par.ests["xi"])
        		plot(qq, sort(upper.excess), xlab = paste("GPD Quantiles, for xi = ", upper.par.ests["xi"]),ylab="Excess over threshold", ...)
       		 	title("Upper Tail")
       		 	qq <- qpareto(ppoints(lower.excess), xi = lower.par.ests["xi"])
        		plot(qq, sort(lower.excess), xlab = paste("GPD Quantiles, for xi = ", lower.par.ests["xi"]),ylab="Excess over threshold", ...)
       		 	title("Lower Tail")
       		 	par(mfrow = c(1, 1))
    		}
	}
	if (tail == "upper")
	{
		out <- new("gpd.upper.tail", 		
			tail = tail,
            	n = length(data),   
            	data = sort(data),
            	upper.exceedances = upper.exceed, 
            	upper.threshold = upper,
            	p.less.upper.threshold= p.less.upper.thresh, 
            	n.upper.exceedances = n.upper.exceed, 
            	upper.method = upper.method, 
            	upper.par.ests = upper.par.ests, 
            	upper.converged = upper.converged)
		if(plot) 
    		{
        			par.orig <- par()
       		 	qq <- qpareto(ppoints(upper.excess), xi = upper.par.ests["xi"])
        		plot(qq, sort(upper.excess), xlab = paste("GPD Quantiles, for xi = ", upper.par.ests["xi"]),ylab="Excess over threshold", ...)
       		 	title("Upper Tail")
    		}
	}
	if (tail == "lower")
	{
		out <- new("gpd.lower.tail", 		
			tail = tail,
            	n = length(data),   
            	data = sort(data),
            	lower.exceedances = lower.exceed, 
            	lower.threshold = lower, 
            	p.larger.lower.threshold =  p.larger.lower.thresh, 
            	n.lower.exceedances = n.lower.exceed, 
            	lower.method = lower.method,
            	lower.par.ests = lower.par.ests, 
            	lower.converged = lower.converged) 
		if(plot) 
    		{
        		par.orig <- par()
       		 	qq <- qpareto(ppoints(lower.excess), xi = lower.par.ests["xi"])
        		plot(qq, sort(lower.excess), xlab = paste("GPD Quantiles, for xi = ", lower.par.ests["xi"]),ylab="Excess over threshold", ...)
       		 	title("Lower Tail")
    		}
	}
   	 return(out)
}

############################
###### Methods for the gpd class & Standard Computations for the GPD distributions

gpd.1p <- function(obj, x, linear = TRUE)
{
    if(is(x, "timeSeries")) {
        x <- seriesData(x)
    }
    if(is.data.frame(x)) {
        x <- as.matrix(x)
    }
    x.orig <- x
    x <- sort(x)
    if(class(obj) != "gpd.upper.tail") {
        stop("Wrong object. Parameter obj has to be of class gpd.upper.tail")
    }
    n <- length(x)
    k <- obj@upper.par.ests["xi"]
    if(!SHAPE.XI) {
        k <-  - k
    }
    ndata <- obj@n
    u <- obj@upper.threshold
    pp <- (ppoints(obj@data))
    small <- x <= u
    val <- vector(length = n, mode = "numeric")
    xsm <- as.double(x[small])
    nsm <- as.integer(length(xsm))
    lsmallpts <- as.integer(sum(obj@data <= u) + 1)
    smallpts <- as.double(obj@data[1:lsmallpts])
    oldind <- vector(length = nsm, mode = "numeric")

#    oldind[1:nsm] <- -1
#   indB <- .C("empirfunc",
#   xsm,
#  	smallpts,
# 	nsm,
#	lsmallpts,
#	as.integer(oldind))[[5]]

if (nsm>0)
{
	OB <- outer(rep(1,nsm),smallpts,"*")
	OX <- outer(xsm,rep(1,lsmallpts),"*")
	indB <- pmin(apply(OX>OB,1,sum),nsm-1)  

    indB <- indB + 1
    lvalB <- obj@data[indB]
    indL <- indB - 1
    indL[indL == 0] <- NA
    lvalS <- obj@data[indL]
    lvalS[is.na(lvalS)] <- (x[small])[is.na(lvalS)]
    lvalB[is.na(lvalB)] <- 0
    pvalB <- pp[indB]
    pvalS <- pp[indL]
    #lvalS[is.na(lvalS)] <- 0
    #lvalB[is.na(lvalB)] <- 0
    pvalS[is.na(pvalS)] <- 0
    pvalB[is.na(pvalB)] <- 0
    if(linear) {
        val[small] <- pvalS + ((pvalB - pvalS) * (x[small] - lvalS))/(lvalB -
            lvalS)
    }
    else {
        val[small] <- pvalS
    }
}
    # this is the estimate of F at u:
    pu <- pp[lsmallpts - 1] + ((pp[lsmallpts] - pp[lsmallpts - 1]) * (u - obj@data[lsmallpts - 1]))/(obj@data[lsmallpts] - obj@data[lsmallpts - 1])
    #Nu <- lsmallpts
    valsm <- 1 - (1 - pu) * (1 + (k * (x[!small] - u))/obj@upper.par.ests["lambda"])^
        (-1/k)
    valsm[((k * (x[!small] - u))/obj@upper.par.ests["lambda"]) <= -1] <- 1
    val[!small] <- valsm
    val.orig <- val
    val.orig[sort.list(x.orig)] <- val
    val.orig
}

gpd.1q <- function(obj, p, linear = TRUE)
{
    x.orig <- p
    p <- sort(p)
    if(class(obj) != "gpd.upper.tail") {
        stop("Wrong object: the parameter obj has to be of class gpd.upper.tail")
    }
    n <- length(p)
    val <- vector(length = n, mode = "numeric")
    goodp <- p >= 0 & p <= 1
    val[!goodp] <- NA
    k <- obj@upper.par.ests["xi"]
    if(!SHAPE.XI) {
        k <-  - k
    }
    ndata <- obj@n
    u <- obj@upper.threshold
    pp <- (ppoints(obj@data))
    lsmallpts <- as.integer(sum(obj@data <= u) + 1)
    pu <- pp[lsmallpts - 1] + ((pp[lsmallpts] - pp[lsmallpts - 1]) * (u - obj@data[lsmallpts - 1]))/(obj@data[lsmallpts] - obj@data[lsmallpts - 1])
    small <- (p < pu) & goodp
    psm <- as.double(p[small])
    nsm <- as.integer(length(psm))
    smallpts <- as.double(pp[1:lsmallpts])
    oldind <- vector(length = nsm, mode = "numeric")
    oldind[1:nsm] <- -1
#cat("psm = ",psm,"\n")
#cat("nsm = ",nsm,"\n")

#    lowInd <- .C("empirfunc",
#        psm,
#        smallpts,
#        nsm,
#        lsmallpts,
#        as.integer(oldind))[[5]]
OB <- outer(rep(1,nsm),smallpts,"*")
OX <- outer(psm,rep(1,lsmallpts),"*")
lowInd <- pmin(apply(OX>OB,1,sum),nsm-1) 


    highInd <- lowInd + 1
    lowInd[lowInd <= 0] <- NA
    lowP <- pp[lowInd]
    lowP[is.na(lowP)] <- 0
    highP <- pp[highInd]
    highVal <- obj@data[highInd]
    lowVal <- obj@data[lowInd]
    lowVal[is.na(lowVal)] <- 0
    if(linear) {
        val[small] <- lowVal + ((highVal - lowVal) * (psm - lowP))/(highP - lowP)
    }
    else {
        val[small] <- lowVal
    }
    quant <- u + (obj@upper.par.ests["lambda"] * (((1 - p[!small & goodp])/(1 - pu))^( - k) - 1))/k
    val[(!small & goodp)] <- quant
    val.orig <- val
    val.orig[sort.list(x.orig)] <- val
    val.orig
}

gpd.2p <- function(obj, x, linear = TRUE)
{
    if(is(x, "timeSeries")) {
        x <- seriesData(x)
    }
    if(is.data.frame(x)) {
        x <- as.matrix(x)
    }
    x.orig <- x
    x <- sort(x)
    if(class(obj) != "gpd.two.tails") {
        stop("Wrong object. Parameter obj has to be of class gpd.two.tails")
    }
    N <- length(x)
    val <- vector(length = N)
    n <- obj@n
    meadX <- obj@lower.threshold <= x & x <= obj@upper.threshold
    pp <- ppoints(obj@data)
    xsm <- as.double(x[meadX])
    nsm <- as.integer(length(xsm))
    mdata.first.ind <- sum(obj@data < obj@lower.threshold)
    mdata.last.ind <- sum(obj@data <= obj@upper.threshold) + 1
    meadpts <- as.double(obj@data[mdata.first.ind:mdata.last.ind])
    lmeadpts <- as.integer(length(meadpts))
    oldind <- vector(length = nsm, mode = "numeric")
#cat("xsm = ",xsm,"\n")
#cat("nsm = ",nsm,"\n")
#    oldind[1:nsm] <- -1
#   indB <- .C("empirfunc",
#        xsm,
#        meadpts,
#        nsm,
#        lmeadpts,
#        as.integer(oldind))[[5]] + mdata.first.ind
if(nsm>0)
{
	OB <- outer(rep(1,nsm),meadpts,"*")
	OX <- outer(xsm,rep(1,lmeadpts),"*")
	indB <- pmin(apply(OX>OB,1,sum),nsm-1)  + mdata.first.ind
 
    lvalB <- obj@data[indB]
    indL <- indB - 1
    indL[indL == 0] <- NA
    lvalS <- obj@data[indL]
    lvalS[is.na(lvalS)] <- (x[meadX])[is.na(lvalS)]
    lvalB[is.na(lvalB)] <- 0
    pvalB <- pp[indB]
    pvalS <- pp[indL]
    #lvalS[is.na(lvalS)] <- 0
    #lvalB[is.na(lvalB)] <- 0
    pvalS[is.na(pvalS)] <- 0
    pvalB[is.na(pvalB)] <- 0
    if(linear) {
        val[meadX] <- pvalS + ((pvalB - pvalS) * (x[meadX] - lvalS))/(lvalB - lvalS)
    }
    else  val[meadX] <- pvalS
}  
    # this is the estimate of F at upper threshold:
    p.upper <- pp[mdata.last.ind - 1] + ((pp[mdata.last.ind] - pp[mdata.last.ind -
        1]) * (obj@upper.threshold - obj@data[mdata.last.ind - 1]))/
        (obj@data[mdata.last.ind] - obj@data[mdata.last.ind -1])
    p.lower <- pp[mdata.first.ind] + ((pp[mdata.first.ind + 1] - pp[mdata.first.ind
        ]) * (obj@lower.threshold - obj@data[mdata.first.ind]))/
        (obj@data[mdata.first.ind + 1] - obj@data[mdata.first.ind])
    uper.tail.x <- x > obj@upper.threshold
    lower.tail.x <- x < obj@lower.threshold
    k <- obj@upper.par.ests["xi"]
    if(!SHAPE.XI) k <-  - k
    a <- obj@upper.par.ests["lambda"]
    b <- obj@upper.threshold
    val[uper.tail.x] <- 1 - (1 - p.upper) * (1 + (k * (x[uper.tail.x] - b))/a)^(-1/k)
    if(k < 0 & sum(x > b - a/k) > 0) {
        val[x > b - a/k] <- 1.
    }
    k <- obj@lower.par.ests["xi"]
    if(!SHAPE.XI)  k <-  - k
 
    a <- obj@lower.par.ests["lambda"]
    b <- obj@lower.threshold
    val[lower.tail.x] <- p.lower * (1 - (k * (x[lower.tail.x] - b))/a)^(-1/k)
    if(k < 0 & sum(x < b + a/k) > 0) {
        val[x < b + a/k] <- 0
    }
    val.orig <- val
    val.orig[sort.list(x.orig)] <- val
    val.orig
}

gpd.2q <- function(obj, p, linear = TRUE)
{
    x.orig <- p
    p <- sort(p)
    if(class(obj) != "gpd.two.tails") {
        stop("Wrong object. Parameter obj has to be of class gpd.two.tails")
    }
#    if(is.na(obj@lower.thres) || is.na(obj@lower.par.ests["xi"])) {
#       stop("Object must have 2 tails estimated")
#   }
    N <- length(p)
    val <- vector(length = N, mode = "numeric")
    goodp <- p >= 0 & p <= 1
    val[!goodp] <- NA
    n <- obj@n
    mdata.first.ind <- sum(obj@data < obj@lower.threshold)
    mdata.last.ind <- sum(obj@data <= obj@upper.threshold) + 1
    pp <- ppoints(obj@data)
    p.upper <- pp[mdata.last.ind - 1] + ((pp[mdata.last.ind] - pp[mdata.last.ind -1]) * (obj@upper.threshold - obj@data[mdata.last.ind - 1]))/(obj@data[mdata.last.ind] - obj@data[mdata.last.ind -1])
    p.lower <- pp[mdata.first.ind] + ((pp[mdata.first.ind + 1] - pp[mdata.first.ind]) * (obj@lower.threshold - obj@data[mdata.first.ind]))/
 (obj@data[mdata.first.ind + 1] - obj@data[mdata.first.ind])
    meadX <- p.lower <= p & p <= p.upper
    xsm <- as.double(p[meadX])
    nsm <- as.integer(length(xsm))
    meadpts <- as.double(pp[mdata.first.ind:mdata.last.ind])
    lmeadpts <- as.integer(length(meadpts))
    oldind <- vector(length = nsm, mode = "numeric")
#    oldind[1:nsm] <- -1
    if (nsm==0) indB <- -1 + mdata.first.ind
#    else indB <- .C("empirfunc",
#        			xsm,
#        			meadpts,
#        			nsm,
#        			lmeadpts,
#        			as.integer(oldind))[[5]] + mdata.first.ind
   else {
OB <- outer(rep(1,nsm),meadpts,"*")
OX <- outer(xsm,rep(1,lmeadpts),"*")
indB <- pmin(apply(OX>OB,1,sum),nsm-1) + mdata.first.ind
}
    lvalB <- pp[indB]
    indL <- indB - 1
    indL[indL == 0] <- NA
    lvalS <- pp[indL]
    lvalS[is.na(lvalS)] <- 0
    lvalB[is.na(lvalB)] <- 0
    pvalB <- obj@data[indB]
    pvalS <- obj@data[indL]
    pvalS[is.na(pvalS)] <- 0
    pvalB[is.na(pvalB)] <- 0
    if(linear) {
        val[meadX] <- pvalS + ((pvalB - pvalS) * (p[meadX] - lvalS))/(lvalB - lvalS)
    }
    else  val[meadX] <- pvalS
   
    # this is the estimate of F at upper threshold:
    uper.tail.x <- p > p.upper & goodp
    lower.tail.x <- p < p.lower & goodp
    k <- obj@upper.par.ests["xi"]
    if(!SHAPE.XI) {k <-  - k}
    a <- obj@upper.par.ests["lambda"]
    b <- obj@upper.threshold
    val[uper.tail.x] <- b + (a * (((1 - p[uper.tail.x])/(1 - p.upper))^( - k) -1))/k
    k <- obj@lower.par.ests["xi"]
    if(!SHAPE.XI) k <-  - k
    a <- obj@lower.par.ests["lambda"]
    b <- obj@lower.threshold
    val[lower.tail.x] <- b - (a * (((p[lower.tail.x])/(p.lower))^( - k) - 1))/k
    val.orig <- val
    val.orig[sort.list(x.orig)] <- val
    val.orig
}


setGeneric("pgpd", function(dist, x) standardGeneric("pgpd"))
#setMethod("pgpd","gpd",
#	function(dist, x) 
#	{
#		tail <- dist@tail
#		if ( tail == "two")
#			p <- gpd.2p(dist,x)
#		if ( tail == "upper")
#			p <- gpd.1p(dist,x)
#		if ( tail == "lower")
#		{
#			tmpdist <- new("gpd",
#					tail = "upper",
#          			n = dist@n,   
#            			data = - dist@data,
#            			upper.exceedances = - dist@lower.exceedances, 
#            			lower.exceedances = NA, 
#            			upper.threshold = - dist@lower.threshold,
#            			lower.threshold = NA, 
#            			p.less.upper.threshold= dist@p.larger.lower.threshold, 
#            			p.larger.lower.threshold =  NA, 
#            			n.upper.exceedances = dist@n.lower.exceedances, 
#            			n.lower.exceedances = NA, 
#            			upper.method = dist@lower.method, 
#            			lower.method = NA,
#            			upper.par.ests = dist@lower.par.ests, 
#            			lower.par.ests = NA, 
#            			upper.converged = dist@lower.converged, 
#            			lower.converged = NA) 
#			p <- 1 - gpd.1p(tmpdist, -x)
#		}			
#		return(p)
#	})
setMethod("pgpd","gpd.two.tails",
	function(dist, x) {gpd.2p(dist,x)})
setMethod("pgpd","gpd.upper.tail",
	function(dist, x) {gpd.1p(dist,x)})
setMethod("pgpd","gpd.lower.tail",
	function(dist, x) 
	{
			tmpdist <- new("gpd.upper.tail",
					tail = "upper",
            			n = dist@n,   
            			data = - dist@data,
            			upper.exceedances = - dist@lower.exceedances, 
            			upper.threshold = - dist@lower.threshold,
            			p.less.upper.threshold= dist@p.larger.lower.threshold, 
            			n.upper.exceedances = dist@n.lower.exceedances, 
            			upper.method = dist@lower.method, 
            			upper.par.ests = dist@lower.par.ests, 
            			upper.converged = dist@lower.converged) 
			p <- 1 - gpd.1p(tmpdist, -x)
		return(p)
	})


setGeneric("qgpd", function(dist, x) standardGeneric("qgpd"))
setMethod("qgpd","gpd.two.tails",
	function(dist, x) {gpd.2q(dist,x)})
setMethod("qgpd","gpd.upper.tail",
	function(dist, x) {gpd.1q(dist,x)})
setMethod("qgpd","gpd.lower.tail",
	function(dist, x) 
	{
			tmpdist <- new("gpd.upper.tail",
					tail = "upper",
            			n = dist@n,   
            			data = - dist@data,
            			upper.exceedances = - dist@lower.exceedances, 
            			upper.threshold = - dist@lower.threshold,
            			p.less.upper.threshold= dist@p.larger.lower.threshold, 
            			n.upper.exceedances = dist@n.lower.exceedances, 
            			upper.method = dist@lower.method, 
            			upper.par.ests = dist@lower.par.ests, 
            			upper.converged = dist@lower.converged)
			q <- gpd.1q(tmpdist, 1-x)			
		return(q)
	})

setGeneric("dgpd", function(dist, x) standardGeneric("dgpd"))
setMethod("dgpd","gpd.two.tails",
	function(dist, x) 
	{
		deltax <- max(x) - min(x)
		if (deltax == 0)
		{
			maxx <- max(x+1)
			minx <- min(x-1)
		}
		else
		{
			maxx <- max(x) + deltax/2 
			minx <- min(x) - deltax/2
		}
		n <- 1000
		divx <- seq(from=minx, to=maxx,length=n+1)
		cdfx <- pgpd(dist,divx)
		cdf.spl <- smooth.spline(divx, cdfx)
		densx <- predict(cdf.spl, x, deriv = 1)$y
		densx <- kmooth(x,pmax(densx,0))$y
		return(densx)
	})
setMethod("dgpd","gpd.upper.tail",
	function(dist, x) 
	{
		deltax <- max(x) - min(x)
		if (deltax == 0)
		{
			maxx <- max(x+1)
			minx <- min(x-1)
		}
		else
		{
			maxx <- max(x) + deltax/2 
			minx <- min(x) - deltax/2
		}
		n <- 1000
		divx <- seq(from=minx, to=maxx,length=n+1)
		cdfx <- pgpd(dist,divx)
		cdf.spl <- smooth.spline(divx, cdfx)
		densx <- predict(cdf.spl, x, deriv = 1)$y
		densx <- kmooth(x,pmax(densx,0))$y
		return(densx)
	})
setMethod("dgpd","gpd.lower.tail",
	function(dist, x) 
	{
		deltax <- max(x) - min(x)
		if (deltax == 0)
		{
			maxx <- max(x+1)
			minx <- min(x-1)
		}
		else
		{
			maxx <- max(x) + deltax/2 
			minx <- min(x) - deltax/2
		}
		n <- 1000
		divx <- seq(from=minx, to=maxx,length=n+1)
		cdfx <- pgpd(dist,divx)
		cdf.spl <- smooth.spline(divx, cdfx)
		densx <- predict(cdf.spl, x, deriv = 1)$y
		densx <- kmooth(x,pmax(densx,0))$y
		return(densx)
	})


setGeneric("rgpd", function(dist, n = 50) standardGeneric("rgpd"))
setMethod("rgpd","gpd",
	function(dist, n = 50) 
	{
		u <- runif(n)
		x <- qgpd(dist,u)
		return(x)
	})
setMethod("rgpd","gpd.two.tails",
	function(dist, n = 50) 
	{
		u <- runif(n)
		x <- qgpd(dist,u)
		return(x)
	})
setMethod("rgpd","gpd.upper.tail",
	function(dist, n = 50) 
	{
		u <- runif(n)
		x <- qgpd(dist,u)
		return(x)
	})
setMethod("rgpd","gpd.lower.tail",
	function(dist, n = 50) 
	{
		u <- runif(n)
		x <- qgpd(dist,u)
		return(x)
	})

############################
###### Plotting functions

tailplot.one.tail <- function(gpd.obj,  optlog = NA, extend = 1.5, labels = TRUE, ...)
{
	tail <- gpd.obj@tail
    	if (tail == "upper")
    	{
        	data <- as.numeric(gpd.obj@upper.exceedances)
        	threshold <- gpd.obj@upper.threshold
        	xi <- gpd.obj@upper.par.ests["xi"]
        	lambda <- gpd.obj@upper.par.ests["lambda"]
    	}
    	if (tail == "lower")
    	{
        	data <- -as.numeric(gpd.obj@lower.exceedances)
        	threshold <- -gpd.obj@lower.threshold
        	xi <- gpd.obj@lower.par.ests["xi"]
        	lambda <- gpd.obj@lower.par.ests["lambda"]
    	}
    	plotmin <- threshold
    	if(extend <= 1)
        	stop("extend must be > 1")
    	plotmax <- max(data) * extend
    	x <- qpareto(seq(from = 0.001, to = 0.999, length = 999),threshold, lambda,xi)
    	x <- pmin(x, plotmax)
    	x <- pmax(x, plotmin)
    	ypoints <- ppoints(sort(data))
    	y <- ppareto(x,  m=threshold, lambda=lambda,xi=xi)
    	type <- "tail"
    	if(!is.na(optlog))
        	alog <- optlog
    	else alog <- "xy"
    	if (tail == "upper") prob <- gpd.obj@p.less.upper.threshold
    	if (tail == "lower") prob <- gpd.obj@p.larger.lower.threshold

    	ypoints <- (1 - prob) * (1 - ypoints)
    	y <- (1 - prob) * (1 - y)
    	shape <- xi
    	scale <- lambda * (1 - prob)^xi
    	location <- threshold - (scale * ((1 - prob)^( - xi) -1))/xi
    	plot(sort(data), ypoints, xlim = range(plotmin, plotmax), 
            ylim = range(ypoints, y, na.rm = TRUE), xlab = "", ylab = "", log = alog, axes = TRUE, ...)
    	lines(x[y >= 0], y[y >= 0])
    	if(labels) 
    	{
        	PlotType <- switch(alog,
             		x = "log scale for x only",
             		xy = "log - log scale",
             		yx = "log - log scale",
             		"natural scale"
             	)
        	xxlab <- switch(tail,
            		upper = "x",
            		lower = "-x"
        	)
        	yylab <- switch(tail,
            		upper = "1-F(x)",
            		lower = "F(x)"
        	)
        	title(main = paste("Plot of ",tail," tail in ", PlotType,sep=""),xlab = xxlab, ylab = yylab)
    	}
    	lastcurve <- list(lastfit = gpd.obj, type = type, dist = "gpd", plotmin = plotmin, plotmax = plotmax,
        	alog = alog, location = as.numeric(location), shape = as.numeric(shape), scale = as.numeric(scale))
    	assign("lastcurve", lastcurve, pos = 1)
    	invisible()
} 

setGeneric("tailplot", function(gpd.obj, optlog = NA, extend = 1.5, labels=TRUE, ...) standardGeneric("tailplot"))
setMethod("tailplot","gpd.two.tails",
	function(gpd.obj, optlog = NA, extend = 1.5, labels = TRUE, ...) 
	{
      		par(mfrow=c(2,1))
		tmpobj <- gpd.obj
		tmpobj@tail <- "upper"
		tailplot.one.tail(tmpobj)
		tmpobj@tail <- "lower"
		tailplot.one.tail(tmpobj)  
		par(mfrow=c(1,1))
    	})

setMethod("tailplot","gpd.upper.tail",
	function(gpd.obj, optlog = NA, extend = 1.5, labels = TRUE, ...) tailplot.one.tail(gpd.obj)) 

setMethod("tailplot","gpd.lower.tail",
	function(gpd.obj, optlog = NA, extend = 1.5, labels = TRUE, ...) tailplot.one.tail(gpd.obj)) 


upper.shape.plot <- function (data, method = "ml", from = 0.5, to = 0.98, nint = 30) 
{
    	if (is(data, "series")) {
       		data <- seriesData(data)
    	}
	if (is.data.frame(data)) {
       		data <- as.matrix(data)
   	 }
    	data <- sort(unclass(data))
    	assign("tempData", data, pos = 1)
    	estFun <- get(paste("gpd", method, sep = "."))
    	assign("tempEstFun", estFun, pos = 1)
    	n <- length(data)
    	l1 <- data[trunc(from * n)]
    	l2 <- data[trunc(to * n)]
    	x <- pretty(c(l1, l2), n = nint)

    	one.y <- function(u) {
        	xx <- tempData[tempData > u]
        	excess <- xx - u
        	gpd.est <- tempEstFun(sample = excess, location = 0)$param.est
		c(gpd.est[3], length(xx)/length(tempData))
   	 }
    	iii <- apply(as.matrix(x), 1, one.y)
    	yy <- iii[1, ]
    	ylim <- range(yy)
    	ylim[1] <- ylim[1] - 0.5
    	ylim[2] <- ylim[2] + 0.5
   	 t1 <- "Estimate of k"
   	 if (SHAPE.XI) 
       		 t1 <- "Estimate of xi"
        plot(x, yy, type = "l", xlab = "Threshold", ylab = t1, ylim = ylim)
        nl <- length(pretty(x))
        xx <- pretty(x)
#cat("x",x,"\n")
#cat("xx",xx,"\n")
#        indB <- .C("empirfunc", as.double(xx), as.double(x), 
#            as.integer(length(xx)), as.integer(length(x)), as.integer(1:length(xx)))[[5]] + 1
OB <- outer(rep(1,as.integer(length(xx))),as.double(x),"*")
OX <- outer(as.double(xx),rep(1,length(x)),"*")
indB <- pmin(apply(OX>OB,1,sum),length(x)-1) +1

#cat("indB",indB,"\n")
#cat("x[indB]",x[indB],"\n")
#cat("dim(iii)",dim(iii),"\n")
#cat("iii[2, indB] * 100",iii[2, indB] * 100,"\n")
  
        axis(3, at = x[indB], lab = paste(format(round(iii[2, indB] * 1000)/10)))
}

shape.plot <- function (data, tail = "upper", method = "ml", from = 0.5, to = 0.98, nint = 30) 
{
	if (tail == "upper")
	{
		upper.shape.plot(data)
		title("Percent Data Points above Threshold\n")
	}
	if (tail == "lower") {
        	data <- -data
		upper.shape.plot(data)
        	title("Percent Data Points below -Threshold\n")
    	}
	if (tail == "two")
	{
		par(mfrow=c(2,1))
		upper.shape.plot(data)
		title("Percent Data Points above Threshold\n")
        	data <- -data
		upper.shape.plot(data)
        	title("Percent Data Points below -Threshold\n")
		par(mfrow=c(1,1))
    	}
}

##############################################
########### Unclassified Misc


block.max <- function (data, overlap = 0, nb.blocks = NA, block.size = 100, 
    method = "ml") 
{
    if ((!is.vector(data)) || (!is.numeric(data))) 
        stop("data should be a numeric vector")
    if (overlap > 50) 
        stop("overlap should be an integer smaller than 50!")
    if (block.size < 100) 
        stop("block.size should be an integer greater than 100!")
    n <- length(data)
    overlap <- as.integer(overlap)
    block <- block.size - overlap
    NB <- as.integer(n/block)
    if (is.na(nb.blocks)) 
        nb.blocks <- NB
    nb.blocks <- min(nb.blocks, NB)
    block.index <- 0:(nb.blocks - 1)
    B <- rev(n - block.index * block)
    A <- B + 1 - block.size
    GOOD <- A > 0
    nb.blocks <- sum(GOOD)
    A <- A[GOOD]
    B <- B[GOOD]
    M <- rep(0, nb.blocks)
    for (I in 1:nb.blocks) M[I] <- max(data[A[I]:B[I]])
    if (method == "lmom") 
        par <- gev.lmom(M)$param.est
    else if (method == "ml") 
        par <- gev.ml(M)$param.est
    else stop("method should be either lmom or ml!")
    out <- list(n = length(data), data = sort(data), method = method, 
        nb.blocks = nb.blocks, overlap = overlap, interval.a = A, 
        interval.b = B, param.est = par)
    oldClass(out) <- "block.max"
    out
}

ploting.positions <- function(sample, gamma = -0.35, delta = 0 ) 
{
       x <- sort(sample)
       n <- length(x)

       pp.pos <- function(i,n, gamma, delta) {
           (i + gamma)/(n + delta)
        }
       
       pp <- pp.pos(c(1:n),n,gamma,delta)
       ell1 <- sum(pp*x)/n
       ell2 <- sum((2*pp-1)*x)/n
       ell3 <- sum((6*pp*pp - 6*pp+1)*x)/n
       ell4 <- sum((20*pp^3-30*pp*pp +12*pp-1)*x)/n
      
       tau3 <- ell3/ell2
       tau4 <- ell4/ell2

       val <- c(ell1,ell2,tau3,tau4)
       names(val) <- c("ell_1", "ell_2", "tau_3", "tau_4")       
       val
}


qqexp <- function (x, nq = 50) 
{
    values <- seq(from = 1e-04, to = 0.99, length = nq)
    qvalues <- qexp(values)
    plot(qvalues, quantile(x, probs = values), xlab = "Theoretical Quantiles", 
        ylab = "Sample Quantiles")
    title("Exponential Q-Q Plot")
    abline(lmrob(quantile(x, probs = values) ~ qvalues))
}

qqnorm <- function (x) 
{
    x <- x[!is.na(x)]
    v <- stats::qqnorm(x)
    abline(lmrob(v$y ~ v$x))
}

eda.shape <- function (x,xname=deparse(substitute(x))) 
{
	cat(xname)
	par(mfrow = c(2, 2))
	hist(x,xlab="", main="",col="blue")
	boxplot(x)
	iqd <- summary(x)[5] - summary(x)[2]
	plot(density(x, width = 2 * iqd), xlab = "", ylab = "", 
						type = "l",main="")
	x <- x[!is.na(x)]    
	v <- stats::qqnorm(x,main="",xlab="")
	abline(lmrob(v$y ~ v$x))
	par(mfrow = c(1, 1))
	title(main=paste("EDA of ",xname,sep=""))
}

