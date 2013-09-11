##################################
### Defining Classes & Methods ###
##################################

##########
# Defining class copula

setClass("copula", representation(parameters = "numeric", param.names = "character", param.lowbnd = "numeric", param.upbnd = "numeric", message = "character"))
   
############
## Remove methods if they already exist
  
if(isGeneric("pcopula")) removeGeneric("pcopula")
if(isGeneric("dcopula")) removeGeneric("dcopula")
if(isGeneric("rcopula")) removeGeneric("rcopula")

if(isGeneric("contour.pcopula")) removeGeneric("contour.pcopula")
if(isGeneric("contour.dcopula")) removeGeneric("contour.dcopula")
if(isGeneric("persp.pcopula")) removeGeneric("persp.pcopula")
if(isGeneric("persp.dcopula")) removeGeneric("persp.dcopula")

if(isGeneric("contour.plot")) removeGeneric("contour.plot")

if(isGeneric("CondCinverse")) removeGeneric("CondCinverse")

if(isGeneric("A")) removeGeneric("A")
if(isGeneric("AsecondDer")) removeGeneric("AsecondDer")
if(isGeneric("AfirstDer")) removeGeneric("AfirstDer")

if(isGeneric("Hderiv")) removeGeneric("Hderiv")
if(isGeneric("PhiDer")) removeGeneric("PhiDer")
if(isGeneric("Phi")) removeGeneric("Phi")
if(isGeneric("InvPhi")) removeGeneric("InvPhi")
if(isGeneric("InvPhisecondDer")) removeGeneric("InvPhisecondDer")
if(isGeneric("InvPhifirstDer")) removeGeneric("InvPhifirstDer")
if(isGeneric("DerPhi.minus1")) removeGeneric("DerPhi.minus1")
if(isGeneric("tail.index")) removeGeneric("tail.index")
if(isGeneric("dcdx")) removeGeneric("dcdx")
if(isGeneric("lambda")) removeGeneric("lambda")

if(isGeneric("Kendalls.tau")) removeGeneric("Kendalls.tau")
if(isGeneric("Spearmans.rho")) removeGeneric("Spearmans.rho")

##############
# Defining class empirical copula
    
setClass("empirical.copula", representation(x = "numeric", y = "numeric"))


classNames <- c("normal","frank","kimeldorf.sampson",
                 "gumbel","galambos","husler.reiss","tawn",
                 "bb5","normal.mix",
                 "bb1","bb2","bb3","bb6","bb7","joe",
                 "bb4","empirical")

############################################
### Set generic methods for class copula ###
############################################

setGeneric("pcopula", function(copula, u, v) standardGeneric("pcopula"))
setGeneric("dcopula", function(copula, u, v) standardGeneric("dcopula"))
setGeneric("rcopula", function(copula, n) standardGeneric("rcopula"))

setGeneric("contour.pcopula", function(cop,n=50, ... ) standardGeneric("contour.pcopula"))
setGeneric("contour.dcopula", function(cop,n=50, ... ) standardGeneric("contour.dcopula"))
setGeneric("persp.pcopula", function(cop,n=50, ... ) standardGeneric("persp.pcopula"))
setGeneric("persp.dcopula", function(cop,n=50, ... ) standardGeneric("persp.dcopula"))

# setGeneric("contour.plot", function(object,n = 100, nlevels = 10, add = FALSE,  ... ) standardGeneric("contour.plot"))

setGeneric("CondCinverse", function(copula, q, u ) standardGeneric("CondCinverse"))

setGeneric("A", function(copula, t) standardGeneric("A"))
setGeneric("AsecondDer", function(copula, t) standardGeneric("AsecondDer"))
setGeneric("AfirstDer", function(copula, t) standardGeneric("AfirstDer"))
setGeneric("Hderiv", function(copula, z) standardGeneric("Hderiv"))

setGeneric("Phi", function(copula, t) standardGeneric("Phi"))
setGeneric("PhiDer", function(copula, t) standardGeneric("PhiDer"))
setGeneric("InvPhi", function(copula, s) standardGeneric("InvPhi"))
setGeneric("InvPhisecondDer", function(copula, s) standardGeneric("InvPhisecondDer"))
setGeneric("InvPhifirstDer", function(copula, s) standardGeneric("InvPhifirstDer"))
setGeneric("DerPhi.minus1", function(copula, s) standardGeneric("DerPhi.minus1"))

setGeneric("Kendalls.tau", function(copula, tol = 1e-5) standardGeneric("Kendalls.tau"))
setGeneric("Spearmans.rho", function(copula, tol = 1e-5) standardGeneric("Spearmans.rho"))

setGeneric("lambda", function(copula, t) standardGeneric("lambda"))
setGeneric("tail.index", function(copula, ...) standardGeneric("tail.index"))
setGeneric("dcdx", function(copula, u, v) standardGeneric("dcdx"))

 
#############################################
## copula classes implemented
#############################################

setClass("normal.copula","copula")
setClass("frank.copula","copula")
setClass("kimeldorf.sampson.copula","copula")
setClass("gumbel.copula","copula")
setClass("galambos.copula","copula")
setClass("husler.reiss.copula","copula")
setClass("tawn.copula","copula")
setClass("bb5.copula","copula")

setClass("normal.mix.copula","copula")

setClass("bb1.copula","copula")
setClass("bb2.copula","copula")
setClass("bb3.copula","copula")
setClass("bb4.copula","copula")
setClass("bb6.copula","copula")
setClass("bb7.copula","copula")
setClass("joe.copula","copula")


setMethod("show","copula", function(object) { cat("    ", object@message,".\n    Parameters : \n", sep = "")
     for (i in (1:length(object@parameters))) cat("      ", object@param.names[i]," = ",object@parameters[i],"\n")   
   })



evalFunc <- function(x, fun, param) 
{
   n <- length(param)
   call <- switch(n,
                    expression(fun(x,param[1])),
                expression(fun(x,param[1],param[2])),
                 expression(fun(x,param[1],param[2],param[3])),
                 expression(fun(x,param[1],param[2],param[3],param[4])))
   eval(call)
}



##################################
### Empirical  Copulas Methods ###
##################################
    
setMethod("Kendalls.tau","empirical.copula", 
function(copula, tol = 1e-5) 
{   
   tau <- cor(c(copula@x),c(copula@y),method="k")
   tau          
})  
    
setMethod("Spearmans.rho","empirical.copula", function(copula, tol = 1e-5) 
{
   tau <- cor(c(copula@x),c(copula@y),method="s")
   tau          
}   )
    

setMethod("A","empirical.copula", function(copula, t) {

   u <- copula@x
   v <- copula@y

   p <- function (t) {
       1 - t
   }
    
    n <- length(u)
    N <- length(t)
    
    Z <- sort(log(u)/log(u*v))
    Q <- cumprod((Z/(1-Z))^(1/n))
    val <- vector(length = N, mode = "numeric")
    val[t < 0 | t > 1] <- 1
    
    block1 <- t >=0 & t <= Z[1]
    block2 <- t > Z[1] & t < Z[n]
    block3 <- t >= Z[n] & t <= 1
    
    
    Qn <- Q[n]
    val[block1] <- (1 - t[block1])*Qn^(1 - p(t[block1]))
    val[block3] <- (t[block3])*Qn^( - p(t[block3]))
    
    tt <- t[block2]
    nt <- length(tt)
    tsort <- sort(tt)
    N2 <- length(tsort)
    v2 <- vector(length = N2, mode = "numeric")
    v2 [1:N2] <- -1
#    outC <- .C("empirfunc",as.double(tsort),as.double(Z),as.integer(N2),
#          as.integer(n),as.integer(v2))[[5]]
OB <- outer(rep(1,N2),Z,"*")
OX <- outer(tsort,rep(1,n),"*")
indB <- pmin(apply(OX>OB,1,sum),N2-1)  
 
    indback <- c(1: N2)
    indback[sort.list(tt)] <- c(1: N2)
    i <- outC[indback] 
    iovn <- i/n
    v2 <- (tt^(iovn))*((1-tt)^(1-iovn))*Qn^(1 - p(tt))/(Q[i])

    val[block2] <- v2
    
    #check fo consistency with 0.5 - 1 range. change to that range if need, but 
    #issue warning
    badval <- (val < 0.5 | val > 1)
    if (sum(badval) > 0) {
           warning(paste("Some of the estimated values of A(t) are not within [0.5,1] range",
           "They are chaged to fit this range", sep = "  "));
           val[val < 0.5] <- 0.5;
           val[val > 1] <- 1.0
        }
    
    val
})


################################################################################
# Fitting copula to data                                                       #
# Generic MLE method                                                           #
################################################################################

fit.copula <- function(data,family="normal", plotPicture=TRUE, init.est=NA, epsilon=1e-5, ...)
{
     if (!inherits(data,"empirical.copula")) {data <- empirical.copula(data)}
     families.implemented <- c("normal","frank","kimeldorf.sampson",
                      "gumbel","galambos","husler.reiss",
                     "bb1","tawn","bb2","bb3","bb5","bb6","bb7","normal.mix","bb4", "joe")

     w <- charmatch(family, families.implemented, nomatch = -1)
     if (w <= 0) {
           message <- paste("Family ",family," is not implemented in fit.copula function")
           message <- paste(message," Valid family names are ",families.implemented)
           stop(message)
     }

     if (is.na(init.est)) {
     Kt <- Kendalls.tau(empirical.copula(data))
     if (Kt <= 0) {
           message <- paste("Negative dependence is detected in fit.copula function")
           message <- paste(message," Negative dependence is not implemented yet.")
           stop(message)
     }
    #starting point of parameters!
    Ktrd <- round(Kt, digits = 1)
    if (Ktrd == 0) {
          Ktrd <- 0.1
          message <- paste("Very small dependence is detected in fit.copula function")
          warning(message)
    }
    if (Ktrd == 1) Ktrd <- 0.9
    nline <- Ktrd *10 + 1
    x0 <- switch(w, Kend.table[nline,2],
          Kend.table[nline,3],Kend.table[nline,4],Kend.table[nline,5],Kend.table[nline,6],
          Kend.table[nline,7],  c(Kend.table[nline,9], Kend.table[nline,8]),
          c(Kend.table[nline,11],Kend.table[nline,11], Kend.table[nline,10]),
          c(Kend.table[nline,12],Kend.table[nline,13]),
          c(Kend.table[nline,14],Kend.table[nline,15]),
          c(Kend.table[nline,16],Kend.table[nline,17]),
          c(Kend.table[nline,18],Kend.table[nline,19]),
          c(Kend.table[nline,20],Kend.table[nline,21]),
          c(0.5, Kend.table[nline,2],Kend.table[nline,2]),
          c(Kend.table[nline,22],Kend.table[nline,23]),
          c(2.0) )
    }
    else{ x0 <- init.est }
    copula <- switch(w, normal.copula(x0), frank.copula(x0), kimeldorf.sampson.copula(x0),
            gumbel.copula(x0), galambos.copula(x0), husler.reiss.copula(x0),
            bb1.copula(x0[1],x0[2]), tawn.copula(x0[1],x0[2],x0[3]),
            bb2.copula(x0[1],x0[2]),
            bb3.copula(x0[1],x0[2]),
            bb5.copula(x0[1],x0[2]),
            bb6.copula(x0[1],x0[2]),
            bb7.copula(x0[1],x0[2]),
            normal.mix.copula(x0[1],x0[2],x0[3]),
            bb4.copula(x0[1],x0[2]), joe.copula(x0[1]))
    copTemp <- copula

    neg.log.lik <- function(param, copula, data) {
               copula@parameters <- param
               val <- dcopula(copula, data@x, data@y)
               nll <- sum(log(val))
               -nll
    }
    fit <- nlminb(x0, neg.log.lik, copula=copula, data=data, lower=copula@param.lowbnd, upper=copula@param.upbnd)
    cat("Maximum Likelihood Estimation of copula:\n")
    cat("MLE terminated due to ",fit$message,"\n")
    copula@parameters <- fit$par
    names(copula@parameters) <- names(copTemp@parameters)

    if (plotPicture) {
		n <- 50
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
    		zmat <- pcopula(data,  xmat,  ymat)
    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    		contour(val, col="red",lty=4)
    		title("Copula Contour Plots", xlab="u",ylab="v")
    		zmat <- pcopula(copula,  xmat,  ymat)
    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat))
    		contour(val, col="blue",axes=FALSE,add=TRUE,drawlabels=FALSE)
	title(sub="Empirical (red) and Fitted (blue) Copula Levet Sets")
    }
    copula
}



#######################
### Copula Families ###
#######################

######### Family B1 ######################
# Bivariate Normal                       #
##########################################
 
normal.copula <- function(delta) {
      if (is.na(delta) | delta >= 1 | delta < 0)
          stop("Parameter of the bivariate normal.copula should be in [0,1)")

      val <- new("normal.copula", parameters = delta, param.names = "delta",
               param.lowbnd = 0, param.upbnd = 1, message = "Normal copula family")
     val
}

##################################
### Mixtures of Normal copulas ###
##################################

normal.mix.copula <- function(p, delta1, delta2) {
      if (delta1 > 1 | delta1 < 0|
          p > 1 | p < 0 | delta2 > 1 | delta2 < 0)
          stop("All parameters for normal.mix.copula should be between 0 and 1")

     val <- new("normal.mix.copula", parameters = c(p,delta1, delta2),
            param.names =  c("p","delta1", "delta2"),
             param.lowbnd = c(0,0,0), param.upbnd = c(1,1,1), 
             message = "Mixed Normal copula family")
    val    
}
    


#################################################
#### Family B4 (Kimerdorf and Sampson 91975)) ###
#################################################

kimeldorf.sampson.copula <- function(delta) {
   if (is.na(delta) | delta < 0)
          stop("Parameter delta for kimeldorf.sampson.copula should be >= 0")

   val <- new("kimeldorf.sampson.copula", parameters = delta,  param.names = "delta",     
               param.lowbnd = 0, param.upbnd = Inf, 
               message = "Kimeldorf and Sampson copula family") 
   val
}

setMethod("CondCinverse","kimeldorf.sampson.copula", function(copula, q, u) {
   delta <- as.vector(copula@parameters[1])
   ((q^(-delta/(delta +1)) -1)*u^(-delta)+1)^(-1/delta)
})

setMethod("Phi","kimeldorf.sampson.copula", function(copula, t) {
    delta <- as.vector(copula@parameters[1])
    t^(-delta) -1
})

setMethod("PhiDer","kimeldorf.sampson.copula", function(copula, t) {
     delta <- as.vector(copula@parameters[1])
     (-delta)* t^(-delta - 1)
})

##################################################################
##################   EXTREME Value copulas #######################
#
# NOTE: all extreme value copulas are 
# created in the archimax setting
# general method for generation for the EV copulas:


####################
# Family   B6  #####
####################

gumbel.copula <- function(delta) 
{
   if (is.na(delta) | delta < 1) 
      stop("The parameter delta of the Gumbel family should be >= 1")
   val <- new("gumbel.copula", parameters = delta, param.names = "delta",
              param.lowbnd = 1, param.upbnd = Inf,
              message = "Gumbel copula family; Extreme value copula")
   val
}

setMethod("A","gumbel.copula", function(copula, t) {
  delta <- as.vector(copula@parameters[1])
  (t^delta+(1-t)^delta)^(1/delta)
})

setMethod("Hderiv","gumbel.copula", function(copula, z) {
   d <- as.vector(copula@parameters[1])
  (d*(-((-1 + z)*z))^(-1 + d))/((1 - z)^d + z^d)^2
})

setMethod("AfirstDer","gumbel.copula", function(copula, t) {
    delta <- as.vector(copula@parameters[1])
        ((-(delta*(1 - t)^(-1 + delta)) + delta*t^(-1 + delta))*
    ((1 - t)^delta + t^delta)^(-1 + delta^(-1)))/delta
    })

setMethod("AsecondDer","gumbel.copula",function(copula, t) {
   d <-  as.vector(copula@parameters[1])
   (-1 + d)*(-((-1 + t)*t))^(-2 + d)*((1 - t)^d + t^d)^
      (-2 + d^(-1))
})

######################################
## Family b7                         #
######################################

galambos.copula <- function(delta) {
   if (is.na(delta) | delta < 0) {
      stop("Parameter delta >=  0 for B7 copula family") }

   val <- new("galambos.copula", parameters = delta, param.names = "delta", 
         param.lowbnd = 0, param.upbnd = Inf,
         message = "Galambos copula family; Extreme value copula")
 
   val  
}

setMethod("A","galambos.copula", function(copula, t) { 
  delta <- as.vector(copula@parameters[1])
  1 - (t^(-delta)+(1-t)^(-delta))^(-1/delta)
})

setMethod("AfirstDer","galambos.copula", function(copula,t) {
   delta <- as.vector(copula@parameters[1])
  ((1 - t)^(-1 - delta) - t^(-1 - delta))/
    ((1 - t)^(-delta) + t^(-delta))^((1 + delta)/delta)
})

setMethod("AsecondDer","galambos.copula", function(copula, t) {
   delta <-  as.vector(copula@parameters[1])
((1 + delta)*(-((-1 + t)*t))^(-2 + delta))/
  (((1 - t)^(-delta) + t^(-delta))^delta^(-1)*
    ((1 - t)^delta + t^delta)^2)
})

######################################
## Family b8                         #
######################################

husler.reiss.copula <- function(delta) {

    if (is.na(delta) | delta < 0) {
      stop("Parameter delta >=  0 for Husler and Reiss copula family") }
  
  val <- new("husler.reiss.copula", parameters = delta, param.names = "delta",
             param.lowbnd = 0, param.upbnd = Inf,
           message = "Husler and Reiss copula family; Extreme value copula")
    val 
    
}

setMethod("A","husler.reiss.copula", function(copula, t) {

  d <- as.vector(copula@parameters[1])  
  ft1 <- pnorm(1/d + 1/2*d*log(t/(1 - t)))
  ft2 <- pnorm(1/d + 1/2*d*log((1-t)/t))
  
  t*ft1 + (1 - t) *ft2
})

setMethod("AfirstDer","husler.reiss.copula", function(copula, t) {
    d <- as.vector(copula@parameters[1])
   ft1 <- pnorm(1/d + 1/2*d*log(t/(1 - t)))
   ft2 <- pnorm(1/d + 1/2*d*log((1-t)/t))

   dn1 <- dnorm(1/d + 1/2*d*log(t/(1 - t))) * d/2/(1 - t)
   dn2 <- dnorm(1/d - 1/2*d*log(t/(1 - t)))  *d/2/t
    
    dn1 + ft1 - dn2 - ft2
    
})

setMethod("AsecondDer","husler.reiss.copula",function(copula, t) {
   d <- as.vector(copula@parameters[1])
   x1 <- 1/d + 1/2*d*log(t/(1 - t))
   x2 <- 1/d - 1/2*d*log(t/(1 - t))
   dterm <- d/2/t/(1-t)

  dnorm(x1)*dterm/(1-t)*(1 - x1*d/2) +
    dnorm(x2)*dterm/t*(1 - x2*d/2)
})


######################################
## Family Tawn                       #
######################################

tawn.copula <- function(a,b,r) 
{
   if (a < 0 | r < 1 |b > 1)
      stop("Parameter a and b >=  0 for Tawn family")

    val <- new("tawn.copula", parameters = c(a, b, r), param.names = c("alpha", "lambda", "r"),
              param.lowbnd = c(0,0,1), param.upbnd = c(1, 1, Inf),  
            message = "Tawn copula family; Extreme Value copula")
    val  
}


setMethod("A","tawn.copula", function(copula, t) {
  a <- as.vector(copula@parameters[1])
  b <- as.vector(copula@parameters[2])
  r <- as.vector(copula@parameters[3])
  1 - b + (b - a)*t + ((a* t)^r + (b*(1 - t))^r)^(1/r)
})

setMethod("AfirstDer","tawn.copula", function(copula, t) 
{
  a <- as.vector(copula@parameters[1])
  b <- as.vector(copula@parameters[2])
  r <- as.vector(copula@parameters[3])
-a + b + ((-(b*r*(b*(1 - t))^(-1 + r)) +
       a*r*(a*t)^(-1 + r))*
     ((b*(1 - t))^r + (a*t)^r)^(-1 + r^(-1)))/r
})

setMethod("Hderiv","tawn.copula",function(copula, z) 
{
   a <- as.vector(copula@parameters[1])
  b <- as.vector(copula@parameters[2])
  r <- as.vector(copula@parameters[3])
1 + ((1 - 2*z)*(-a + b +
    ((-(b^r*r*(1 - z)^(-1 + r)) + a^r*r*z^(-1 + r))*
      (b^r*(1 - z)^r + a^r*z^r)^(-1 + r^(-1)))/r))/
  (1 + b*(-1 + z) - a*z + (b^r*(1 - z)^r + a^r*z^r)^
    r^(-1)) +
 ((1 - z)*z*
   (-(-a + b + ((-(b^r*r*(1 - z)^(-1 + r)) +
          a^r*r*z^(-1 + r))*(b^r*(1 - z)^r + a^r*z^r)^
          (-1 + r^(-1)))/r)^2 + a^r*b^r*(-1 + r)*
     (-((-1 + z)*z))^(-2 + r)*(b^r*(1 - z)^r + a^r*z^r)^
      (-2 + r^(-1))*(1 + b*(-1 + z) - a*z +
      (b^r*(1 - z)^r + a^r*z^r)^r^(-1))))/
  (1 + b*(-1 + z) - a*z + (b^r*(1 - z)^r + a^r*z^r)^r^(-1))^2
})

setMethod("AsecondDer","tawn.copula", function(copula, t) 
{ 
  a <- as.vector(copula@parameters[1])
  b <- as.vector(copula@parameters[2])
  r <- as.vector(copula@parameters[3])
  a^r*b^r*(-1 + r)*(-((-1 + t)*t))^(-2 + r)*
     (b^r*(1 - t)^r + a^r*t^r)^(-2 + r^(-1))
})


###########################################################
#### Family B3 (Frank (1979)) #############################

frank.copula <- function(delta) {
    if (is.na(delta) | delta < 0)
          stop("Parameter delta for frank.copula should be >= 0")

   val <- new("frank.copula", parameters = delta,  param.names = "delta",     
               param.lowbnd = 0, param.upbnd = Inf, 
               message = "Frank copula family") 
   val
}

setMethod("CondCinverse","frank.copula", function(copula, q, u) {
   delta <- as.vector(copula@parameters[1])
   eta <-  1 - exp(-delta)
   -1/delta*log(1-eta/((1/q - 1)*exp(-delta*u) + 1))
})

setMethod("Phi","frank.copula", function(copula, t) {
    delta <- as.vector(copula@parameters[1])
    -log((exp(-delta*t) - 1)/(exp(-delta) -1))
})


setMethod("PhiDer","frank.copula", function(copula, t) {
     delta <- as.vector(copula@parameters[1])
     delta/(exp((delta*t))*(-1 + exp(-(delta*t))))
})

############################################
#############BB1 copula ####################

archm.copula <- function(family = "BB1", param = c(0.8,1.43)) {
     family <- casefold(family, upper=FALSE)

     families.implemented <- c("gumbel","frank",
         "kimeldorf.sampson","joe","bb1","bb2","bb3","bb6","bb7")
       
       w <- charmatch(family, families.implemented, nomatch = -1)
       
       if (w <= 0) {
           message <- paste("Family ",family," is not implemented in function")
        message <- paste(message," Valid family names are ",families.implemented)
            print(message)
           stop(message)
        }
        
       copula <- switch(w, gumbel.copula(param), frank.copula(param), 
               joe.copula(param), bb1.copula(param[1],param[2]),
               bb2.copula(param[1],param[2]),
               bb3.copula(param[1],param[2]),
              bb6.copula(param[1],param[2]),
              bb7.copula(param[1],param[2]))
       copula   
    
}

bb1.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 1 | theta < 0)
      stop("Parameter delta >= 1 and theta >  0 for BB1 family")
   val <- new("bb1.copula", parameters = c(theta,delta),
            param.names = c("theta","delta"),
            param.lowbnd = c(0,1), param.upbnd = c(Inf, Inf),
   message = "BB1 copula family; Archimedean Copula")
   val
}

setMethod("Phi","bb1.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   delta <- as.vector(copula@parameters[2])
   (t ^(-theta) - 1)^(delta)
})


setMethod("PhiDer","bb1.copula",function(copula,t) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
     - delta*theta*(t ^(-theta) - 1)^(delta - 1)* t ^(-theta-1)
})

setMethod("InvPhi","bb1.copula",function(copula,s) {
     theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
        ( 1 + s^(1/delta) ) ^(-1/theta)
})

setMethod("InvPhifirstDer","bb1.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
        (-1/theta)*( 1 + s^(1/delta) ) ^(-1/theta -1)*(1/delta)*s^(1/delta-1)
})

setMethod("InvPhisecondDer","bb1.copula", function(copula, s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
      (s^(-2 + delta^(-1))*(1 + s^delta^(-1))^
     (-2 - theta^(-1))*
    ((-1 + delta)*theta +
      s^delta^(-1)*(1 + delta*theta)))/(delta^2*theta^2)
})


###################################################################
###  Family   BB2                                               ###
###################################################################


bb2.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 0 | theta < 0)
      stop("Parameter delta >= 0 and theta >  0 for BB2 family")

    val <- new("bb2.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
            param.lowbnd = c(0,0), param.upbnd = c(Inf, Inf), 
            message = "BB2 copula family; Archimedean Copula") 
      val
}

setMethod("Phi","bb2.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   delta <- as.vector(copula@parameters[2])
   -1 + exp(delta*(-1 + t^(-theta)))
})

setMethod("PhiDer","bb2.copula", function(copula,t) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
    -(delta*exp(delta*(-1 + t^(-theta)))*t^(-1 - theta)*theta)
})

setMethod("InvPhi","bb2.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
     delta <- as.vector(copula@parameters[2])
  (1 + log(1 + s)/delta)^(-theta^(-1))
})

setMethod("InvPhisecondDer","bb2.copula", function(copula, s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
     (1 + theta + delta*theta + theta*log(1 + s))/
  ((1 + s)^2*theta^2*(delta + log(1 + s))^2*
    ((delta + log(1 + s))/delta)^theta^(-1))
})


####################################################
## Family BB3  #####################################
####################################################

bb3.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 0 | theta < 1)
      stop("Parameter delta >= 0 and theta >  1 for BB3 family")

    val <- new("bb3.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
       param.lowbnd = c(1,0), param.upbnd = c(Inf, Inf), 
     message = "BB3 copula family; Archimedean Copula") 
    val
}

setMethod("Phi","bb3.copula",function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   delta <- as.vector(copula@parameters[2])
    -1 + exp(delta*(-log(t))^theta)

 })

setMethod("PhiDer","bb3.copula",function(copula,t) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
    (delta*exp(delta*(-log(t))^theta)*theta*(-log(t))^theta)/
                   (t*log(t))
})

setMethod("InvPhi","bb3.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
     delta <- as.vector(copula@parameters[2])
   exp(-(log(1 + s)/delta)^theta^(-1))
})

setMethod("InvPhisecondDer","bb3.copula",  function(copula, s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
  ((log(1 + s)/delta)^theta^(-1)*
    (-1 + theta + theta*log(1 + s) +
      (log(1 + s)/delta)^theta^(-1)))/
  (exp(log(1 + s)/delta)^theta^(-1)*(1 + s)^2*theta^2*
    log(1 + s)^2)
})


###################################################################
# Archimax Copulas #################################################
###################################################################
    
setMethod("Hderiv","bb4.copula", function(copula, z) {
       ad0 <- A(copula, z)
       ad1 <- AfirstDer(copula, z)
       ad2 <- AsecondDer(copula, z)

       1 + (1 - 2*z)*ad1/ad0 + z*(1 - z) *(ad2/ad0 - ad1^2/ad0^2)
    }   )

#################################################################
###   BB4 copula                                              ###
################################################################# 

bb4.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 0 | theta < 0)
      stop("Parameter delta > 0 and theta >=  0 for BB4 family")

   val <- new("bb4.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
     param.lowbnd = c(0,0), param.upbnd = c(Inf, Inf), 
     message = "BB4 copula family; Archimax copula") 
   
   val
}

setMethod("Phi","bb4.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   t^(-theta) - 1
})

setMethod("PhiDer","bb4.copula",function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   (-theta)*t^(-theta - 1)
})

setMethod("InvPhi","bb4.copula",function(copula,s) {
   theta <- as.vector(copula@parameters[1])
   (s + 1)^(-1/theta)
})

setMethod("InvPhifirstDer","bb4.copula",function(copula,s) {
   theta <- as.vector(copula@parameters[1])
   (-1/theta)*(s + 1)^(-1/theta - 1)
})

setMethod("InvPhisecondDer","bb4.copula",function(copula,s) {
   theta <- as.vector(copula@parameters[1])
   (-1/theta)*(-1/theta - 1)*(s + 1)^(-1/theta - 2)
})

setMethod("DerPhi.minus1","bb4.copula",function(copula,s) {
   theta <- as.vector(copula@parameters[1])
    (-(s/theta))^(-1 - theta)^(-1)
})

setMethod("A","bb4.copula", function(copula, t) {
    delta <- as.vector(copula@parameters[2])
    1 - (t^(-delta)+(1-t)^(-delta))^(-1/delta)
})

setMethod("AfirstDer","bb4.copula", function(copula,t) {
   delta <- as.vector(copula@parameters[2])
   ((1 - t)^(-1 - delta) - t^(-1 - delta))/
    ((1 - t)^(-delta) + t^(-delta))^((1 + delta)/delta)
})

setMethod("AsecondDer","bb4.copula", function(copula, t) {
   delta <-  as.vector(copula@parameters[2])
 ((1 + delta)*(-((-1 + t)*t))^(-2 + delta))/
  (((1 - t)^(-delta) + t^(-delta))^delta^(-1)*
    ((1 - t)^delta + t^delta)^2)
})

###################################################################
### Family BB5  ##################################################


bb5.copula <- function(theta,delta) {
   if (is.na(delta) | delta < 0 | theta < 1) {
      stop("Parameter delta >  0 and theta >= 1 for BB5 copula family") }

   val <- new("bb5.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
            param.lowbnd = c(1,0), param.upbnd = c(Inf,Inf),
         message = "BB5 copula family; Extreme value copula")
    val  
}


setMethod("A","bb5.copula", function(copula, t) { 
  theta <- as.vector(copula@parameters[1])
  delta <- as.vector(copula@parameters[2])

  ((1 - t)^theta + t^theta - 
    ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
     (-delta^(-1)))^theta^(-1)
})

setMethod("AfirstDer","bb5.copula" ,function(copula,t) {
  theta <- as.vector(copula@parameters[1])
  delta <- as.vector(copula@parameters[2])
  ((1 - t)^theta + t^theta -
     ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
      (-delta^(-1)))^(-1 + theta^(-1))*
  (-(1 - t)^(-1 + theta) + t^(-1 + theta) +
    ((1 - t)^(-1 - delta*theta) - t^(-1 - delta*theta))/
     ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
      ((1 + delta)/delta))

})

 
setMethod("AsecondDer","bb5.copula", function(copula, t) {
   theta <- as.vector(copula@parameters[1])
  delta <- as.vector(copula@parameters[2])

((1 - t)^theta + t^theta - 
     ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
      (-delta^(-1)))^(-2 + theta^(-1))*
  ((-(1 - t)^(-1 + theta) + t^(-1 + theta) + 
        ((1 - t)^(-1 - delta*theta) - 
           t^(-1 - delta*theta))/
         ((1 - t)^(-(delta*theta)) + 
            t^(-(delta*theta)))^((1 + delta)/delta))^2*
     (-1 + theta^(-1))*theta + 
    ((1 - t)^theta + t^theta - 
       ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
        (-delta^(-1)))*((1 - t)^(-2 + theta)*
        (-1 + theta) + t^(-2 + theta)*(-1 + theta) - 
       ((1 + delta)*(-(1 - t)^(delta*theta) + 
             (1 - t)^(delta*theta)*t + 
             t^(1 + delta*theta))^2*theta)/
        ((-1 + t)^2*t^2*
          ((1 - t)^(-(delta*theta)) + 
             t^(-(delta*theta)))^delta^(-1)*
          ((1 - t)^(delta*theta) + t^(delta*theta))^2) + 
       ((-(1 - t)^(-2 - delta*theta) - 
            t^(-2 - delta*theta))*(-1 - delta*theta))/
        ((1 - t)^(-(delta*theta)) + t^(-(delta*theta)))^
         ((1 + delta)/delta)))

})


####################################################
## Family BB6  #####################################

bb6.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 1 | theta < 1)
      stop("Parameter delta >= 1 and theta >=  1 for BB6 family")

    val <- new("bb6.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
            param.lowbnd = c(1,1), param.upbnd = c(Inf, Inf), 
            message = "BB6 copula family; Archimedean Copula") 
     val
}

setMethod("Phi","bb6.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   delta <- as.vector(copula@parameters[2])
 (-log(1 - (1 - t)^theta))^delta
 })

setMethod("PhiDer","bb6.copula",function(copula,t) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
    -((delta*(1 - t)^(-1 + theta)*theta*
      (-log(1 - (1 - t)^theta))^(-1 + delta))/
    (1 - (1 - t)^theta))
})

setMethod("InvPhi","bb6.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
     delta <- as.vector(copula@parameters[2])
   1 - (1 - exp(-s^(1/delta)))^(1/theta)
})

setMethod("InvPhisecondDer","bb6.copula",function(copula, s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])

    s1d <- exp(-s^delta^(-1))
  
  ((1 - s1d)^theta^(-1)*
    s^(-2 + delta^(-1))*
    ((-1 + delta)*(-1 + s1d)*theta + 
      s^delta^(-1)*(-1 + theta/s1d)))/
  (delta^2*(-1 + 1/s1d)^2*theta^2)
})

#####################################################################
####### Family BB7   ################################################

bb7.copula <- function(theta,delta) {
    if (is.na(delta) | is.na(theta) | delta < 0 | theta < 1)
      stop("Parameter delta > 0 and theta >=  1 for BB7 family")

   val <- new("bb7.copula", parameters = c(theta,delta), 
            param.names = c("theta","delta"),
            param.lowbnd = c(1,0), param.upbnd = c(Inf, Inf), 
            message = "BB7 copula family; Archimedean Copula") 
   val
}

setMethod("Phi","bb7.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
   delta <- as.vector(copula@parameters[2])
   -1 + (1 - (1 - t)^theta)^(-delta)
})

setMethod("PhiDer","bb7.copula", function(copula,t) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
   -(delta*(1 - (1 - t)^theta)^(-1 - delta)*
    (1 - t)^(-1 + theta)*theta)
})

setMethod("InvPhi","bb7.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
        1 - (1 - (1 + s)^(- 1/delta))^( 1/theta)
})

setMethod("InvPhifirstDer","bb7.copula", function(copula,s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
      - ( 1/theta)* (1 - (1 + s)^(- 1/delta))^( 1/theta - 1)*(1 + s)^(- 1/delta - 1)*(1/delta)
})

setMethod("InvPhisecondDer","bb7.copula", function(copula, s) {
     theta <- as.vector(copula@parameters[1])
        delta <- as.vector(copula@parameters[2])
((1 - (1 + s)^(-delta^(-1)))^theta^(-1)*
    (-1 + (1 + s)^delta^(-1)*theta +
      delta*(-1 + (1 + s)^delta^(-1))*theta))/
  (delta^2*(1 + s)^2*(-1 + (1 + s)^delta^(-1))^2*
    theta^2)
})

###################################################################
###   Joe Copula ##################################################

joe.copula <- function(theta) {
    if (  is.na(theta) | theta < 1)
      stop("Parameter delta >= 1 for Joe family")

   val <- new("joe.copula", parameters = c(theta), param.names = "theta",
       param.lowbnd = c(1), param.upbnd = c(Inf), 
         message = "Joe copula family; Archimedean Copula") 
   val
}

setMethod("Phi","joe.copula", function(copula,t) {
   theta <- as.vector(copula@parameters[1])
    (-log(1 - (1 - t)^theta))
 })

setMethod("PhiDer","joe.copula", function(copula,t) {
         theta <- as.vector(copula@parameters[1])
            -((1 - t)^(-1 + theta)*theta)/(1 - (1 - t)^theta)
})

setMethod("InvPhi","joe.copula", function(copula,s) {
    theta <- as.vector(copula@parameters[1])
    1 - (1 - exp(-s))^(1/theta)
})

setMethod("InvPhisecondDer","joe.copula", function(copula, s) {
     theta <- as.vector(copula@parameters[1])
     s1d <- exp(-s)

  ((1 - s1d)^theta^(-1)*
    s^(-2 + 1)*
    ( s*(-1 + theta/s1d)))/((-1 + 1/s1d)^2*theta^2)
})

###################
# Tail Index Computations

setMethod("tail.index","bb1.copula", function(copula,...) 
{   
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   LI <- 2^(-1/delta/theta)
    UI <- 2 - 2^(1/delta)
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
})

setMethod("tail.index","bb2.copula", function(copula,...) 
{
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   LI <- 1
    UI <- 0
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
})

 setMethod("tail.index","bb3.copula", function(copula,...) 
{   
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   if(theta == 1) {LI <- 2^(-1/delta)}
   else{LI <- 1}
    UI <- 2 - 2^(1/theta)
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
})

setMethod("tail.index","bb4.copula", function(copula,...) 
{   
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   LI <- (2 - 2^(-1/delta))^(-1/theta)
   UI <- 2^(-1/delta)
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
})

setMethod("tail.index","bb6.copula", function(copula,...) 
{   
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   LI <- 0
   UI <- 2 - 2^(1/delta/theta)
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
})

setMethod("tail.index","bb7.copula", function(copula,...) 
{   
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
   LI <- 2^(-1/delta)
   UI <- 2 - 2^(1/theta)
   val <- c(LI,UI)
   names(val) <- c("lower.tail","upper.tail")
   val
  })

setMethod("tail.index","empirical.copula", function(copula,...) 
{   
      n <- length(cop$x)
      alpha <- c(trunc(n/5):n)/n
      uind <- (1 - 2*alpha + pcopula(cop,alpha,alpha))/(1 - alpha)
      lind <- (pcopula(cop,1 - alpha,1 - alpha))/(1 - alpha)

      par(mfrow = c(2,1))
      plot(alpha,lind, ylab = "Lower Index")
      plot(1 - alpha,uind, ylab = "Upper Index", xlab = "alpha")
      par(mfrow = c(1,1))     
})


#################
# Functions for Comparing Multiple Copulas

cop.dist.L2.norm <- function(copula1, copula2, n = 100) {

    xmat <- rep(c(1:(n-1))/(n), n -1)
    ymat <- rep(c(1:(n-1))/n, each = n - 1) 
    zmat1 <- matrix(ncol = n - 1, nrow = n - 1, pcopula(copula1, xmat, ymat), byrow = FALSE)    
    zmat2 <- matrix(ncol = n - 1, nrow = n - 1, pcopula(copula2, xmat, ymat), byrow = FALSE)
    val <- sum((zmat2 - zmat1)^2)/n/n
    val 
}

A.dist.L1.norm <- function(copula1, copula2, n = 100) {

    inputPOS1 <- inherits(copula1,"empirical.copula") & 
                          hasMethod("A",class(copula2))
    inputPOS2 <- inherits(copula2,"empirical.copula") & 
                         hasMethod("A",class(copula1))
    if (!(inputPOS1 |  inputPOS2)) {
          warning("One of the copulas should be empirical copula, and another should be EV copula")
          return(NA)
    }
    
    t <- c(0:100)/100
    
    a1 <- A(copula1,t)
    a2 <- A(copula2,t)
    val <- sum(abs(a1 - a2))/100
    val 
}

Pearson.test <- function(empir.cop, param.copula, Ncells = 36, est.param = 0) 
{
    n <- floor(sqrt(Ncells + 1))
    data <- list(x = empir.cop@x, y = empir.cop@y)
    xll <- c(0:(n-1))/n
    xur <- c(1:n)/n
    yll <- c(0:(n-1))/n
    yur <- c(1:n)/n

    obs.x <- data$x
    obs.y <- data$y
    slistX <- sort.list(obs.x)
    slistY <- sort.list(obs.y)
    xsort <- as.double(obs.x[slistX])
    ysort <- as.double(obs.y[slistX])
    
    ind <- as.integer(match(slistX, slistY) - 1)
   N <- as.integer(n)
    nsamp <- as.integer(length(xsort))
        val1 <- as.double(rep(0, N*N))
    outC <- .C("jointDensity",
        xsort,
        ysort,
        nsamp,
       as.double(xur), as.double(yur),
        N,N,val1,as.double(1/n),as.double(1/n))
    val1 <- outC[[8]]
    emp.mat <- matrix(nrow = n, ncol = n, data = val1, byrow = TRUE)
    
    copula <- param.copula
    x <- rep(xll,n)
    y <- rep(yll,each = n)
    v00 <- pcopula(copula,x,y)
    v00[ x==0 | y == 0] <- 0
    v00[ x == 1] <- y[x ==1]
    v00[ y == 1] <- x[y ==1]
     
    x <- rep(xur,n)
    y <- rep(yll,each = n)
    v01 <- pcopula(copula,x,y)
    v01[ x==0 | y == 0] <- 0
    v01[ x == 1] <- y[x ==1]
    v01[ y == 1] <- x[y ==1]
     
    x <- rep(xll,n)
    y <- rep(yur,each = n)
    v10 <- pcopula(copula,x,y)
    v10[ x==0 | y == 0] <- 0
    v10[ x == 1] <- y[x ==1]
    v10[ y == 1] <- x[y ==1]
     
    x <- rep(xur,n)
    y <- rep(yur,each = n)
    v11 <- pcopula(copula,x,y)
    v11[ x==0 | y == 0] <- 0
    v11[ x == 1] <- y[x ==1]
    v11[ y == 1] <- x[y ==1]
     
    vcop <-  matrix(nrow = n, ncol = n, data = (v11 - v01 -v10 +v00), byrow = TRUE)
    
    #pull together all cells with p < 0.05!
    small.cells <- vcop < 0.005
    total.pr <- sum(vcop[small.cells])
    total.count <- sum(emp.mat[small.cells])
    
    #cat(sum(small.cells), total.pr)
    vcop[small.cells] <- 0
    emp.mat[small.cells] <- 0
    
    real.mat <-  vcop*nsamp     

    Chi <- sum((real.mat[!small.cells] - emp.mat[!small.cells])^2/real.mat[!small.cells])          
    
    if (total.pr > 0) Chi <- Chi + (total.count - total.pr *nsamp)^2/total.pr/nsamp
    
    Df <- as.integer(n*n) - 1 - est.param
    
    res <- list(Chi.Statistics = Chi, deegrees.freedom = Df, p.value = 1 - pchisq(Chi,Df))
    as.data.frame(res)
}

############################
# K estimators for copulas
##KEST

Wis <- function(x,y) {

   tf <- function(i,x,y) {
      tt <- x < x[i] & y < y[i]
      sum(tt)
   }
    
   ww <- matrix(ncol = length(x), nrow = 1, data = c(1:length(x)))
   ww <- apply(ww,2,tf,x = x, y = y)
   
  as.vector(ww)/(length(x) -1)
}

#lambda.archm.copula <- function(copula, t) {
#       Phi(copula,t)/PhiDer(copula,t)
#   }
    
    
K.Kern.est <- function(w,t, D) {
       b <- 1/length(w)*D
       tf <- function(t,w,b) {
           sum(pnorm(t-w,0, b))
       }

         ww <- matrix(ncol = length(t), nrow = 1, data = t)
        ww <- apply(ww,2,tf,w = w, b = b)
        ww[t ==0] <- 0
        ww[t ==1] <- length(w)         
      as.vector(ww)/(length(w))
    }
    
K.ast.Kern.est <- function(w,t, D) {
       b <- 1/length(w)*D
       tf <- function(t,w,b) {
           sum(pnorm(t-w,0, b))
       }

         ww <- matrix(ncol = length(t), nrow = 1, data = t)
        ww <- apply(ww,2,tf,w = w, b = b)
        ww[t ==0] <- 0
        ww[t ==0] <- 1
        val <- as.vector(ww)/(length(w))
        cummax(val)
    }
    
setMethod("lambda","empirical.copula", function(copula, t) {
         w <- Wis(copula@x, copula@y)
         t -  K.Kern.est(w,t, 0.5)
})




########################## 
### Functions dcopula  ### 


setMethod("dcopula","normal.copula",  
function(copula, u, v) {
       delta <- copula@parameters[1]
	   if (is.na(delta) | delta >= 1 | delta < 0)
          stop("Parameter of the bivariate the normal.copula should be in [0,1)")

      x <- qnorm(u)
       y <- qnorm(v)
       val <- (1 - delta^2)^(-1/2)*exp(-1/2*(1 - delta^2)^(-1)*
                      (x^2 + y^2 - 2* delta*x*y))*exp(1/2*(x^2 + y^2))
       val[ u <= 0 | v <= 0 ] <- 0
       val[ u >= 1 | v >= 1] <- 0
      val
} )


setMethod("dcopula","normal.mix.copula", 
function(copula, u, v) {
    p <- as.vector(copula@parameters[1])
    delta1 <- as.vector(copula@parameters[2])
   delta2 <- as.vector(copula@parameters[3])

      x <- qnorm(u)
       y <- qnorm(v)
       val1 <- (1 - delta1^2)^(-1/2)*exp(-1/2*(1 - delta1^2)^(-1)*
                      (x^2 + y^2 - 2* delta1*x*y))*exp(1/2*(x^2 + y^2)) 
    
       val2 <- (1 - delta2^2)^(-1/2)*exp(-1/2*(1 - delta2^2)^(-1)*
                      (x^2 + y^2 - 2* delta2*x*y))*exp(1/2*(x^2 + y^2)) 
    
       val <- val1*p + (1 - p)*val2
    
       val[ u <= 0 | v <= 0 ] <- 0
       val[ u >= 1 | v >= 1] <- 0
      
      val
})

setMethod("dcopula","kimeldorf.sampson.copula",
function(copula, u, v) {
   delta <- as.vector(copula@parameters[1])
   (1 + delta)*(u*v)^(-(delta + 1))*
       (u^(-delta) + v^(-delta) -1) ^(-2-1/delta)
})


setMethod("dcopula","husler.reiss.copula", 
function(copula, u, v) {
    d <- as.vector(copula@parameters[1])
    z <- log(u)/log(v)

    ff <- pnorm(1/d + 1/2*d*log(1/z)) *  pnorm(1/d + 1/2*d*log(z)) 
   
   1/u/v * pcopula(copula, u,v) * (ff + 1/2 * d/(-log(v))*dnorm(1/d + 1/2*d*log(z)))    
})


setMethod("dcopula","frank.copula", 
function(copula, u, v) {
     delta <- as.vector(copula@parameters[1])

     eta <- 1 - exp(-delta)
     delta *eta *exp( - delta *(u+v))/(eta 
        -(1 - exp(-delta*u))*(1 - exp(-delta*v)))^2   
})


setMethod("dcopula","bb4.copula", 
function(copula, u, v) {
       pu <- Phi(copula,u)
       pv <- Phi(copula,v)

       t <- pu/(pu+pv)
       ad0 <- A(copula, t)
       ad1 <- AfirstDer(copula, t)
       ad2 <- AsecondDer(copula, t)
       t2 <- ad0*(pu + pv)
       val <-  ((-(ad2*InvPhifirstDer(copula,t2)*pu*pv) + 
    InvPhisecondDer(copula,t2)*(pu+pv)*
       (ad0*(ad0 - ad1)*pu^2 + 
         (2*ad0^2 - ad1^2)*pu*pv + 
         ad0*(ad0 + ad1)*pv^2))*
    PhiDer(copula,u)*PhiDer(copula,v))/
  (pu+pv)^3

       val[is.na(val)] <- 0
       val
})



setMethod("dcopula","gumbel.copula",
function(copula, u, v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    t <- log(x)/log(x * y)
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    ad2 <- AsecondDer(copula, t)
    lxy <- log(x * y)
    lx <- log(x)
    val <- ad0^2/(x * y) + (ad0 * ad1)/(x * y) + (ad2 * lx^2)/(x * y * lxy^3) - (ad2 * lx)/(x * y * lxy^2) +
        (ad1^2 * lx^2)/(x * y * lxy^2) - (2 * ad0 * ad1 * lx)/(x * y * lxy) - (ad1^2 * lx)/(x * y * lxy)
    val <- val * pcopula(copula, x, y)
    val[u == 0 | v == 0 | v == 1 | u == 1] <- 0
    #   val[is.na(val)] <- NA
    val
} )

setMethod("dcopula","galambos.copula",
function(copula, u, v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    t <- log(x)/log(x * y)
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    ad2 <- AsecondDer(copula, t)
    lxy <- log(x * y)
    lx <- log(x)
    val <- ad0^2/(x * y) + (ad0 * ad1)/(x * y) + (ad2 * lx^2)/(x * y * lxy^3) - (ad2 * lx)/(x * y * lxy^2) +
        (ad1^2 * lx^2)/(x * y * lxy^2) - (2 * ad0 * ad1 * lx)/(x * y * lxy) - (ad1^2 * lx)/(x * y * lxy)
    val <- val * pcopula(copula, x, y)
    val[u == 0 | v == 0 | v == 1 | u == 1] <- 0
    #   val[is.na(val)] <- NA
    val
} )

setMethod("dcopula","tawn.copula",
function(copula, u, v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    t <- log(x)/log(x * y)
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    ad2 <- AsecondDer(copula, t)
    lxy <- log(x * y)
    lx <- log(x)
    val <- ad0^2/(x * y) + (ad0 * ad1)/(x * y) + (ad2 * lx^2)/(x * y * lxy^3) - (ad2 * lx)/(x * y * lxy^2) +
        (ad1^2 * lx^2)/(x * y * lxy^2) - (2 * ad0 * ad1 * lx)/(x * y * lxy) - (ad1^2 * lx)/(x * y * lxy)
    val <- val * pcopula(copula, x, y)
    val[u == 0 | v == 0 | v == 1 | u == 1] <- 0
    #   val[is.na(val)] <- NA
    val
} )

setMethod("dcopula","bb5.copula",
function(copula, u, v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    t <- log(x)/log(x * y)
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    ad2 <- AsecondDer(copula, t)
    lxy <- log(x * y)
    lx <- log(x)
    val <- ad0^2/(x * y) + (ad0 * ad1)/(x * y) + (ad2 * lx^2)/(x * y * lxy^3) - (ad2 * lx)/(x * y * lxy^2) +
        (ad1^2 * lx^2)/(x * y * lxy^2) - (2 * ad0 * ad1 * lx)/(x * y * lxy) - (ad1^2 * lx)/(x * y * lxy)
    val <- val * pcopula(copula, x, y)
    val[u == 0 | v == 0 | v == 1 | u == 1] <- 0
    #   val[is.na(val)] <- NA
    val
} )


setMethod("dcopula","bb1.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )

setMethod("dcopula","bb2.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )

setMethod("dcopula","bb3.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )

setMethod("dcopula","bb6.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )

setMethod("dcopula","bb7.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )

setMethod("dcopula","joe.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    val <- InvPhisecondDer(copula, val)
    val <- val * PhiDer(copula, u) * PhiDer(copula, v)
    val
} )


########################## 
### Functions pcopula  ### 
########################## 


setMethod("pcopula","normal.copula", 
function(copula, u, v) {
       delta <- as.vector(copula@parameters[1])
      if (is.na(delta) | delta >= 1 | delta < 0)
          stop("Parameter of the bivariate the normal.copula should be in [0,1)")

      good <- u < 1 & u > 0 & v < 1 & v > 0      
      ug <- u[good]
      vg <- v[good]      
      if (sum(good) > 0) 
      {
      	tlen <- max(length(ug),length(vg))
       	x <- rep(qnorm(ug),length.out = tlen)
      	 y <- rep(qnorm(vg),length.out = tlen)           
       	XX <- matrix(ncol = 2, nrow= length(x), data = c(x,y), byrow = FALSE)     
#       valg <- p2Dnorm(XX, rho = delta)
       	ff <- function(x) p2Dnorm(x,rho=delta)
       	valg <- apply(XX,1,ff)    
       }
       val <- vector(length =   length(good), mode = "numeric")
       val[1:length(good)] <- NA    
       if (sum(good) > 0) { val[good] <- valg}
       val[ u <= 0 | v <= 0 ] <- 0
       val[ u >= 1 & v >= 1] <- 1
       val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
       val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
       val
} )


setMethod("pcopula","frank.copula", 
function(copula, u, v) {
   delta <- as.vector(copula@parameters[1])
   eta <- 1 - exp(-delta)
   -1/delta*log((eta - (1 - exp(-delta*u))*(1 - exp(-delta*v)))/eta)
})


setMethod("pcopula","kimeldorf.sampson.copula", 
function(copula, u, v) {
   delta <- as.vector(copula@parameters[1])
   (u^(-delta) + v^(-delta) -1) ^(-1/delta)
})

setMethod("pcopula","normal.mix.copula", 
function(copula, u, v) {
    p <- as.vector(copula@parameters[1])
    delta1 <- as.vector(copula@parameters[2])
   delta2 <- as.vector(copula@parameters[3])

      tlen <- max(length(u),length(v))
      u <- rep(u,length.out = tlen)
      v <- rep(v,length.out = tlen)
    
      good <- u < 1 & u > 0 & v < 1 & v > 0      
      ug <- u[good]        
      vg <- v[good]      

      x <- qnorm(ug)       
      y <- qnorm(vg)     
      
       XX <- matrix(ncol = 2, nrow= length(x), data = c(x,y), byrow = FALSE)
       
     #  cat(ug,"\n")       
     #  cat(vg,"\n")       
     #  cat(sum(is.na(XX)),"\n")

#       valg1 <- p2Dnorm(XX, rho = delta1)           
       ff <- function(x) p2Dnorm(x,rho=delta1)
       valg1 <- apply(XX,1,ff)    
#       valg2 <- p2Dnorm(XX, rho = delta2)
       ff <- function(x) p2Dnorm(x,rho=delta2)
       valg2 <- apply(XX,1,ff)    
       
       valg <- p * valg1 + (1 -p) * valg2
    
       val <- vector(length = length(good), mode = "numeric")
       val[1:length(good)] <- NA    
       val[good] <- valg
       val[ u <= 0 | v <= 0 ] <- 0
       val[ u >= 1 & v >= 1] <- 1
       val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
       val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
       val
})


setMethod("pcopula","gumbel.copula",
function(copula, u, v) {
    t <- log(u)/log(u * v)
    val <- exp(log(u * v) * A(copula, t))
    val[u <= 0 | v <= 0] <- 0
    val[u >= 1 & v >= 1] <- 1
    val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
    val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
    val
} )

setMethod("pcopula","galambos.copula",
function(copula, u, v) {
    t <- log(u)/log(u * v)
    val <- exp(log(u * v) * A(copula, t))
    val[u <= 0 | v <= 0] <- 0
    val[u >= 1 & v >= 1] <- 1
    val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
    val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
    val
} )


setMethod("pcopula","husler.reiss.copula",
function(copula, u, v) {
  d <- as.vector(copula@parameters[1])
  uu <- -log(u)
  vv <- -log(v)
  ft1 <- pnorm(1/d + 1/2*d*log(uu/vv))
  ft2 <- pnorm(1/d + 1/2*d*log(vv/uu))

  exp(-uu*ft1 - vv *ft2)
 })

setMethod("pcopula","tawn.copula",
function(copula, u, v) {
    t <- log(u)/log(u * v)
    val <- exp(log(u * v) * A(copula, t))
    val[u <= 0 | v <= 0] <- 0
    val[u >= 1 & v >= 1] <- 1
    val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
    val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
    val
} )

setMethod("pcopula","bb5.copula",
function(copula, u, v) {
    t <- log(u)/log(u * v)
    val <- exp(log(u * v) * A(copula, t))
    val[u <= 0 | v <= 0] <- 0
    val[u >= 1 & v >= 1] <- 1
    val[is.na(val) & u >= 1] <- v[is.na(val) & u >= 1]
    val[is.na(val) & v >= 1] <- u[is.na(val) & v >= 1]
    val
} )

setMethod("pcopula","bb1.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","bb2.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","bb3.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","bb6.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","bb7.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","joe.copula",
function(copula, u, v) {
    val <- Phi(copula, u) + Phi(copula, v)
    InvPhi(copula, val)
} )

setMethod("pcopula","bb4.copula",
function(copula, u,v) {
      pu <- Phi(copula,u)
      pv <- Phi(copula,v)
      val <- A(copula,pu/(pu+pv))
     InvPhi(copula,(pu+pv)*val)
})

       
setMethod("pcopula","empirical.copula", 
function (copula,u,v) 
{
data <- copula
    obs.x <- data@x
    obs.y <- data@y
x.mat <- outer(obs.x,as.vector(u),"<=")
y.mat <- outer(obs.y,as.vector(v),"<=")
    val <- apply(x.mat&y.mat,2,mean)
val
})

#function(copula, u, v) 
#{
#    data <- copula
#    obs.x <- c(0,data@x,1)
#    obs.y <- c(0,data@y,1)
#    slistX <- sort.list(obs.x)
#   slistY <- sort.list(obs.y)
#   xsort <- as.double(obs.x[slistX])
#   ysort <- as.double(obs.y[slistY])
#   ind <- as.integer(match(slistX,slistY) - 1)
#   nsamp <- as.integer(length(xsort))

#   N.in <- length(u)
#    good.in.values <- (u >= 0) & (u <= 1) & (v >= 0) & (v <= 1)
#    val <- vector (length = N.in, mode = "numeric")
#    val[!good.in.values] <- NA
#    xx <- as.double(u[good.in.values])
#    yy <- as.double(v[good.in.values]) 
#    N <- as.integer(length(xx))
#    val1 <- as.double(rep(-1,N))
       
#   outC <- .C("interpcopR", xsort,ysort,ind, nsamp,xx,yy,
#                N,val1)
#   val[good.in.values] <- outC[[8]]
#    val
#})


########################## 
### Functions rcopula  ### 
########################## 


setMethod("rcopula", "normal.copula" , function(copula, n) {
       delta <-  copula@parameters[1]
	   if (is.na(delta) | delta >= 1 | delta < 0)
          stop("Parameter of the bivariate the normal.copula should be in [0,1)")

       bvsn <- rmvnorm(n,cov=matrix(c(1,delta,delta,1),byrow=TRUE,ncol=2))
       x <- pnorm(bvsn[,1])
       y <- pnorm(bvsn[,2])
       sim <- list(x = x, y = y)
       sim
})

setMethod("rcopula","frank.copula", function(copula, n) {
   # inverse EXISTS!!
   Usamp <- runif(n)
#   Vsamp <- CondCinverse.b3.copula(copula,runif(n),Usamp)
   Vsamp <- CondCinverse(copula,runif(n),Usamp)
   val <- list(x = Usamp, y = Vsamp)
   val
})

setMethod("rcopula","kimeldorf.sampson.copula", function(copula, n) {
   # inverse EXISTS!!
   Usamp <- runif(n)
   Vsamp <- CondCinverse(copula,runif(n),Usamp)
   val <- list(x = Usamp, y = Vsamp)
   val
})

setMethod("rcopula","gumbel.copula",
     function(copula, n) {
    # generate z from distribution h
    # using rejection method
    zrand <- vector(length = n, mode = "numeric")
    xgenerated <- vector(length = n, mode = "logical")
    lg <- n
    #set c to 1/2 true for symmetric copula
    cc <- Hderiv(copula, 1/2)
    stopnow <- FALSE
    while(!stopnow) {
        usamp <- runif(lg)
        ysamp <- runif(lg)
        ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
        #indexes of xrand vector to store "good" generated values
        toind <- (c(1:n)[!xgenerated])[ykeep]
        zrand[toind] <- ysamp[ykeep]
        xgenerated[toind] <- TRUE
        #number of rvs left to generate:
        lg <- n - sum(xgenerated)
        if(lg == 0)
            stopnow <- TRUE
    }
    z <- zrand
    u <- runif(n)
    # cat(length(z), length(AsecondDer(copula,z)))
    pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
    w <- NULL
    w[1:n] <- 0
    nn <- sum(u <= pz)
    w[u <= pz] <- runif(nn)
    w[u > pz] <- runif(n - nn) * runif(n - nn)
    xx <- exp((z * log(w))/A(copula, z))
    yy <- exp(((1 - z) * log(w))/A(copula, z))
    val <- list(x = xx, y = yy)
    val
} )

setMethod("rcopula","galambos.copula",
     function(copula, n) {
    # generate z from distribution h
    # using rejection method
    zrand <- vector(length = n, mode = "numeric")
    xgenerated <- vector(length = n, mode = "logical")
    lg <- n
    #set c to 1/2 true for symmetric copula
    cc <- Hderiv(copula, 1/2)
    stopnow <- FALSE
    while(!stopnow) {
        usamp <- runif(lg)
        ysamp <- runif(lg)
        ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
        #indexes of xrand vector to store "good" generated values
        toind <- (c(1:n)[!xgenerated])[ykeep]
        zrand[toind] <- ysamp[ykeep]
        xgenerated[toind] <- TRUE
        #number of rvs left to generate:
        lg <- n - sum(xgenerated)
        if(lg == 0)
            stopnow <- TRUE
    }
    z <- zrand
    u <- runif(n)
    # cat(length(z), length(AsecondDer(copula,z)))
    pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
    w <- NULL
    w[1:n] <- 0
    nn <- sum(u <= pz)
    w[u <= pz] <- runif(nn)
    w[u > pz] <- runif(n - nn) * runif(n - nn)
    xx <- exp((z * log(w))/A(copula, z))
    yy <- exp(((1 - z) * log(w))/A(copula, z))
    val <- list(x = xx, y = yy)
    val
} )

setMethod("rcopula","husler.reiss.copula",
     function(copula, n) {
    # generate z from distribution h
    # using rejection method
    zrand <- vector(length = n, mode = "numeric")
    xgenerated <- vector(length = n, mode = "logical")
    lg <- n
    #set c to 1/2 true for symmetric copula
    cc <- Hderiv(copula, 1/2)
    stopnow <- FALSE
    while(!stopnow) {
        usamp <- runif(lg)
        ysamp <- runif(lg)
        ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
        #indexes of xrand vector to store "good" generated values
        toind <- (c(1:n)[!xgenerated])[ykeep]
        zrand[toind] <- ysamp[ykeep]
        xgenerated[toind] <- TRUE
        #number of rvs left to generate:
        lg <- n - sum(xgenerated)
        if(lg == 0)
            stopnow <- TRUE
    }
    z <- zrand
    u <- runif(n)
    # cat(length(z), length(AsecondDer(copula,z)))
    pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
    w <- NULL
    w[1:n] <- 0
    nn <- sum(u <= pz)
    w[u <= pz] <- runif(nn)
    w[u > pz] <- runif(n - nn) * runif(n - nn)
    xx <- exp((z * log(w))/A(copula, z))
    yy <- exp(((1 - z) * log(w))/A(copula, z))
    val <- list(x = xx, y = yy)
    val
} )

setMethod("rcopula","tawn.copula",
     function(copula, n) {
    # generate z from distribution h
    # using rejection method
    zrand <- vector(length = n, mode = "numeric")
    xgenerated <- vector(length = n, mode = "logical")
    lg <- n
    #set c to 1/2 true for symmetric copula
    cc <- Hderiv(copula, 1/2)
    stopnow <- FALSE
    while(!stopnow) {
        usamp <- runif(lg)
        ysamp <- runif(lg)
        ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
        #indexes of xrand vector to store "good" generated values
        toind <- (c(1:n)[!xgenerated])[ykeep]
        zrand[toind] <- ysamp[ykeep]
        xgenerated[toind] <- TRUE
        #number of rvs left to generate:
        lg <- n - sum(xgenerated)
        if(lg == 0)
            stopnow <- TRUE
    }
    z <- zrand
    u <- runif(n)
    # cat(length(z), length(AsecondDer(copula,z)))
    pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
    w <- NULL
    w[1:n] <- 0
    nn <- sum(u <= pz)
    w[u <= pz] <- runif(nn)
    w[u > pz] <- runif(n - nn) * runif(n - nn)
    xx <- exp((z * log(w))/A(copula, z))
    yy <- exp(((1 - z) * log(w))/A(copula, z))
    val <- list(x = xx, y = yy)
    val
} )

setMethod("rcopula","bb5.copula",
     function(copula, n) {
    # generate z from distribution h
    # using rejection method
    zrand <- vector(length = n, mode = "numeric")
    xgenerated <- vector(length = n, mode = "logical")
    lg <- n
    #set c to 1/2 true for symmetric copula
    cc <- Hderiv(copula, 1/2)
    stopnow <- FALSE
    while(!stopnow) {
        usamp <- runif(lg)
        ysamp <- runif(lg)
        ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
        #indexes of xrand vector to store "good" generated values
        toind <- (c(1:n)[!xgenerated])[ykeep]
        zrand[toind] <- ysamp[ykeep]
        xgenerated[toind] <- TRUE
        #number of rvs left to generate:
        lg <- n - sum(xgenerated)
        if(lg == 0)
            stopnow <- TRUE
    }
    z <- zrand
    u <- runif(n)
    # cat(length(z), length(AsecondDer(copula,z)))
    pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
    w <- NULL
    w[1:n] <- 0
    nn <- sum(u <= pz)
    w[u <= pz] <- runif(nn)
    w[u > pz] <- runif(n - nn) * runif(n - nn)
    xx <- exp((z * log(w))/A(copula, z))
    yy <- exp(((1 - z) * log(w))/A(copula, z))
    val <- list(x = xx, y = yy)
    val
} )


setMethod("rcopula","bb1.copula",
     function(copula, n) {
    #figure out family:
    famstr <- substring(class(copula)[1], 1, 3)
    implem.fam <- c("bb1", "bb2", "bb3", "bb6", "bb7")
    famID <- c(1, 2, 3, 6, 7)
    ind <- charmatch(famstr, implem.fam)
    family <- as.integer(famID[ind])
    #use Genest MacKay 1996algorithm
    u <- runif(n)
    t <- runif(n)
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
    # use C code to compute the inverse of secon deriv. of Phi
    x <- PhiDer(copula, u)/t
    # this is temporary!!!
    infx <- is.infinite(x)
    while(sum(infx) > 0) {
        warning(paste(" The derivative of Phi is Infinity in ", sum(infx), " points",
            " Simulations results maybe unreliable", sep = ""))
        tinf <- sum(infx)
        x[infx] <- PhiDer(copula, runif(tinf)/runif(tinf))
        infx <- is.infinite(x)
    }
    Xtemp <- as.double(sort(x))
    result <- vector(length = n)
    result[1:n] <- as.double(-999)
    epsX <- as.double(1e-010)
    epsY <- as.double(1e-010)
    conv <- as.integer(rep(-9, n))
    maxit <- as.integer(100)
    cout <- .C("inverseDerivPhiBB",
        Xtemp,
        as.integer(n),
        delta,
        theta,
        result,
        family,
        epsX,
        epsY,
        maxit,
        conv)
    result <- cout[[5]]
    conv <- cout[[10]]
    if(sum(conv <= 0) != 0)
        warning("Inverting derivative of Phi is unsuccessful. Simulation results might be unreliable")
    w <- result
    w[sort.list(x)] <- result
    Vtemp <- Phi(copula, w) - Phi(copula, u)
    v <- InvPhi(copula, Vtemp)
    val <- list(x = u, y = v)
    val
} )

setMethod("rcopula","bb2.copula",
     function(copula, n) {
    #figure out family:
    famstr <- substring(class(copula)[1], 1, 3)
    implem.fam <- c("bb1", "bb2", "bb3", "bb6", "bb7")
    famID <- c(1, 2, 3, 6, 7)
    ind <- charmatch(famstr, implem.fam)
    family <- as.integer(famID[ind])
    #use Genest MacKay 1996algorithm
    u <- runif(n)
    t <- runif(n)
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
    # use C code to compute the inverse of secon deriv. of Phi
    x <- PhiDer(copula, u)/t
    # this is temporary!!!
    infx <- is.infinite(x)
    while(sum(infx) > 0) {
        warning(paste(" The derivative of Phi is Infinity in ", sum(infx), " points",
            " Simulations results maybe unreliable", sep = ""))
        tinf <- sum(infx)
        x[infx] <- PhiDer(copula, runif(tinf)/runif(tinf))
        infx <- is.infinite(x)
    }
    Xtemp <- as.double(sort(x))
    result <- vector(length = n)
    result[1:n] <- as.double(-999)
    epsX <- as.double(1e-010)
    epsY <- as.double(1e-010)
    conv <- as.integer(rep(-9, n))
    maxit <- as.integer(100)
    cout <- .C("inverseDerivPhiBB",
        Xtemp,
        as.integer(n),
        delta,
        theta,
        result,
        family,
        epsX,
        epsY,
        maxit,
        conv)
    result <- cout[[5]]
    conv <- cout[[10]]
    if(sum(conv <= 0) != 0)
        warning("Inverting derivative of Phi is unsuccessful. Simulation results might be unreliable")
    w <- result
    w[sort.list(x)] <- result
    Vtemp <- Phi(copula, w) - Phi(copula, u)
    v <- InvPhi(copula, Vtemp)
    val <- list(x = u, y = v)
    val
} )

setMethod("rcopula","bb3.copula",
     function(copula, n) {
    #figure out family:
    famstr <- substring(class(copula)[1], 1, 3)
    implem.fam <- c("bb1", "bb2", "bb3", "bb6", "bb7")
    famID <- c(1, 2, 3, 6, 7)
    ind <- charmatch(famstr, implem.fam)
    family <- as.integer(famID[ind])
    #use Genest MacKay 1996algorithm
    u <- runif(n)
    t <- runif(n)
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
    # use C code to compute the inverse of secon deriv. of Phi
    x <- PhiDer(copula, u)/t
    # this is temporary!!!
    infx <- is.infinite(x)
    while(sum(infx) > 0) {
        warning(paste(" The derivative of Phi is Infinity in ", sum(infx), " points",
            " Simulations results maybe unreliable", sep = ""))
        tinf <- sum(infx)
        x[infx] <- PhiDer(copula, runif(tinf)/runif(tinf))
        infx <- is.infinite(x)
    }
    Xtemp <- as.double(sort(x))
    result <- vector(length = n)
    result[1:n] <- as.double(-999)
    epsX <- as.double(1e-010)
    epsY <- as.double(1e-010)
    conv <- as.integer(rep(-9, n))
    maxit <- as.integer(100)
    cout <- .C("inverseDerivPhiBB",
        Xtemp,
        as.integer(n),
        delta,
        theta,
        result,
        family,
        epsX,
        epsY,
        maxit,
        conv)
    result <- cout[[5]]
    conv <- cout[[10]]
    if(sum(conv <= 0) != 0)
        warning("Inverting derivative of Phi is unsuccessful. Simulation results might be unreliable")
    w <- result
    w[sort.list(x)] <- result
    Vtemp <- Phi(copula, w) - Phi(copula, u)
    v <- InvPhi(copula, Vtemp)
    val <- list(x = u, y = v)
    val
} )

setMethod("rcopula","bb6.copula",
     function(copula, n) {
    #figure out family:
    famstr <- substring(class(copula)[1], 1, 3)
    implem.fam <- c("bb1", "bb2", "bb3", "bb6", "bb7")
    famID <- c(1, 2, 3, 6, 7)
    ind <- charmatch(famstr, implem.fam)
    family <- as.integer(famID[ind])
    #use Genest MacKay 1996algorithm
    u <- runif(n)
    t <- runif(n)
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
    # use C code to compute the inverse of secon deriv. of Phi
    x <- PhiDer(copula, u)/t
    # this is temporary!!!
    infx <- is.infinite(x)
    while(sum(infx) > 0) {
        warning(paste(" The derivative of Phi is Infinity in ", sum(infx), " points",
            " Simulations results maybe unreliable", sep = ""))
        tinf <- sum(infx)
        x[infx] <- PhiDer(copula, runif(tinf)/runif(tinf))
        infx <- is.infinite(x)
    }
    Xtemp <- as.double(sort(x))
    result <- vector(length = n)
    result[1:n] <- as.double(-999)
    epsX <- as.double(1e-010)
    epsY <- as.double(1e-010)
    conv <- as.integer(rep(-9, n))
    maxit <- as.integer(100)
    cout <- .C("inverseDerivPhiBB",
        Xtemp,
        as.integer(n),
        delta,
        theta,
        result,
        family,
        epsX,
        epsY,
        maxit,
        conv)
    result <- cout[[5]]
    conv <- cout[[10]]
    if(sum(conv <= 0) != 0)
        warning("Inverting derivative of Phi is unsuccessful. Simulation results might be unreliable")
    w <- result
    w[sort.list(x)] <- result
    Vtemp <- Phi(copula, w) - Phi(copula, u)
    v <- InvPhi(copula, Vtemp)
    val <- list(x = u, y = v)
    val
} )

setMethod("rcopula","bb7.copula",
     function(copula, n) {
    #figure out family:
    famstr <- substring(class(copula)[1], 1, 3)
    implem.fam <- c("bb1", "bb2", "bb3", "bb6", "bb7")
    famID <- c(1, 2, 3, 6, 7)
    ind <- charmatch(famstr, implem.fam)
    family <- as.integer(famID[ind])
    #use Genest MacKay 1996algorithm
    u <- runif(n)
    t <- runif(n)
    theta <- as.vector(copula@parameters[1])
    delta <- as.vector(copula@parameters[2])
    # use C code to compute the inverse of secon deriv. of Phi
    x <- PhiDer(copula, u)/t
    # this is temporary!!!
    infx <- is.infinite(x)
    while(sum(infx) > 0) {
        warning(paste(" The derivative of Phi is Infinity in ", sum(infx), " points",
            " Simulations results maybe unreliable", sep = ""))
        tinf <- sum(infx)
        x[infx] <- PhiDer(copula, runif(tinf)/runif(tinf))
        infx <- is.infinite(x)
    }
    Xtemp <- as.double(sort(x))
    result <- vector(length = n)
    result[1:n] <- as.double(-999)
    epsX <- as.double(1e-010)
    epsY <- as.double(1e-010)
    conv <- as.integer(rep(-9, n))
    maxit <- as.integer(100)
    cout <- .C("inverseDerivPhiBB",
        Xtemp,
        as.integer(n),
        delta,
        theta,
        result,
        family,
        epsX,
        epsY,
        maxit,
        conv)
    result <- cout[[5]]
    conv <- cout[[10]]
    if(sum(conv <= 0) != 0)
        warning("Inverting derivative of Phi is unsuccessful. Simulation results might be unreliable")
    w <- result
    w[sort.list(x)] <- result
    Vtemp <- Phi(copula, w) - Phi(copula, u)
    v <- InvPhi(copula, Vtemp)
    val <- list(x = u, y = v)
    val
} )

setMethod("rcopula","joe.copula", function(copula, n) {
       theta <- as.vector(copula@parameters[1])
       rcopula(bb6.copula(theta,1.0), n)
    })

    
setMethod("rcopula","bb4.copula", function(copula, n) {
   # generate z from distribution h
   # using rejection method

   zrand <- vector(length = n, mode = "numeric")
   xgenerated <- vector(length = n, mode = "logical")
   lg <- n
   #set c to 1/2 true for symmetric copula
   cc <- Hderiv(copula,1/2)
   stopnow <- FALSE
   while(!stopnow) {
     usamp <- runif(lg)
     ysamp <- runif(lg)
     ykeep <- (usamp <= Hderiv(copula,ysamp)/cc)
     #indexes of xrand vector to store "good" generated values
     toind <- (c(1:n)[!xgenerated])[ykeep]
     zrand[toind] <- ysamp[ykeep]
     xgenerated[toind] <- TRUE
     #number of rvs left to generate:
     lg <- n - sum(xgenerated)
     if (lg == 0)  stopnow <- TRUE
   }
   
   z <- zrand 
   u <- runif(n)
   
 # cat(length(z), length(AsecondDer(copula,z)))
   pz <- z*(1-z)*AsecondDer(copula,z)/Hderiv(copula,z)/A(copula,z)

   w <- NULL
   w[1:n] <- 0
   nn <- sum(u <=  pz)
   w[u <=  pz] <- runif(nn)

   # we need to generate n - nn r.v from 
   # archimedian copula with phi
      PHIU <- runif(n - nn)
      PHIT <- runif(n - nn)
      PHIW <- DerPhi.minus1(copula, PhiDer(copula,PHIU)/PHIT)

   w[u > pz] <- PHIW
   
   pw <- Phi(copula,w)
   xx <- InvPhi(copula, z*pw/A(copula,z))
   yy <- InvPhi(copula, (1 - z)*pw/A(copula,z))
   
   val <- list(x = xx, v = yy)
   val   
})


setMethod("rcopula","normal.mix.copula", function(copula, n) {
          p <- as.vector(copula@parameters[1])
    delta1 <- as.vector(copula@parameters[2])
   delta2 <- as.vector(copula@parameters[3])

       bvsn1 <- rmvnorm(n, d = 2, rho = delta1)
       bvsn2 <- rmvnorm(n, d = 2, rho = delta2)

       uu <- runif(n)
      
      bvsn <- bvsn1
      bvsn[uu > p,] <- bvsn2[uu > p,]      

       x <- pnorm(bvsn[,1])
       y <- pnorm(bvsn[,2])
       sim <- list(x = x, y = y)
       sim  
})



######################### 
### Functions Hderiv  ### 
#########################

setMethod("Hderiv","galambos.copula",
     function(copula, z) {
    ad0 <- A(copula, z)
    ad1 <- AfirstDer(copula, z)
    ad2 <- AsecondDer(copula, z)
    1 + ((1 - 2 * z) * ad1)/ad0 + z * (1 - z) * (ad2/ad0 - ad1^2/ad0^2)
} )

setMethod("Hderiv","husler.reiss.copula",
     function(copula, z) {
    ad0 <- A(copula, z)
    ad1 <- AfirstDer(copula, z)
    ad2 <- AsecondDer(copula, z)
    1 + ((1 - 2 * z) * ad1)/ad0 + z * (1 - z) * (ad2/ad0 - ad1^2/ad0^2)
} )

setMethod("Hderiv","bb5.copula",
     function(copula, z) {
    ad0 <- A(copula, z)
    ad1 <- AfirstDer(copula, z)
    ad2 <- AsecondDer(copula, z)
    1 + ((1 - 2 * z) * ad1)/ad0 + z * (1 - z) * (ad2/ad0 - ad1^2/ad0^2)
} )




############################### 
### Functions Kendalls.tau  ### 
###############################

setMethod("Kendalls.tau","gumbel.copula", function(copula, tol = 1e-5)
 {
    f <- function(t, copula)
    {
        (t * (1 - t) * AsecondDer(copula, t))/A(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    t1
} )

setMethod("Kendalls.tau","galambos.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        (t * (1 - t) * AsecondDer(copula, t))/A(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    t1
} )

setMethod("Kendalls.tau","husler.reiss.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        (t * (1 - t) * AsecondDer(copula, t))/A(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    t1
} )

setMethod("Kendalls.tau","tawn.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        (t * (1 - t) * AsecondDer(copula, t))/A(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    t1
} )

setMethod("Kendalls.tau","bb5.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        (t * (1 - t) * AsecondDer(copula, t))/A(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    t1
} )

setMethod("Kendalls.tau","normal.copula", function(copula, tol = 1e-5) {
	  delta <- as.vector(copula@parameters[1])
     return(2*asin(delta)/pi)	
})

setMethod("Kendalls.tau","kimeldorf.sampson.copula", function(copula, tol = 1e-5) {
	  delta <- as.vector(copula@parameters[1])
     return(delta/(delta+2))	
})
 
setMethod("Kendalls.tau","bb1.copula", function(copula, tol = 1e-5) {
	  theta <- as.vector(copula@parameters[1])
	  delta <- as.vector(copula@parameters[2])
	  return(1 - 2/delta/(theta + 2))	
})

setMethod("Kendalls.tau","bb6.copula", function(copula, tol = 1e-5) {
	  delta <- as.vector(copula@parameters[1])
     return((delta-1)/(delta))	
})

setMethod("Kendalls.tau","bb4.copula", function(copula, tol = 1e-5){
	  theta <- as.vector(copula@parameters[1])
	  delta <- as.vector(copula@parameters[2])
	  tauphi <- 1 - 4/(4 + 2*theta)
	   
	   f <- function(t, copula) 
	   {	   
		    t*(1 - t)*AsecondDer(copula,t)/A(copula,t)	
	}   		
	   t1 <- integrate(f,subdivisions=100, lower = 0.0, 
                        upper = 1.0, copula = copula)$value      
	   return(t1+(1-t1)*tauphi)	
})

setMethod("Kendalls.tau","frank.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","bb2.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","bb3.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","bb6.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","bb7.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","joe.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        Phi(copula, t)/PhiDer(copula, t)
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    1 + 4 * t1
} )

setMethod("Kendalls.tau","normal.mix.copula", 
function(copula, tol = 1e-2) 
{          
		   assign("tempCop", copula, frame = 1)		
		   Cv <- function(x,y) {
			        x1 <- pcopula(tempCop,x,y)
              	 x2 <- dcopula(tempCop,x,y)
                  return(x1*x2)  
	      }
	  fun <- function(t, integrand, tol)
           { 
           	 val <- unlist(lapply(t, function(t, integrand, tol)
                 { integrate(f = integrand,
                    lower = 0, upper = 1, rel.tol = tol,
                    y = t)$value }
                , integrand, tol = tol))
               val
            }      
		   t1 <- integrate(f = fun,  subdivisions=100, lower = 0, rel.tol = tol,
                        upper = 1, integrand = Cv, tol = tol )$value
          val2 <- 4*t1 - 1
	    val2	
})


################################ 
### Functions Spearmans.rho  ### 
################################

setMethod("Spearmans.rho","gumbel.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        1/(A(copula, t) + 1)^2
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","galambos.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        1/(A(copula, t) + 1)^2
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","husler.reiss.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        1/(A(copula, t) + 1)^2
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","tawn.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        1/(A(copula, t) + 1)^2
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb5.copula",
     function(copula, tol = 1e-5) {
    f <- function(t, copula)
    {
        1/(A(copula, t) + 1)^2
    }
    t1 <- integrate(f, subdivisions = 100, lower = 0., upper = 1., copula = copula)$value
    12 * t1 - 3
} )


setMethod("Spearmans.rho","normal.copula", function(copula, tol = 1e-5) {
	  delta <- as.vector(copula$parameters[1])
     return(6*asin(delta/2.0)/pi)	
})

setMethod("Spearmans.rho","normal.mix.copula", function(copula, tol = 1e-5) 
{    
		   assign("tempCop", copula, frame = 1)
		
	 Cv <- function(x,y) {
			        x1 <- pcopula(tempCop,x,y)
			        return(x1)  
	      }	 
	 fun <- function(t, integrand)
           {  val <- unlist(lapply(t, function(t, integrand)
                 { integrate(f = integrand,
                    lower = 0, upper = 1, 
                    y = t)$value }
                , integrand))
               val
            } 
       
	 t1 <- integrate(f = fun,  subdivisions=100, lower = 0, upper = 1, integrand = Cv )$value
        12*t1 - 3	
})

setMethod("Spearmans.rho","bb4.copula", function(copula, tol = 1e-5) 
{    
		   assign("tempCop", copula, frame = 1)
		
	 Cv <- function(x,y) {
			        x1 <- pcopula(tempCop,x,y)
			        return(x1)  
	      }
 
		 
	 fun <- function(t, integrand)
           {  val <- unlist(lapply(t, function(t, integrand)
                 { integrate(f = integrand,
                    lower = 0, upper = 1, 
                    y = t)$value }
                , integrand))
               val
            } 
	 t1 <- integrate(f = fun,  subdivisions=100, lower = 0, upper = 1, integrand = Cv )$value
          12*t1 - 3	
})


setMethod("Spearmans.rho","frank.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","kimeldorf.sampson.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb1.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb2.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb3.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb6.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","bb7.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )

setMethod("Spearmans.rho","joe.copula",
     function(copula, tol = 1e-5) {
    assign("tempCop", copula, frame = 1)
    Cv <- function(x, y)
    {
        x1 <- pcopula(tempCop, x, y)
        return(x1)
    }
    fun <- function(t, integrand)
    {
        val <- unlist(lapply(t, function(t, integrand)
        {
            integrate(f = integrand, lower = 0, upper = 1, y = t)$value
        }
        , integrand))
        val
    }
    t1 <- integrate(f = fun, subdivisions = 100, lower = 0, upper = 1, integrand = Cv)$value
    12 * t1 - 3
} )



############################# 
### Functions tail.index  ### 
#############################

setMethod("tail.index","gumbel.copula",
     function(copula, ...) {
    LI <- 0
    UI <- 2 - 2 * A(copula, 0.5)
    val <- c(LI, UI)
    names(val) <- c("lower.tail", "upper.tail")
    val
} )

setMethod("tail.index","galambos.copula",
     function(copula, ...) {
    LI <- 0
    UI <- 2 - 2 * A(copula, 0.5)
    val <- c(LI, UI)
    names(val) <- c("lower.tail", "upper.tail")
    val
} )

setMethod("tail.index","husler.reiss.copula",
     function(copula, ...) {
    LI <- 0
    UI <- 2 - 2 * A(copula, 0.5)
    val <- c(LI, UI)
    names(val) <- c("lower.tail", "upper.tail")
    val
} )

setMethod("tail.index","tawn.copula",
     function(copula, ...) {
    LI <- 0
    UI <- 2 - 2 * A(copula, 0.5)
    val <- c(LI, UI)
    names(val) <- c("lower.tail", "upper.tail")
    val
} )

setMethod("tail.index","bb5.copula",
     function(copula, ...) {
    LI <- 0
    UI <- 2 - 2 * A(copula, 0.5)
    val <- c(LI, UI)
    names(val) <- c("lower.tail", "upper.tail")
    val
} )


############################# 
### Method Functions dcdx ###
#############################


setMethod("dcdx","gumbel.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx + ly
    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0 * lxy) * (ad0 * lxy + ad1 * ( - lx + lxy)))/(x * lxy)
    val[u == 0 | u == 1] <- 1
    val[v == 0 & (u >= 0) & (u <= 1)] <- 1
    val[v == 1 & (u >= 0) & (u <= 1)] <- 0
    val
} )

setMethod("dcdx","galambos.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx + ly
    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0 * lxy) * (ad0 * lxy + ad1 * ( - lx + lxy)))/(x * lxy)
    val[u == 0 | u == 1] <- 1
    val[v == 0 & (u >= 0) & (u <= 1)] <- 1
    val[v == 1 & (u >= 0) & (u <= 1)] <- 0
    val
} )

setMethod("dcdx","husler.reiss.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx + ly
    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0 * lxy) * (ad0 * lxy + ad1 * ( - lx + lxy)))/(x * lxy)
    val[u == 0 | u == 1] <- 1
    val[v == 0 & (u >= 0) & (u <= 1)] <- 1
    val[v == 1 & (u >= 0) & (u <= 1)] <- 0
    val
} )

setMethod("dcdx","tawn.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx + ly
    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0 * lxy) * (ad0 * lxy + ad1 * ( - lx + lxy)))/(x * lxy)
    val[u == 0 | u == 1] <- 1
    val[v == 0 & (u >= 0) & (u <= 1)] <- 1
    val[v == 1 & (u >= 0) & (u <= 1)] <- 0
    val
} )

setMethod("dcdx","bb5.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx + ly
    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0 * lxy) * (ad0 * lxy + ad1 * ( - lx + lxy)))/(x * lxy)
    val[u == 0 | u == 1] <- 1
    val[v == 0 & (u >= 0) & (u <= 1)] <- 1
    val[v == 1 & (u >= 0) & (u <= 1)] <- 0
    val
} )


setMethod("dcdx","frank.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","kimeldorf.sampson.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","bb1.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","bb2.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","bb3.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","bb6.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","bb7.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

setMethod("dcdx","joe.copula",
     function(copula,u,v) {
    val <- vector(length = length(u))
    good <- u > 0 & u < 1
    x <- u[good]
    y <- v[good]
    iphip <- InvPhifirstDer(copula, Phi(copula, x) + Phi(copula, y))
    phip <- PhiDer(copula, x)
    val[good] <- iphip * phip
    val[u == 0] <- 1
    val[u == 1] <- v
    val
} )

dcdx.ev.copula <- function(copula,u,v) {
   val <- vector (length = length(u))
    good <- u > 0 & u < 1 & v > 0 & v < 1
    x <- u[good]
    y <- v[good]
    lx <- log(x)
    ly <- log(y)
    lxy <- lx +ly

    t <- lx/lxy
    ad0 <- A(copula, t)
    ad1 <- AfirstDer(copula, t)
    val[good] <- 1 - (exp(ad0*lxy)*(ad0*lxy + ad1*(-lx + lxy)))/(x*lxy)
 
   val[u == 0 | u == 1] <- 1    
   val[v == 0 & (u >= 0) & (u <= 1)] <- 1   
   val[v == 1 & (u >= 0) & (u <= 1)] <- 0        
   val
}

dcdx.archm.copula <- function(copula,u,v) {
  val <- vector (length = length(u))
    good <- u > 0 & u < 1

    x <- u[good]
    y <- v[good]
      
    iphip <- InvPhifirstDer(copula, Phi(copula,x)+Phi(copula,y))
    phip <- PhiDer(copula,x)
    val[good] <- iphip * phip   
  
   val[u == 0] <- 1
   val[u == 1] <- v     
   val
}


######################################
### Plotting Functions for Copulas ###
######################################


setMethod("contour.pcopula","copula",
	function(cop, n=50, ... )
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
    		zmat <- pcopula(cop,  xmat,  ymat)
    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    		contour(val)
    		title("Contour Plot of the Copula", xlab="u",ylab="v")
    		invisible(val)
	})


setMethod("contour.pcopula","empirical.copula",
	function(cop, n=50, ... )
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
    		zmat <- pcopula(cop,  xmat,  ymat)
    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    		contour(val)
    		title("Contour Plot of the Copula", xlab="u",ylab="v")
    		invisible(val)
	})

setMethod("contour.dcopula","copula",
function(cop, n=50, ... ) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
		zmat <- dcopula(cop,  xmat,  ymat) 
    		val <- list(x=divis, y=divis, z=matrix(ncol=n, nrow=n, byrow=FALSE, data=zmat)   )
    		contour(val)
    		title("Contour Plot of the Copula Density", xlab="u",ylab="v")
    		invisible(val)
	})

setMethod("contour.dcopula","empirical.copula",
function(cop, n=50, ... ) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
		bigobj <- kdest(cop@x,cop@y,n=n,lims=c(range(divis),range(divis)))
            	zmat <-   bigobj$z
    		val <- list(x=divis, y=divis, z=matrix(ncol=n, nrow=n, byrow=FALSE, data=zmat)   )
    		contour(val)
    		title("Contour Plot of the Copula Density", xlab="u",ylab="v")
    		invisible(val)
	})


setMethod("persp.dcopula","empirical.copula",
  function(cop, n = 50, ... ) 
  {
     	divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    	divis <- divis[2:(n+1)]
      	xmat <- rep(divis, n )
    	ymat <- rep(divis, each = n)
	bigobj <- kdest(c(cop@x),c(cop@y),N= n, lims = c(range(divis),range(divis)))
        zmat <-   bigobj$z
	val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    	persp3d(divis, divis, zmat, aspect=c(1, 1, 0.5), col = "lightblue", box=FALSE, xlab = "x", ylab = "y", zlab = "")
    	invisible(val)
})

setMethod("persp.pcopula","empirical.copula",
  function(cop, n = 50, ... ) 
{
divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    divis <- divis[2:(n+1)]
    xmat <- rep(divis, n )
    ymat <- rep(divis, each = n)
    zmat <- pcopula(cop,  xmat,  ymat)   
    val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    persp3d(divis, divis, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE, xlab = "x", ylab = "y", zlab = "")
    invisible(val)
})

setMethod("persp.pcopula","copula",
function(cop, n = 50, ... ) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
        		xmat <- rep(divis, n )
    		ymat <- rep(divis, each = n)
    		zmat <- pcopula(cop,  xmat,  ymat)   
    		val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
    		persp3d(divis, divis, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE, xlab = "x", ylab = "y", zlab = "")
    		invisible(val)
	})

setMethod("persp.dcopula","copula",
function(cop, n = 50, ... ) 
{
	divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
	divis <- divis[2:(n+1)]
	xmat <- rep(divis, n )
	ymat <- rep(divis, each = n)
	zmat <- dcopula(cop,  xmat,  ymat)   
	val <- list(x = divis, y = divis, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)   )
	persp3d(divis, divis, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE, xlab = "x", ylab = "y", zlab = "")
	invisible(val)
})



######################### 
### Functions lambda  ### 
#########################

setMethod("lambda","frank.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","kimeldorf.sampson.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","bb1.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","bb2.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","bb3.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","bb6.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","bb7.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )

setMethod("lambda","joe.copula",
     function(copula, t) {
    Phi(copula, t)/PhiDer(copula, t)
} )


############################################# 
### Functions for Bivariate Distributions ###
#############################################


###########################
### Defining class bivd ###
###########################

setClass("bivd", representation(copula="copula", Xmarg = "character", Ymarg = "character",   
     param.Xmarg = "numeric", param.Ymarg = "numeric"))


if(isGeneric("pbivd")) removeGeneric("pbivd")
if(isGeneric("dbivd")) removeGeneric("dbivd")
if(isGeneric("rbivd")) removeGeneric("rbivd")
if(isGeneric("contour.pbivd")) removeGeneric("contour.pbivd")
if(isGeneric("contour.dbivd")) removeGeneric("contour.dbivd")
if(isGeneric("persp.pbivd")) removeGeneric("persp.pbivd")
if(isGeneric("persp.dbivd")) removeGeneric("persp.dbivd")


setClass("empirical.bivd", representation(x = "numeric", y = "numeric"))
setClass("gpd.bivd", representation(copula="copula", Xmarg = "gpd", Ymarg = "gpd"))   


bivd <- function(cop, marginX = "unif", marginY = marginX, 
            param.marginX = c(0,1), param.marginY = param.marginX) 
{
     yyy <- class(cop)
     cl <- substring(yyy,1,nchar(yyy) - 7)
 
#     val <- new( paste(cl,"bivd",sep = "."), copula=cop, Xmarg=marginX, Ymarg=marginY,
     	val <- new( "bivd", copula=cop, Xmarg=marginX, Ymarg=marginY,
         param.Xmarg = param.marginX, param.Ymarg = param.marginY)
       
     val  
}

###############################
### Bivariate Distributions ###
###############################
    
#"gpdjoint.1p" <- function(x, y, copula, x.est, y.est)
#{
#    h <- gpd.1p(x, x.est)
#    g <- gpd.1p(y, y.est)
#    val <- pcopula(copula, h, g)
#    val
#}

#"gpdjoint.2p" <- function(x, y, copula, x.est, y.est)
#{
#    h <- gpd.2p(x, x.est)
#    g <- gpd.2p(y, y.est)
#    val <- pcopula(copula, h, g)
#    val
#}

setGeneric("pbivd", function(dist, x, y) standardGeneric("pbivd"))
setMethod("pbivd","bivd",
	function(dist, x, y) 
	{
   		xmar <- evalFunc(x,get(paste("p",dist@Xmarg,sep = "")),dist@param.Xmarg)
   		ymar <- evalFunc(y,get(paste("p",dist@Ymarg,sep = "")),dist@param.Ymarg )
   		return(pcopula(dist@copula,xmar,ymar))
	})
setMethod("pbivd","gpd.bivd",
	function(dist, x, y) 
	{
   		xmarg <- pgpd(x,dist@Xmarg)
   		ymarg <- pgpd(y,dist@Ymarg )
   		return(pcopula(dist@copula,xmarg,ymarg))
	})

setGeneric("dbivd", function(dist, x, y) standardGeneric("dbivd"))
setMethod("dbivd","bivd",
	function(dist, x, y) 
	{
       		xmar <- evalFunc(x,get(paste("p",dist@Xmarg,sep = "")),dist@param.Xmarg)
       		ymar <- evalFunc(y,get(paste("p",dist@Ymarg,sep = "")),dist@param.Ymarg )
       		dmarx <- evalFunc(x,get(paste("d",dist@Xmarg,sep = "")),dist@param.Xmarg)
       		dmary <- evalFunc(y,get(paste("d",dist@Ymarg,sep = "")),dist@param.Ymarg)
       		return(dcopula(dist@copula,xmar,ymar)*dmarx *dmary)
	})
setMethod("dbivd","gpd.bivd",
	function(dist, x, y) 
	{
   		xmarg <- pgpd(x,dist@Xmarg)
   		ymarg <- pgpd(y,dist@Ymarg )
       		dmargx <- dgpd(x,dist@Xmarg)
       		dmargy <- dgpd(y,dist@Ymarg)
       		return(dcopula(dist@copula,xmarg,ymarg)*dmargx *dmargy)
	})

setGeneric("rbivd", function(dist, n) standardGeneric("rbivd"))
setMethod("rbivd","bivd",
	function(dist,n) 
	{
   		uv <- rcopula(dist@copula,n)
   		Xsim <- evalFunc(uv$x, get(paste("q",dist@Xmarg,sep = "")),dist@param.Xmarg)
   		Ysim <- evalFunc(uv$y, get(paste("q",dist@Ymarg,sep = "")),dist@param.Ymarg)
   		val <- list(x = Xsim, y = Ysim)
   		data.frame(val)
	})
setMethod("rbivd","gpd.bivd",
	function(dist,n) 
	{
   		uv <- rcopula(dist@copula,n)
   		Xsim <- qgpd(dist@Xmarg,uv$x)
   		Ysim <- qgpd(dist@Ymarg,uv$y)
   		val <- list(x = Xsim, y = Ysim)
   		data.frame(val)
	})

#########################
### Graphic functions ###
#########################

setGeneric("contour.dbivd", function(dist, n = 50, xlim = NA, ylim = NA, ...) standardGeneric("contour.dbivd"))
setMethod("contour.dbivd","bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
    		if ( (is.na(xlim[1])) | (length(xlim) != 2)) 
			{ x <- evalFunc(divis,get(paste("q",dist@Xmarg,sep = "")),dist@param.Xmarg)} 
		else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) 
			{y <- evalFunc(divis,get(paste("q",dist@Ymarg,sep = "")),dist@param.Ymarg)} 
		else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- dbivd(dist, xmat, ymat)
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
    		contour(val)
    		title("Contour Plot of the Density", xlab="x",ylab="y")
    		invisible(val)
	})
setMethod("contour.dbivd","gpd.bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
    		if ( (is.na(xlim[1])) | (length(xlim) != 2)) 
			{ x <- qgpd(dist@Xmarg,divis)} 
		else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) 
			{y <- qgpd(dist@Ymarg,divis)} 
		else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- dbivd(dist, xmat, ymat)
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
    		contour(val)
    		title("Contour Plot of the Density", xlab="x",ylab="y")
    		invisible(val)
	})

setGeneric("contour.pbivd", function(dist, n = 50, xlim = NA, ylim = NA, ...) standardGeneric("contour.pbivd"))
setMethod("contour.pbivd","bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
    		if ( (is.na(xlim[1])) | (length(xlim) != 2)) 
			{x <- evalFunc(divis,get(paste("q",dist@Xmarg,sep="")),dist@param.Xmarg)} 
		else {x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) 
			{y <- evalFunc(divis,get(paste("q",dist@Ymarg,sep = "")),dist@param.Ymarg)} 
		else {y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- pbivd(dist, xmat, ymat)
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
    		contour(val)
    		title("Contour Plot of the CDF", xlab="x",ylab="y")
    		invisible(val)
	})
setMethod("contour.pbivd","gpd.bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
    		if ( (is.na(xlim[1])) | (length(xlim) != 2)) 
			{x <- qgpd(dist@Xmarg,divis)} 
		else {x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) 
			{y <- qgpd(dist@Ymarg,divis)} 
		else {y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- pbivd(dist, xmat, ymat)
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
    		contour(val)
    		title("Contour Plot of the CDF", xlab="x",ylab="y")
    		invisible(val)
	})


setGeneric("persp.dbivd", function(dist, n = 50, xlim = NA, ylim = NA, ...) standardGeneric("persp.dbivd"))
setMethod("persp.dbivd","bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
   		if ( (is.na(xlim[1])) | (length(xlim) != 2)) {  
      			x <- evalFunc(divis,get(paste("q",dist@Xmarg,sep = "")),dist@param.Xmarg)} 
		else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) {
     			y <- evalFunc(divis,get(paste("q",dist@Ymarg,sep = "")),dist@param.Ymarg)} 
		else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- dbivd(dist, xmat, ymat)
    		persp3d(x, y, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE,xlab = "x", ylab = "y", zlab = "")
    
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)) 
    		invisible(val)               
	})
setMethod("persp.dbivd","gpd.bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
   		if ( (is.na(xlim[1])) | (length(xlim) != 2)) 
			{x <- qgpd(dist@Xmarg,divis)} 
		else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2) 
			{y <- qgpd(dist@Ymarg,divis)} 
		else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- dbivd(dist, xmat, ymat)
    		persp3d(x, y, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE,xlab = "x", ylab = "y", zlab = "")
    
    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat)) 
    		invisible(val)               
	})

setGeneric("persp.pbivd", function(dist, n = 50, xlim = NA, ylim = NA, ...) standardGeneric("persp.pbivd"))
setMethod("persp.pbivd","bivd",
function(dist, n = 50, xlim = NA, ylim = NA, ...) 
{
	divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
	divis <- divis[2:(n+1)]
    
	if ( (is.na(xlim[1])) | (length(xlim) != 2))   
	{x <- evalFunc(divis,get(paste("q",dist@Xmarg,sep = "")),dist@param.Xmarg)} 
	else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }
	
	if (is.na(ylim[1]) | length(ylim) != 2) {
		y <- evalFunc(divis,get(paste("q",dist@Ymarg,sep = "")),dist@param.Ymarg)} 
	else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
	
	xmat <- rep(x, n )
	ymat <- rep(y, each = n )
	zmat <- pbivd(dist, xmat, ymat)
	persp3d(x, y, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE,xlab = "x", ylab = "y", zlab = "")
	
	val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
	invisible(val)               
})
setMethod("persp.pbivd","gpd.bivd",
	function(dist, n = 50, xlim = NA, ylim = NA, ...) 
	{
    		divis <- seq(from = 0.001, to = 0.999, length = (n+2)) 
    		divis <- divis[2:(n+1)]
    
    		if ( (is.na(xlim[1])) | (length(xlim) != 2))  
			{x <- qgpd(dist@Xmarg,divis)} 
		else {  x <- seq(from = xlim[1], to = xlim[2], length = n) }

    		if (is.na(ylim[1]) | length(ylim) != 2)
			{y <- qgpd(dist@Ymarg,divis)} 
		else { y <- seq(from = ylim[1], to = ylim[2], length = n) }
 
    		xmat <- rep(x, n )
    		ymat <- rep(y, each = n )
    		zmat <- pbivd(dist, xmat, ymat)
		persp3d(x, y, zmat, aspect=c(1, 1, 0.5), col = "lightblue",box=FALSE,xlab = "x", ylab = "y", zlab = "")

    		val <- list(x = x, y = y, z = matrix(ncol = n, nrow = n, byrow = FALSE, data = zmat   ))
    		invisible(val)               
	})



#################################
### VaR Examples Computations ###
#################################

fbar.exp.portf <- function(a, copula, x.est, y.est, lambda1, lambda2, q = NA,
       subdivisions = 1000, tol.f = 5e-5) 
{
    a <- -a
    BE <- exp(a)/lambda2
    UL <- gpd.2p(a - log(lambda1) , x.est)
    ftem <- function(x,BE,copula, x.est, y.est, lambda1, lambda2)
    {
        val <- vector(length = length(x))
        u <- gpd.2q(x,x.est)
        uuu <- BE - lambda1/lambda2*exp(u)
        val[uuu <= 0] <- 1
        val[uuu >= exp(1)] <- 1
        v <- log(uuu[uuu>0 & uuu <exp(1)])
        px <- gpd.2p(v , y.est)
       val[uuu>0 & uuu <exp(1) ] <- dcdx(copula,x[uuu>0 & uuu <exp(1)] ,px)
     val
    }

   ret <- integrate(ftem, 0, UL, abs.tol= tol.f, subdivisions = subdivisions, 
BE=BE,copula=copula,x.est=x.est,y.est=y.est,lambda1=lambda1,lambda2=lambda2)$value
   if (!is.na(q)) ret <- ret - q
   ret
}

fbar.lin.portf <- function(Q, copula, x.est, y.est, lambda1, lambda2, 
        subdivisions = 1000, tol.f = 5e-5)  
{   
    ftem <- function(x,Q,copula, x.est, y.est, lambda1, lambda2)
    {
        qq <- (Q - gpd.1q(x,y.est)*lambda2)/lambda1
        px <- vector(length = length(qq))
        px[is.infinite(qq) & qq < -1] <- 0
        px[is.infinite(qq) & qq > 1] <- 1
        px[!is.infinite(qq) ] <- gpd.1p(qq[!is.infinite(qq) ] , x.est)
        val <- vector(length = length(x))
        val[ px >=1] <- 1
        val[ px <= 0] <- 1
        val[px > 0 & px < 1] <- dcdx.max.copula(copula,x[px > 0 & px < 1] ,px[px > 0 & px < 1] )
        val
    }
   integrate(ftem, 0, 1,Q = Q,copula = copula, x.est = x.est, y.est = y.est,
       lambda1 = lambda1, lambda2 = lambda2, subdivisions = 1000, abs.tol = 5e-5)$value
}

VaR.exp.portf <- function(Q, copula, x.est, y.est, lambda1, lambda2,
      range, subdivisions = 1000, tol.f = 5e-5)   
{     
     uniroot(fbar.exp.portf, interval=range, copula = copula,x.est = x.est, y.est = y.est,
       lambda1 = lambda1, lambda2 = lambda2, q = Q,
          subdivisions = subdivisions,tol.f=tol.f)$root
}

VaR.sim <- function(n, p, copula, x.est, y.est, lambda1, lambda2 ) 
{
    sim.uv <- rcopula(copula, n)
    Xsim <- qgpd(x.est,sim.uv$x)
    Ysim <- qgpd(y.est,sim.uv$y)
    ww <- log(lambda1*exp(Xsim) + lambda2*exp(Ysim))
    varv <- -as.vector(quantile(ww,p))
    esf <- vector(length = length(p), mode = "numeric")
    for (i in 1:length(p))  esf[i] <- -mean(ww[ww < -varv[i]])

    val <- c(n,varv,esf)
    nn <-  vector(length = length(p)*2 + 1, mode = "character")
    nn[1] <- "Simulation size"
    nn[2:(length(p)+1)] <- paste("VaR_p = ",p,sep = "")
    nn[(length(p)+2):(2*length(p)+1)] <- paste("ES_p = ",p,sep = "")
    names(val) <- nn
    val
}
    
    
#####################
### Kendall Table ###
#####################

"Kend.table" <- 
structure(.Data = list("tau" = c(0., 0.10000000000000001, 0.20000000000000001, 0.29999999999999999, 0.40000000000000002,
    0.5, 0.59999999999999998, 0.69999999999999996, 0.80000000000000004, 
    0.90000000000000002, 1.)
, "B1" = c(0., 0.1564344650402309, 0.3090169943749474, 0.45399049973954669, 0.58778525229247314,
    0.70710678118654746, 0.80901699437494745, 0.89100652418836779, 0.95105651629515364,
    0.98768834059513777, 1.)
, "B3" = c(0., 0.91000000000000003, 1.8600000000000001, 2.9199999999999999, 4.1600000000000001,
    5.7400000000000002, 7.9299999999999997, 11.4, 18.199999999999999, 
    20.899999999999999, Inf)
, "B4" = c(0., 0.22, 0.5, 0.85999999999999999, 1.3300000000000001, 2., 3., 4.6699999999999999, 8.,
    18., Inf)
, "B6" = c(1., 1.1111111111111109, 1.25, 1.4285714285714279, 1.666666666666667, 2., 2.5, 
    3.333333333333333, 5.0000000000000009, 10., Inf)
, "B7" = c(0., 0.34000000000000002, 0.51000000000000001, 0.69999999999999996, 0.94999999999999996,
    1.28, 1.79, 2.6200000000000001, 4.29, 9.3000000000000007, Inf)
, "B8" = c(0., 0.66000000000000003, 0.87, 1.1100000000000001, 1.4099999999999999, 1.8, 
    2.3900000000000001, 3.3399999999999999, 5.2400000000000002, 10.9, Inf)
, "BB1.delta" = c(1., 1.1000000000000001, 1.2, 1.3, 1.3999999999999999, 1.6499999999999999, 2., 
    2.6499999999999999, 3., 5., Inf)
, "bb1.theta" = c(0., 0.020202020202019891, 0.083333333333333481, 0.1978021978021976, 0.38095238095238099,
    0.42424242424242431, 0.5, 0.51572327044025146, 1.3333333333333339, 
    2.0000000000000009, Inf)
, "Twan.copula" = c(0., 1.1258349999999999, 1.290394, 1.5154939999999999, 1.843656, 2.3701919999999999, 
    3.3637899999999998, 5.997922, 37.203530000000001, 172.10230000000001, Inf)
, "TC.a.b" = c(0.90000000000000002, 0.90000000000000002, 0.90000000000000002, 0.90000000000000002, 
    0.90000000000000002, 0.90000000000000002, 0.90000000000000002, 0.90000000000000002,
    0.90000000000000002, 0.94999999999999996, Inf)
, "BB2.theta" = c(0., 0.14999999999999999, 0.20000000000000001, 0.29999999999999999, 0.40000000000000002,
    0.5, 0.59999999999999998, 0.69999999999999996, 0.78000000000000003, 
    0.81000000000000005, Inf)
, "BB2.delta" = c(0., 0.41547420000000002, 1.2451000000000001, 1.4251780000000001, 1.6769670000000001,
    2.0480670000000001, 2.6324589999999999, 3.6496490000000001, 5.9888640000000004,
    8.5999999999999996, Inf)
, "BB3.theta" = c(1., 1.05, 1.1000000000000001, 1.2, 1.3, 1.3999999999999999, 1.5, 1.8, 2.1000000000000001,
    2.2010000000000001, Inf)
, "BB3.delta" = c(0., 0.1152268, 0.26838640000000002, 0.37031439999999999, 0.54731560000000001, 
    0.84593099999999999, 1.37849, 1.8886909999999999, 3.7887529999999998, 
    13.029999999999999, Inf)
, "BB5.theta" = c(1., 1.05, 1.1000000000000001, 1.2, 1.3, 1.3999999999999999, 1.5, 1.8, 1.8999999999999999,
    1.95, Inf)
, "BB5.delta" = c(0., 0.26372449999999997, 0.37258530000000001, 0.43810700000000002, 0.54184339999999998,
    0.69927689999999998, 0.94607929999999996, 1.1348119999999999, 1.9201589999999999,
    4.4188349999999996, Inf)
, "BB6.theta" = c(1., 1.05, 1.1000000000000001, 1.1499999999999999, 1.2, 1.25, 1.3500000000000001, 1.55,
    1.8, 2., Inf)
, "BB6.delta" = c(1., 1.079931, 1.18201, 1.31548, 1.495727, 1.750732, 2.0865680000000002, 
    2.5480239999999998, 3.4635950000000002, 6.4493419999999997, Inf)
, "BB7.theta" = c(1., 1.05, 1.1000000000000001, 1.25, 1.5, 1.8, 2.2999999999999998, 2.7999999999999998,
    3.5, 9.5, Inf)
, "BB7.delta" = c(0., 0.1649669, 0.3882816, 0.58719399999999999, 0.81623789999999996, 1.2346159999999999,
    1.9072309999999999, 3.8238219999999998, 10.30184, 60.21819, Inf)
, "BB4.theta" = c(0., 0.19, 0.29999999999999999, 0.40000000000000002, 0.5, 0.69999999999999996, 
    0.90000000000000002, 1.3, 1.8999999999999999, 2.2999999999999998, Inf)
, "BB4.delta" = c(0., 0.1750353, 0.30702400000000002, 0.43810700000000002, 0.59777539999999996, 
    0.75476860000000001, 1.004845, 1.3052299999999999, 1.8524259999999999, 
    3.9417589999999998, Inf)
)
, names = c("tau", "B1", "B3", "B4", "B6", "B7", "B8", "BB1.delta", "bb1.theta", "Twan.copula", 
    "TC.a.b", "BB2.theta", "BB2.delta", "BB3.theta", "BB3.delta", "BB5.theta", 
    "BB5.delta", "BB6.theta", "BB6.delta", "BB7.theta", "BB7.delta", "BB4.theta", 
    "BB4.delta")
, row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
, .guiColInfo = list(c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-2")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-3")
, c("TXPROP_ColWidth<-12", "TXPROP_ColJustification<-Right", "TXPROP_ObjectNote<-", 
    "TXPROP_ColFloatFormat<-Decimal", "TXPROP_ColPrecision<-5")
)
, class = "data.frame"
)


##################
# Stand alone (default) empirical.copula function

empirical.copula <- function(x, y = NA) 
{
    data <- NA
    if(inherits(x,"empirical.copula")) {data <- x}
    else if ( is.list(x))
    {
       temp <- data.frame(x)
       data <- new("empirical.copula",x = temp[,1], y = temp[,2])
    }
    else if (!is.na(y[1]) & length(x) == length(y))
    {
       data <- new("empirical.copula",x = x, y = y)
    }
    else
    {
      stop("Can't create empirical copula. Unknown parameters supplied")
    }

    inputsOK <- TRUE
    message <- ""
    if (sum(data@x < 0) > 0 | sum(data@y < 0) > 0) 
    {
        inputsOK <- FALSE
        message <-
           paste(message," At least one of the observed data points has a negative coordinate", sep = " \n")
    }
    if (sum(data@x > 1) > 0 | sum(data@y > 1) > 0) 
    {
        inputsOK <- FALSE
        message <-
           paste(message," At least one of the observed data points has a coordinate greater than 1",
                  sep = " \n")
    }
    if (!inputsOK) {stop(message)}         
    data    
}
 

#########################################
# Multivariate Gaussian distribution

lm.diag <- function (LM) 
{
    list(hat = hatvalues(LM), stdres = rstandard(LM), studres = rstudent(LM), 
        cooks = cooks.distance(LM), covratio = covratio(LM))
}

p2Dnorm <-
function (x, mean = rep(0, 2), std = rep(1, 2), rho = 0) 
{
    if (!missing(mean)) {
        if (!is.vector(mean) || length(mean) != 2) 
            stop("mean must be a 2D-vector")
    }
    if (is.vector(x)) 
        x <- matrix(x, ncol = length(x))
    if (!missing(rho)) {
        if (length(rho) != 1) 
            stop("rho must be a number (i.e. have length 1)")
        if (any(rho < -1 | rho > 1)) 
            stop("rho must be between -1 and +1")
    }
    if (!missing(std)) {
        if (any(std <= 0)) 
            stop("Negative or zero std found")
        if (!is.vector(std) || length(std) != 2) 
            stop("std must be a 2D-vector")
    }
    ff <- function(z) {
        m21 <- mean[2] + rho * std[2] * (z - mean[1])/std[1]
        s21 <- std[2] * sqrt(1 - rho * rho)
        pnorm(x[2], mean = m21, sd = s21) * dnorm(z, mean = mean[1], 
            sd = std[1])
    }
    integrate(ff, lower = -Inf, upper = x[1])$value
}

pmvnorm <- function (x, mean = rep(0, d), cov = diag(d), sd, rho, d = 2, 
    sigma = cov, log = FALSE) 
{
    if (!missing(cov) && (d.cov != d.mean)) 
        stop("The dimensions of cov and mean do not match")
    if (missing(mean)) 
        mean <- rep(0, d)
    if (!is.vector(mean)) 
        stop("mean must be a vector")
    d.mean <- length(mean)
    if (is.vector(x) | dim(x)[2] != d) 
        stop("Error, x should be a matrix with d columns")
    d.x <- dim(x)[2]
    if ((!missing(cov) || !missing(sigma)) && (!missing(sd) || 
        !missing(rho))) 
        warning("Either (cov or sigma) or (sd and rho) should be given, NOT BOTH. cov will be used")
    if (!missing(cov) && !missing(sigma)) 
        warning("Either cov or sigma should be given, NOT BOTH. cov will be used")
    if (missing(cov) && missing(sigma)) 
        cov <- diag(d)
    if (missing(cov) && !missing(sigma)) 
        cov <- sigma
    if (!missing(cov)) {
        if (!is.matrix(cov) || diff(dim(cov)) || any(abs(cov - 
            t(cov)) > sqrt(.Machine$double.eps))) 
            stop("cov must be a square symmetric matrix")
        d.cov <- nrow(cov)
    }
    else {
        if (missing(sd) || missing(rho)) 
            stop("Both parameters sd and rho must be present if cov is missing")
        if (!missing(rho)) {
            if (length(rho) != 1) 
                stop("rho must be a number (i.e. have length 1)")
            if (any(rho < -1 | rho > 1)) 
                stop("rho must be between -1 and +1")
        }
        if (!missing(sd)) {
            if (any(sd <= 0)) 
                stop("Negative or zero sd found")
            if (!is.vector(sd)) 
                stop("sd must be a vector")
            d.sd <- length(sd)
        }
    }
    if (!missing(cov)) 
        sigma <- cov
    else sigma <- (1 - rho) * diag(sd * sd) + rho * outer(sd, 
        sd, "*")
    lower <- rep(-Inf, dim(sigma)[1])
    OneOverSigma <- sqrt(diag(sigma))
    upper <- (x - mean) * OneOverSigma
    mean <- rep(0, length(lower))
    corr <- sigma * outer(OneOverSigma, OneOverSigma, "*")
    RET <- mvt(lower = lower, upper = upper, df = 0, corr = corr, 
        delta = mean, maxpts = 25000, abseps = 0.001, releps = 0)
    return(RET$value)
}

qmvnorm <-
function (p, interval = NULL, tail = c("lower.tail", "upper.tail", 
    "both.tails"), mean = 0, corr = NULL, sigma = NULL, algorithm = GenzBretz(), 
    ...) 
{
    if (length(p) != 1 || (p <= 0 || p >= 1)) 
        stop(sQuote("p"), " is not a double between zero and one")
    dots <- dots2GenzBretz(...)
    if (!is.null(dots$algorithm) && !is.null(algorithm)) 
        algorithm <- dots$algorithm
    tail <- match.arg(tail)
    if (tail == "both.tails" && p < 0.5) 
        stop("cannot compute two-sided quantile for p < 0.5")
    dim <- length(mean)
    if (is.matrix(corr)) 
        dim <- nrow(corr)
    if (is.matrix(sigma)) 
        dim <- nrow(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    dim <- length(args$mean)
    pfct <- function(q) {
        switch(tail, both.tails = {
            low <- rep(-abs(q), dim)
            upp <- rep(abs(q), dim)
        }, upper.tail = {
            low <- rep(q, dim)
            upp <- rep(Inf, dim)
        }, lower.tail = {
            low <- rep(-Inf, dim)
            upp <- rep(q, dim)
        }, )
        pmvnorm(lower = low, upper = upp, mean = args$mean, corr = args$corr, 
            sigma = args$sigma, algorithm = algorithm) - p
    }
    if (is.null(interval) || length(interval) != 2) {
        if (mean == 0) {
            cr <- args$corr
            if (!is.null(args$sigma)) 
                cr <- cov2cor(args$sigma)
            interval <- approx_interval(p = p, tail = tail, corr = cr, 
                df = 0)
        }
        else {
            interval <- c(ifelse(tail == "both.tails", 0, -10), 
                10)
        }
    }
    if (is.null(dots$uniroot)) {
        qroot <- uniroot(pfct, interval = interval)
    }
    else {
        qroot <- do.call("uniroot", list(pfct, interval = interval, 
            dots$uniroot))
    }
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}

dmvnorm <- function (x, mean, cov, sd, rho, d = 2, sigma = cov, log = FALSE) 
{
    if (!missing(cov) && (d.cov != d.mean)) 
        stop("The dimensions of cov and mean do not match")
    if (missing(mean)) 
        mean <- rep(0, d)
    if (!is.vector(mean)) 
        stop("mean must be a vector")
    d.mean <- length(mean)
    if (is.vector(x) | dim(x)[2] != d) 
        stop("Error, x should be a matrix with d columns")
    d.x <- dim(x)[2]
    if (d  != d.x) 
        stop("Error, the number of columns of x should be equal to d")
    if ((!missing(cov) || !missing(sigma)) && (!missing(sd) || 
        !missing(rho))) 
        warning("Either (cov or sigma) or (sd and rho) should be given, NOT BOTH. cov will be used")
    if (!missing(cov) && !missing(sigma)) 
        warning("Either cov or sigma should be given, NOT BOTH. cov will be used")
    if (missing(cov) && missing(sigma)) 
        cov <- diag(d)
    if (missing(cov) && !missing(sigma)) 
        cov <- sigma
    if (!missing(cov)) {
        if (!is.matrix(cov) || diff(dim(cov)) || any(abs(cov - 
            t(cov)) > sqrt(.Machine$double.eps))) 
            stop("cov must be a square symmetric matrix")
        d.cov <- nrow(cov)
    }
    else {
        if (missing(sd) || missing(rho)) 
            stop("Both parameters sd and rho must be present if cov is missing")
        if (!missing(rho)) {
            if (length(rho) != 1) 
                stop("rho must be a number (i.e. have length 1)")
            if (any(rho < -1 | rho > 1)) 
                stop("rho must be between -1 and +1")
        }
        if (!missing(sd)) {
            if (any(sd <= 0)) 
                stop("Negative or zero sd found")
            if (!is.vector(sd)) 
                stop("sd must be a vector")
            d.sd <- length(sd)
        }
    }
    if (!missing(cov)) 
        sigma <- cov
    else sigma <- (1 - rho) * diag(sd * sd) + rho * outer(sd, 
        sd, "*")
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if (log) 
        return(logretval)
    exp(logretval)
}

rmvnorm <-
function (n, mean = rep(0, d), cov = diag(d), sd, rho, d = 2, 
    sigma = cov, method = c("svd", "chol")) 
{
    if (!missing(mean)) {
        if (!is.vector(mean)) 
            stop("mean must be a vector")
        d.mean <- length(mean)
    }
     if (missing(mean) && !missing(cov)) {
        d.mean <- dim(cov)[1]
        mean <- rep(0, d.mean)
    }
    if (missing(mean) && !missing(sigma)) {
        d.mean <- dim(sigma)[1]
        mean <- rep(0, d.mean)
    }
 if (!missing(cov)) 
	{
		if (!is.array(cov))
			stop("cov should be an array or a matrix")
		if (dim(cov)[1] != dim(cov)[2])
	         stop("cov should be a square matrix")
	      d.cov <- dim(cov)[2]
	 }
    if (!missing(cov) && (d.cov != d.mean)) 
        stop("The dimensions of cov and mean do not match")
    if ((!missing(cov) || !missing(sigma)) && (!missing(sd) || 
        !missing(rho))) 
        warning("Either (cov or sigma) or (sd and rho) should be given, NOT BOTH. cov will be used")
    if (!missing(cov) && !missing(sigma)) 
        warning("Either cov or sigma should be given, NOT BOTH. cov will be used")
    if (missing(cov) && !missing(sigma)) 
        cov <- sigma
    if (!missing(cov)) {
        if (!is.matrix(cov) || diff(dim(cov)) || any(abs(cov - 
            t(cov)) > sqrt(.Machine$double.eps))) 
            stop("cov must be a square symmetric matrix")
        d.cov <- nrow(cov)
    }
    else {
        if (missing(sd) || missing(rho)) 
            stop("Both parameters sd and rho must be present if cov is missing")
        if (!missing(rho)) {
            if (length(rho) != 1) 
                stop("rho must be a number (i.e. have length 1)")
            if (any(rho < -1 | rho > 1)) 
                stop("rho must be between -1 and +1")
        }
        if (!missing(sd)) {
            if (any(sd <= 0)) 
                stop("Negative or zero sd found")
            if (!is.vector(sd)) 
                stop("sd must be a vector")
            d.sd <- length(sd)
        }
    }
     if (!missing(cov)) 
        sigma <- cov
    else sigma <- (1 - rho) * diag(sd * sd) + rho * outer(sd, 
        sd, "*")
    method <- match.arg(method)
    if (method == "svd") {
        ev <- eigen(sigma, sym = TRUE)$values
        if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1]))) 
            warning("sigma is numerically not positive definite")
        sigsvd <- svd(sigma)
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    retval
}

rmvgaussian.copula <- function (N, d=2, rho, Sigma) 
{
	if ( missing(rho) && missing(Sigma)) stop("Error, rho or Sigma should be provided")
	if (!missing(Sigma)) 
	{
		if (!is.array(Sigma))
			stop("Sigma should be an array or a matrix")
		if (dim(Sigma)[1] != dim(Sigma)[2])
	         stop("Sigma should be a square matrix")
		d <- dim(Sigma)[2]
	}
	if (!missing(rho)) 
	{
		if ( (rho<0) || (rho>1)) stop("Error, rho should be between 0 and 1")
		Sigma <- array(rho, dim=c(d,d))
		diag(Sigma) <- rep(1,d)	
	}
	SDtmp <- rmvnorm(N,mean=rep(0,d),cov=Sigma)
	SD <- apply(SDtmp,2,pnorm)
	SD
}