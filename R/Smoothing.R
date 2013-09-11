kdest <- function (X, Y, H, N = 256, lims = c(range(X), range(Y)), COND = 0, 
    PLOT = TRUE, XNAME = "X", YNAME = "Y") 
{
    LX <- length(X)
    if (length(Y) != LX) 
        stop("Vectors must have same lengths")
    GRIDX <- seq((5 * lims[1] - lims[2])/4, (5 * lims[2] - lims[1])/4, 
        length = N)
    DELTAX <- (3 * (lims[2] - lims[1]))/(2 * 256)
    GRIDY <- seq((5 * lims[3] - lims[4])/4, (5 * lims[4] - lims[3])/4, 
        length = N)
    DELTAY <- (3 * (lims[4] - lims[3]))/(2 * 256)
    if (missing(H)) 
        H <- c(bandwidth.nrd(X), bandwidth.nrd(Y))
    H <- H/4
    DX <- outer(GRIDX, X, "-")/H[1]
    DY <- outer(GRIDY, Y, "-")/H[2]
    Z <- matrix(dnorm(DX), N, LX) %*% t(matrix(dnorm(DY), N, 
        LX))/(LX * H[1] * H[2])
    cat("Bandwidths values used: ", H, "\n")
    MARGX <- outer(DELTAX * DELTAY * apply(Z, 1, sum), rep(1, 
        N), "*")
    if (COND == 1) 
        Z[MARGX > 1e-05] <- Z[MARGX > 1e-05]/MARGX[MARGX > 1e-05]
    MARGY <- outer(rep(1, N), DELTAX * DELTAY * apply(Z, 2, sum), 
        "*")
    if (COND == 2) 
        Z[MARGY > 1e-05] <- Z[MARGY > 1e-05]/MARGY[MARGY > 1e-05]
    if (PLOT == TRUE) {
        persp3d(GRIDX, GRIDY, Z, aspect = c(1, 1, 0.5), col = "lightblue", 
            box="FALSE",xlab = XNAME, ylab = YNAME, zlab = "")
    }
    list(deltax = DELTAX, deltay = DELTAY, gridx = GRIDX, gridy = GRIDY, 
        z = Z)
}

kreg <- function (x, y, xpred = x, kernel = "gaussian", bandwidth) 
{
#cat("is x=xpred?", x== xpred,"\n")    
#cat(" x = ", x,"\n")    
#cat(" xpred = ", xpred,"\n")    

#cat("is.vector(x)", is.vector(x),"\n")    
#cat("is.matrix(x)", is.matrix(x),"\n")    
#cat("is.array(x)", is.array(x),"\n")    

	if (is.vector(x)) {
        n <- length(x)
        p <- 1
        x <- array(x, dim=c(n,p))
    }
    else if (is.matrix(x)) {
        n <- dim(x)[1]
        p <- dim(x)[2]
    }
    else if (is.array(x)) {
        n <- dim(x)[1]
        p <- dim(x)[2]
    }

#cat("n = ",n, "  p = ",p,"\n")    

#cat("is.matrix(x)", is.matrix(x),"\n")    
#cat("is.array(x)", is.array(x),"\n")    

#cat("is.vector(xpred)", is.vector(xpred),"\n")    
#cat("is.matrix(xpred)", is.matrix(xpred),"\n")    
#cat("is.array(xpred)", is.array(xpred),"\n")        
    
    if (is.vector(xpred)) {
        npred <- length(xpred)
        ppred <- 1
        xpred <- array(xpred, dim=c(npred,ppred))
    }
    else if (is.matrix(xpred)) {
        npred <- dim(xpred)[1]
        ppred <- dim(xpred)[2]
    }
    else if (is.array(xpred)) {
        npred <- dim(xpred)[1]
        ppred <- dim(xpred)[2]
    }
    
#cat("is.vector(xpred)", is.vector(xpred),"\n")    
#cat("is.matrix(xpred)", is.matrix(xpred),"\n")    
#cat("is.array(xpred)", is.array(xpred),"\n")        
#cat("dim(xpred) = ",dim(xpred),"\n")    
#cat("npred = ",npred, "  ppred = ",ppred,"\n")    

    if (is.vector(y)) {
        ny <- length(y)
        py <- 1
        y <- as.array(y, dim=c(ny,py))
    }
    else if (is.matrix(y)) {
        ny <- dim(y)[1]
        py <- dim(y)[2]
    }
    if (ny != n) {
        stop("The number of rows of x should be the same\n  as the number of rows of y")
    }
    if (py != 1) {
        stop("y should not have more than one column")
    }
    if (p != ppred) {
        stop("The number of columns of x should be the same\n  as the number of columns of xpred")
    }
    ypred <- rep(0, npred)
   xpp <- array(0,dim = c(npred,n,p))

#cat("dim xpp = ",dim(xpp),"\n")    
   
   for ( k in 1:p) xpp[,,k] <- outer(xpred[,k], x[,k], "-")
   dist2 <- apply(xpp^2,c(1,2),sum)/(bandwidth*bandwidth)

#cat("dim dist2 = ",dim(dist2),"\n")    
#cat("dim xpp = ",dim(xpp),"\n")    

   K <- switch(kernel,
   	gaussian = exp(-dist2),
	box = (dist2 < 1),
	triangle = ((dist2<1)*(1-sqrt(dist2))),
	parzen = (dist2<0.25)*(0.75 - dist2) + (0.25<= dist2)*(dist2<9/4)*((9/8)-1.5*sqrt(dist2)+0.5*dist2)
	)
	
   NUM <- K%*%y
   DEN <- K%*%array(1,dim=c(n,1))
   if (sum(is.na(DEN))>0) stop("Attempt to divide by 0")
    ypred <- NUM/DEN
    return(list(xpred = xpred, ypred = ypred, bandwidth = bandwidth))
}



twoDkreg <- function (X, Y, B, XPRED = X, FITTED = FALSE, N = 256, lims = apply(X, 
    2, range), X1NAME = "X1", X2NAME = "X2", YNAME = "Y", PLOT = TRUE) 
{
    n <- dim(X)[1]
    if (length(Y) != n) 
        stop("Response must have same length as explanatory variables")
    if (missing(B)) 
        B <- c(bandwidth.nrd(X[, 1]), bandwidth.nrd(X[, 2]))/4
    if (length(B) < 2) 
        B <- B * c(bandwidth.nrd(X[, 1]), bandwidth.nrd(X[, 2]))/4
    if (length(B) > 2) 
        stop("Error: when provided, B should be a numeric vector of length 1 or 2", 
            "\n")
    cat("Bandwidths values used: ", B, "\n")
    if (!missing(XPRED)) {
        npred <- dim(XPRED)[1]
        DX1 <- outer(XPRED[, 1], X[, 1], "-")/B[1]
        DX2 <- outer(XPRED[, 2], X[, 2], "-")/B[2]
        NUM <- (outer(rep(1, npred), Y, "*") * matrix(dnorm(DX1), 
            npred, n)) %*% t(matrix(dnorm(DX2), npred, n))
        DEN <- matrix(dnorm(DX1), npred, n) %*% t(matrix(dnorm(DX2), 
            npred, n))
        YPRED <- diag(NUM/DEN)
        OUT <- list(xpred = XPRED, ypred = YPRED)
    }
    else {
        if (FITTED == TRUE) {
            DX1 <- outer(X[, 1], X[, 1], "-")/B[1]
            DX2 <- outer(X[, 2], X[, 2], "-")/B[2]
            NUM <- (outer(rep(1, n), Y, "*") * matrix(dnorm(DX1), 
                n, n)) %*% t(matrix(dnorm(DX2), n, n))
            DEN <- matrix(dnorm(DX1), n, n) %*% t(matrix(dnorm(DX2), 
                n, n))
            Yhat <- diag(NUM/DEN)
            SSR <- sum((Y - Yhat) * (Y - Yhat))
            OUT <- list(fitted = Yhat, ssr = SSR)
        }
        else {
            GRIDX1 <- seq((11 * lims[1, 1] - lims[2, 1])/10, 
                (11 * lims[2, 1] - lims[1, 1])/10, length = N)
            GRIDX2 <- seq((11 * lims[1, 2] - lims[2, 2])/10, 
                (11 * lims[2, 2] - lims[1, 2])/10, length = N)
            DX1 <- outer(GRIDX1, X[, 1], "-")/B[1]
            DX2 <- outer(GRIDX2, X[, 2], "-")/B[2]
            NUM <- (outer(rep(1, N), Y, "*") * matrix(dnorm(DX1), 
                N, n)) %*% t(matrix(dnorm(DX2), N, n))
            DEN <- matrix(dnorm(DX1), N, n) %*% t(matrix(dnorm(DX2), 
                N, n))
            R <- NUM/DEN
            if (PLOT == TRUE) 
                persp3d(GRIDX1, GRIDX2, R, aspect = c(1, 1, 0.5), 
                  col = "lightblue", box=FALSE, xlab = X1NAME, ylab = X2NAME, 
                  zlab = "")
            OUT <- list(gridx1 = GRIDX1, gridx2 = GRIDX2, y = R)
        }
    }
    OUT
}

#########################################
# L1 regression

l1fit <-
function (x, y, intercept = TRUE) 
{
    if (intercept) 
        rq(y ~ x, tau = 0.5)
    else rq(y ~ x - 1, tau = 0.5)
}

