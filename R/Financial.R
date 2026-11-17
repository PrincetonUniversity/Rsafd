bscall <- function (TAU = 0.03, K = 1, S, R = 0.1, SIG = 0.15) 
{
	if (SIG==0) CALL <- max(S-exp(-R*TAU)*K,0)
	else
	{
	    d1 <- log(S/K) + TAU * (R + SIG^2/2)
   	 	d1 <- d1/(SIG * sqrt(TAU))
    	d2 <- d1 - SIG * sqrt(TAU)
    	CALL <- S * pnorm(d1) - K * exp(-R * TAU) * pnorm(d2)
    }
    return(CALL)
}

fns <- function (x, THETA) 
{
    FORWARD <- THETA[1] + (THETA[2] + THETA[3] * x) * exp(-x/THETA[4])
    FORWARD
}

isig <- function (TAU, K, S, R, CALL, TOL = 1e-05, JMAX = 1000) 
{
    SIG1 <- 0
    SIG2 <- 1
	C1 <- max(S-exp(-R*TAU)*K,0)
	C2 <- S
    if ((CALL < C1) | (C2 < CALL) ) 
        return(-1)
    for (J in 1:JMAX) {
        MIDSIG <- (SIG2 + SIG1)/2
        TC <- bscall(TAU, K, S, R, MIDSIG)
        if (TC < CALL) 
            SIG1 <- MIDSIG
        else SIG2 <- MIDSIG
        if (SIG2 - SIG1 < TOL) 
            return(SIG1)
    }
}



bns <- function (COUPON, AI, LIFE, X = 100, THETA = c(0.06, 0, 0, 1)) 
{
    NbBonds <- length(COUPON)
    TT <- THETA[3] * THETA[4]
    TTT <- THETA[4] * (THETA[2] + TT)
    LL <- floor(1 + LIFE)
    PRICE <- rep(0, NbBonds)
    DURATION <- rep(0, NbBonds)
    for (I in 1:NbBonds) {
        x <- seq(to = LIFE[I], by = 1, length = LL[I])
        EX <- exp(-x/THETA[4])
        DISCOUNT <- exp(x * THETA[1] + (TTT * (1 - EX)) - TT * 
            (EX * x))
        CF <- rep((COUPON[I] * X)/100, LL[I])
        CF[LL[I]] <- CF[LL[I]] + X
        PRICE[I] <- sum(CF * DISCOUNT)
        NUM <- sum(x * CF * DISCOUNT)
        DURATION[I] <- NUM/PRICE[I]
    }
    PRICE <- PRICE - AI
    list(price = PRICE, duration = DURATION)
}



yns <- function (x, THETA) 
{
    TT <- THETA[3] * THETA[4]
    TTT <- THETA[4] * (THETA[2] + TT)
    EX <- exp(-x/THETA[4])
    YIELD <- THETA[1] + (TTT * (1 - EX))/x - TT * EX
    YIELD
}
