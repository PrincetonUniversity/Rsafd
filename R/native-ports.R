################################################################################
# Pure-R replacements for compiled routines.                                   #
#                                                                              #
# Earlier releases called the C routines "jointDensity" and                    #
# "inverseDerivPhiBB" through .C(), but their C source was never part of the   #
# package, so those calls could not be resolved. The functions below take the  #
# same positional arguments as the original .C() calls and return a list laid #
# out the same way, so the calling code is unchanged apart from the call name. #
################################################################################

# Two-dimensional histogram of the sample (x, y) on an nx-by-ny grid of cells
# with upper-right corners xur, yur and widths dx, dy. A point is counted in the
# cell (i, j) when xur[i] - dx < x <= xur[i] and yur[j] - dy < y <= yur[j].
# Counts are stored row by row with rows indexing y and columns indexing x,
# matching matrix(val, nrow = ny, ncol = nx, byrow = TRUE) in Pearson.test.
.jointDensity <- function(x, y, nsamp, xur, yur, nx, ny, val, dx, dy) {
  ix <- findInterval(x, xur - dx, left.open = TRUE)
  iy <- findInterval(y, yur - dy, left.open = TRUE)
  ix <- pmin(pmax(ix, 1L), nx)
  iy <- pmin(pmax(iy, 1L), ny)
  counts <- tabulate((iy - 1L) * nx + ix, nbins = nx * ny)
  val[] <- as.double(counts)
  list(x, y, nsamp, xur, yur, nx, ny, val, dx, dy)
}

# For each value x[i], solve PhiDer(copula, w) = x[i] for w in (0, 1), where
# PhiDer is the derivative of the Archimedean generator of the BB family
# identified by `family` (1, 2, 3, 6 or 7) with parameters (theta, delta).
# PhiDer is increasing on (0, 1), so the root is found by bisection, stopping
# when the bracket is narrower than epsX or |PhiDer(w) - x| < epsY.
# conv[i] is the number of iterations used, or -1 if maxit was reached first.
.inverseDerivPhiBB <- function(x, n, delta, theta, result, family,
                               epsX, epsY, maxit, conv) {
  copula <- switch(as.character(family),
    "1" = bb1.copula(theta, delta),
    "2" = bb2.copula(theta, delta),
    "3" = bb3.copula(theta, delta),
    "6" = bb6.copula(theta, delta),
    "7" = bb7.copula(theta, delta),
    stop("Unknown BB copula family: ", family)
  )
  for (i in seq_len(n)) {
    lo <- 0
    hi <- 1
    w <- 0.5
    status <- -1L
    for (it in seq_len(maxit)) {
      w <- (lo + hi) / 2
      f <- PhiDer(copula, w) - x[i]
      if (is.na(f)) f <- -Inf
      if (abs(f) < epsY || (hi - lo) < epsX) {
        status <- as.integer(it)
        break
      }
      if (f < 0) lo <- w else hi <- w
    }
    result[i] <- w
    conv[i] <- status
  }
  list(x, n, delta, theta, result, family, epsX, epsY, maxit, conv)
}
