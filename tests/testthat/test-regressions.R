# Regression tests for defects present in the 1.0.x binary release.

test_that("estimation routines do not write into the global environment", {
  set.seed(1)
  x <- rgev(200, m = 0, lambda = 1, xi = 0.2)
  before <- global_names()
  suppressWarnings(gev.ml(x))
  gev.lmom(sample.LMOM(x))
  gpd.ml(abs(x))
  Spearmans.rho(frank.copula(3))
  # gev.lmom still assigns `t0` globally; changing that is tied to the pending
  # decision tested below, so it is excluded here.
  expect_identical(setdiff(global_names(), "t0"), setdiff(before, "t0"))
})

test_that("gev.lmom solves Hosking's equation outside the approximation range", {
  skip("Pending maintainer decision: gev.lmom solves with t0 = 1 instead of (t3 + 3)/2")
  # For L-skewness outside [-0.1, 0.5] the shape k must satisfy
  # the ratio (1 - 3^-k) / (1 - 2^-k) equal to (3 + t3) / 2.
  for (t3 in c(0.7, -0.3)) {
    k <- suppressWarnings(gev.lmom(c(0, 1, t3)))$param.est[3]
    if (SHAPE.XI) k <- -k
    lhs <- (1 - 3^(-k)) / (1 - 2^(-k))
    expect_equal(unname(lhs), (3 + t3) / 2, tolerance = 1e-4)
  }
})

test_that("rcopula works for the BB Archimedean families (no compiled code needed)", {
  fams <- list(bb1 = bb1.copula(0.5, 1.5), bb2 = bb2.copula(0.3, 0.5),
               bb3 = bb3.copula(1.5, 0.5), bb6 = bb6.copula(1.5, 1.5),
               bb7 = bb7.copula(1.5, 0.5))
  for (nm in names(fams)) {
    set.seed(42)
    s <- suppressWarnings(rcopula(fams[[nm]], 2000))
    expect_length(s$x, 2000)
    expect_true(all(s$x > 0 & s$x < 1 & s$y > 0 & s$y < 1), info = nm)
    expect_gt(stats::ks.test(s$y, "punif")$p.value, 1e-3)
    expect_equal(stats::cor(s$x, s$y, method = "kendall"),
                 Kendalls.tau(fams[[nm]]), tolerance = 0.05, info = nm)
  }
})

test_that("Pickands A() of an empirical copula is computed in pure R", {
  set.seed(3)
  u <- rcopula(gumbel.copula(2), 500)
  ec <- empirical.copula(u$x, u$y)
  a <- suppressWarnings(A(ec, seq(0.05, 0.95, by = 0.05)))
  expect_length(a, 19)
  expect_true(all(a >= 0.5 & a <= 1))
})

test_that("Pickands A() of an empirical copula tracks the source copula", {
  skip(paste("Pending maintainer decision: bin index is clamped at N2 - 1",
             "(evaluation points) rather than n - 1 (sample size)"))
  set.seed(3)
  u <- rcopula(gumbel.copula(2), 500)
  t <- seq(0.05, 0.95, by = 0.05)
  a <- suppressWarnings(A(empirical.copula(u$x, u$y), t))
  expect_lt(max(abs(a - A(gumbel.copula(2), t))), 0.02)
})

test_that("Pearson.test runs without compiled code", {
  set.seed(4)
  cop <- normal.copula(0.5)
  u <- rcopula(cop, 1000)
  res <- Rsafd:::Pearson.test(empirical.copula(u$x, u$y), cop, Ncells = 25)
  expect_s3_class(res, "data.frame")
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_gt(res$p.value, 1e-4) # data drawn from the tested copula
})

test_that("fit.gpd reports an unknown estimation method clearly", {
  set.seed(5)
  expect_error(fit.gpd(rnorm(500), tail = "upper", upper.method = "nope", plot = FALSE),
               "Unknown method")
})

test_that("twoDkreg can plot without the rgl package", {
  set.seed(6)
  X <- cbind(runif(100), runif(100))
  Y <- X[, 1] + X[, 2] + rnorm(100, sd = 0.1)
  expect_no_error(with_null_device(twoDkreg(X, Y, B = c(0.2, 0.2), N = 20, PLOT = TRUE)))
})

test_that("merge.timeSeries applies replacement units", {
  a <- timeSeries(positions = timeDate::timeDate(c("2020-01-01", "2020-01-02")),
                  data = matrix(1:2), units = "a")
  b <- timeSeries(positions = timeDate::timeDate(c("2020-01-02", "2020-01-03")),
                  data = matrix(3:4), units = "b")
  m <- merge.timeSeries(a, b, units = c("x", "y"))
  expect_identical(seriesUnits(m), c("x", "y"))
  expect_identical(colnames(seriesData(m)), c("x", "y"))
})

test_that("as.ts dispatches to the timeSeries method", {
  x <- timeSeries(positions = timeDate::timeDate(c("2020-01-01", "2020-01-02", "2020-01-03")),
                  data = c(1, 2, 3), units = "v")
  expect_s3_class(as.ts(x), "ts")
})

test_that("bscall runs the book's vectorised example on R >= 4.2", {
  # Book example (bscall.Rd): one call over every row of TRGSP.
  bs <- bscall(TAU = TRGSP[, 3] / 252, K = TRGSP[, 2], S = TRGSP[, 1],
               R = TRGSP[, 4], SIG = TRGSP[, 5])
  expect_length(bs, nrow(TRGSP))
  # same values as pricing each row on its own
  one <- function(i) {
    bscall(TAU = TRGSP[i, 3] / 252, K = TRGSP[i, 2], S = TRGSP[i, 1],
           R = TRGSP[i, 4], SIG = TRGSP[i, 5])
  }
  expect_equal(bs[c(1, 50, nrow(TRGSP))], vapply(c(1, 50, nrow(TRGSP)), one, numeric(1)))
})

test_that("timeSeries accepts character positions", {
  x <- timeSeries(positions = c("2020-01-01", "2020-01-02", "2020-01-03"), data = 1:3)
  expect_s4_class(x, "timeSeries")
  expect_equal(as.character(seriesPositions(x)), c("2020-01-01", "2020-01-02", "2020-01-03"))
})
