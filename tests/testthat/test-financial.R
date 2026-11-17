test_that("bscall matches the Black-Scholes closed form", {
  TAU <- 0.5
  K <- 100
  S <- 105
  R <- 0.03
  SIG <- 0.2
  d1 <- (log(S / K) + TAU * (R + SIG^2 / 2)) / (SIG * sqrt(TAU))
  d2 <- d1 - SIG * sqrt(TAU)
  expected <- S * pnorm(d1) - K * exp(-R * TAU) * pnorm(d2)
  expect_equal(bscall(TAU = TAU, K = K, S = S, R = R, SIG = SIG), expected)
  # zero volatility gives the discounted intrinsic value
  expect_equal(bscall(TAU = TAU, K = K, S = S, R = R, SIG = 0), S - exp(-R * TAU) * K)
})

test_that("isig inverts bscall", {
  price <- bscall(TAU = 0.25, K = 1, S = 1.02, R = 0.05, SIG = 0.3)
  expect_equal(isig(TAU = 0.25, K = 1, S = 1.02, R = 0.05, CALL = price), 0.3,
               tolerance = 1e-4)
})

test_that("fns is the Nelson-Siegel forward curve", {
  THETA <- c(0.06, -0.02, 0.01, 2)
  x <- c(0.5, 1, 5, 10)
  expect_equal(fns(x, THETA),
               THETA[1] + (THETA[2] + THETA[3] * x) * exp(-x / THETA[4]))
})

test_that("yns approaches the long rate and is consistent with fns", {
  THETA <- c(0.06, -0.02, 0.01, 2)
  expect_equal(yns(200, THETA), THETA[1], tolerance = 1e-3)
  # the yield is the average forward rate: x * y(x) = integral of f
  x <- 3
  expect_equal(x * yns(x, THETA), integrate(fns, 0, x, THETA = THETA)$value,
               tolerance = 1e-6)
})
