test_that("dmvnorm matches the bivariate normal density", {
  sig <- matrix(c(1, 0.5, 0.5, 2), 2)
  x <- matrix(c(0.3, -0.4, 1, 1), ncol = 2, byrow = TRUE)
  inv <- solve(sig)
  expected <- apply(x, 1, function(z) exp(-0.5 * drop(z %*% inv %*% z)) / (2 * pi * sqrt(det(sig))))
  expect_equal(as.numeric(dmvnorm(x, mean = c(0, 0), cov = sig)), expected)
})

test_that("rmvnorm reproduces the requested covariance", {
  set.seed(20)
  sig <- matrix(c(1, 0.5, 0.5, 2), 2)
  z <- rmvnorm(20000, mean = c(1, -1), cov = sig)
  expect_equal(colMeans(z), c(1, -1), tolerance = 0.05)
  expect_equal(stats::cov(z), sig, tolerance = 0.05)
})

test_that("copula constructors validate parameters", {
  expect_error(bb1.copula(0.5, 0.5), "delta >= 1")
})

test_that("copula cdfs satisfy the boundary conditions", {
  u <- c(0.1, 0.4, 0.8)
  for (cop in list(normal.copula(0.5), frank.copula(3), gumbel.copula(2),
                   kimeldorf.sampson.copula(1.5), bb1.copula(0.5, 1.5))) {
    expect_equal(pcopula(cop, u, rep(1, 3)), u, tolerance = 1e-6, info = class(cop))
    expect_equal(pcopula(cop, rep(1, 3), u), u, tolerance = 1e-6, info = class(cop))
  }
})

test_that("Kendall's tau matches known closed forms", {
  expect_equal(Kendalls.tau(gumbel.copula(2)), 1 - 1 / 2, tolerance = 1e-4)
  expect_equal(Kendalls.tau(normal.copula(0.5)), 2 / pi * asin(0.5), tolerance = 1e-4)
  expect_equal(Kendalls.tau(kimeldorf.sampson.copula(2)), 2 / (2 + 2), tolerance = 1e-4)
})

test_that("Spearman's rho of the normal copula matches its closed form", {
  expect_equal(Spearmans.rho(normal.copula(0.5)), 6 / pi * asin(0.25), tolerance = 1e-6)
})

test_that("rcopula samples reproduce the copula's Kendall's tau", {
  for (cop in list(normal.copula(0.6), frank.copula(4), gumbel.copula(1.8))) {
    set.seed(21)
    s <- rcopula(cop, 3000)
    expect_equal(stats::cor(s$x, s$y, method = "kendall"), Kendalls.tau(cop),
                 tolerance = 0.04, info = class(cop))
  }
})

test_that("fit.copula recovers the parameter of a normal copula", {
  set.seed(22)
  s <- rcopula(normal.copula(0.6), 1000)
  fit <- with_null_device(fit.copula(empirical.copula(s$x, s$y), family = "normal",
                                     plotPicture = FALSE))
  expect_equal(as.numeric(fit@parameters[1]), 0.6, tolerance = 0.15)
})
