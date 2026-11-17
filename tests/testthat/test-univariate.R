test_that("GEV d/p/q functions are mutually consistent", {
  for (xi in c(-0.2, 0, 0.3)) {
    p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
    expect_equal(pgev(qgev(p, m = 1, lambda = 2, xi = xi), m = 1, lambda = 2, xi = xi), p,
                 tolerance = 1e-8, info = paste("xi =", xi))
    q <- qgev(c(0.001, 0.999), m = 1, lambda = 2, xi = xi)
    mass <- integrate(dgev, q[1], q[2], m = 1, lambda = 2, xi = xi)$value
    expect_equal(mass, 0.998, tolerance = 1e-4, info = paste("xi =", xi))
  }
  # xi = 0 is the Gumbel distribution
  expect_equal(pgev(0), exp(-1))
})

test_that("rgev draws follow pgev", {
  set.seed(10)
  x <- rgev(5000, m = 0, lambda = 1, xi = 0.2)
  expect_gt(stats::ks.test(x, function(q) pgev(q, m = 0, lambda = 1, xi = 0.2))$p.value, 1e-3)
})

test_that("Pareto d/p/q functions are mutually consistent", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(ppareto(qpareto(p, m = 0, lambda = 1, xi = 0.25), m = 0, lambda = 1, xi = 0.25),
               p, tolerance = 1e-8)
  set.seed(11)
  x <- rpareto(3000, m = 0, lambda = 1, xi = 0.25)
  expect_gt(stats::ks.test(x, function(q) ppareto(q, m = 0, lambda = 1, xi = 0.25))$p.value, 1e-3)
})

test_that("sample.LMOM matches the L-moments of the uniform distribution", {
  set.seed(12)
  l <- sample.LMOM(runif(20000))
  # lambda1 = 1/2, lambda2 = 1/6, L-skewness = 0
  expect_equal(unname(l[1]), 0.5, tolerance = 0.01)
  expect_equal(unname(l[2]), 1 / 6, tolerance = 0.01)
  expect_equal(unname(l[3]), 0, tolerance = 0.02)
})

test_that("gev.ml recovers parameters of a simulated sample", {
  set.seed(13)
  x <- rgev(3000, m = 2, lambda = 1.5, xi = 0.2)
  est <- suppressWarnings(gev.ml(x))$param.est
  expect_equal(unname(est[1:2]), c(2, 1.5), tolerance = 0.1)
  expect_equal(unname(est[3]), 0.2, tolerance = 0.25)
})

test_that("fit.gpd builds a two-tailed model whose quantiles invert its cdf", {
  set.seed(14)
  x <- rt(3000, df = 4)
  fit <- fit.gpd(x, plot = FALSE)
  expect_s4_class(fit, "gpd.two.tails")
  p <- c(0.01, 0.2, 0.5, 0.8, 0.99)
  expect_equal(pgpd(fit, qgpd(fit, p)), p, tolerance = 1e-3)
  # heavy tails of the Student t are detected
  expect_gt(unname(fit@upper.par.ests["xi"]), 0)
  expect_gt(unname(fit@lower.par.ests["xi"]), 0)
  r <- rgpd(fit, 100)
  expect_length(r, 100)
})

test_that("block.max splits data into blocks and fits a GEV to the maxima", {
  set.seed(16)
  x <- rnorm(250 * 400)
  bm <- suppressWarnings(block.max(x, block.size = 250))
  expect_s3_class(bm, "block.max")
  expect_equal(bm$nb.blocks, 400)
  expect_true(all(bm$interval.b - bm$interval.a + 1 == 250))
  expect_named(bm$param.est, c("m", "lambda", "xi"))
  # the GEV location for maxima of 250 standard normals is close to 2.65
  expect_equal(unname(bm$param.est["m"]), 2.65, tolerance = 0.05)
})

test_that("plotting functions run on a headless device", {
  set.seed(15)
  x <- rt(1000, df = 4)
  expect_no_error(with_null_device({
    qqexp(abs(x))
    eda.shape(x)
    shape.plot(x, tail = "upper")
    tailplot(fit.gpd(x, plot = FALSE))
  }))
})
