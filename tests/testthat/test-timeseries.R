make_ts <- function(values, start = "2021-03-01", units = "v") {
  pos <- timeDate::timeSequence(from = start, length.out = length(values), by = "day")
  timeSeries(positions = pos, data = values, units = units)
}

test_that("timeSeries objects expose their components", {
  x <- make_ts(c(10, 11, 13, 12))
  expect_true(is.timeSeries(x))
  expect_s4_class(x, "timeSeries")
  expect_equal(as.numeric(seriesData(x)), c(10, 11, 13, 12))
  expect_identical(seriesUnits(x), "v")
  expect_length(seriesPositions(x), 4)
})

test_that("arithmetic and Math group generics act on the data", {
  x <- make_ts(c(1, 2, 4))
  expect_equal(as.numeric(seriesData(x * 2)), c(2, 4, 8))
  expect_equal(as.numeric(seriesData(log(x))), log(c(1, 2, 4)))
})

test_that("diff of a timeSeries differences the data", {
  x <- make_ts(c(1, 4, 9, 16))
  expect_equal(as.numeric(seriesData(diff(x))), c(3, 5, 7))
})

test_that("bundled time series datasets load as timeSeries objects", {
  for (nm in c("DSP.ts", "co2.ts", "Charlotte.ts")) {
    e <- new.env()
    utils::data(list = nm, package = "Rsafd", envir = e)
    expect_s4_class(get(nm, e), "timeSeries")
  }
})

test_that("every dataset in data/ loads and is non-empty", {
  ds <- utils::data(package = "Rsafd")$results[, "Item"]
  expect_gt(length(ds), 100)
  for (nm in ds) {
    e <- new.env()
    utils::data(list = nm, package = "Rsafd", envir = e)
    obj <- get(nm, e)
    expect_gt(length(obj), 0, label = nm)
  }
})

test_that("sstl decomposes a seasonal series", {
  set.seed(30)
  n <- 3 * 365
  x <- make_ts(sin(2 * pi * (1:n) / 365) + (1:n) / n + rnorm(n, sd = 0.1))
  res <- sstl(x)
  expect_true(all(c("trend", "sea", "rem") %in% names(res)))
})
