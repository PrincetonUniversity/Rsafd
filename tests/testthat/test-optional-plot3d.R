# plot3D is optional: it loads tcltk (via misc3d), which needs XQuartz on macOS.

test_that("loading Rsafd does not load plot3D or tcltk", {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressMessages(library(Rsafd))",
    "cat(c('plot3D', 'tcltk') %in% loadedNamespaces())"
  ), script)
  out <- system2(file.path(R.home("bin"), "Rscript"), script, stdout = TRUE, stderr = FALSE,
                 env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)))
  expect_identical(utils::tail(out, 1), "FALSE FALSE")
})

test_that("3D plots explain how to install plot3D when it is missing", {
  local_mocked_bindings(.plot3d_available = function() FALSE)
  set.seed(7)
  X <- cbind(runif(50), runif(50))
  Y <- X[, 1] + rnorm(50, sd = 0.1)
  expect_error(
    with_null_device(twoDkreg(X, Y, B = c(0.2, 0.2), N = 10, PLOT = TRUE)),
    "install.packages(\"plot3D\")", fixed = TRUE
  )
  # computations that do not plot still work
  expect_no_error(twoDkreg(X, Y, B = c(0.2, 0.2), N = 10, PLOT = FALSE))
})
