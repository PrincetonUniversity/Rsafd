# Rsafd 1.1.0

## Packaging

* The repository is a standard R source package again: `R/*.R`, `man/*.Rd` and
  `data/*.rda` replace the installed-package snapshot (`R/Rsafd.rdb`, `help/`,
  `Meta/`, ...) shipped in 1.0.x. The source was rebuilt from the 2013 source in
  the git history plus every later change recovered from the 1.0.x binary, and
  verified to be identical to that binary: all functions, S4 generics, methods
  and classes, all 126 data sets, and all 111 help topics.
* Releases now use `X.Y.Z` version numbers and `vX.Y.Z` tags, and are published
  only after `R CMD check --as-cran` passes on macOS, Windows and Ubuntu with
  the current and previous R release. Each release provides the source
  package, Windows and macOS binaries, and a Jupyter image.
* The Jupyter image formerly built from the separate `rsafd-docker` repository
  is now built from this repository (`docker/jupyter/`) and published as
  `ghcr.io/princetonuniversity/rsafd`. The old
  `ghcr.io/princetonuniversity/rsafd-docker` image is no longer updated.
* Dependencies that the code never used (`glasso`, `qgraph`, `quadprog`,
  `scatterplot3d`, `splines`) are no longer imported. The package no longer
  needs compilation.

## Bug fixes

These functions failed with an error in 1.0.x. They now run, with their
computations unchanged:

* `rcopula()` for the BB1, BB2, BB3, BB6 and BB7 copulas, and the internal
  `Pearson.test()`, called C routines whose source was never part of the
  package. They now use R implementations with the same interface.
* `A()` for empirical copulas referred to the result of a removed C call.
* `bscall()` with a vector of volatilities, as in the book's example, failed on
  R >= 4.2 ("the condition has length > 1"). The zero-volatility test now uses
  the first element, which is what earlier versions of R did implicitly.
* `timeSeries()` with character positions called `.whichFormat()`, which the
  timeDate package has since renamed `whichFormat()`.
* `dmvnorm()` failed whenever `cov` was supplied.
* `Spearmans.rho()` for normal copulas used `$` on an S4 object.
* `merge.timeSeries()` with `units =` called a nonexistent replacement function.
* `twoDkreg(PLOT = TRUE)` called `rgl::persp3d`; it now plots with
  `plot3D::persp3D`, as `kdest()` already did.
* `fit.gpd()` with an unknown estimation method now reports that method instead
  of failing on an undefined variable.
* `as.ts()` on a `timeSeries` is registered as an S3 method.
* `gev.ml()`, `gpd.ml()`, `upper.shape.plot()`/`shape.plot()` and the copula
  `Kendalls.tau()`/`Spearmans.rho()` methods no longer create the temporary
  objects `tempX`, `tempN`, `tempCop`, `tempData` and `tempEstFun` in the
  user's workspace.
