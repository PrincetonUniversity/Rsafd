# Rsafd

Rsafd is an R package that accompanies René Carmona's book *Statistical Analysis of Financial Data in R*. It gathers the datasets and reusable functions referenced throughout the text so readers can reproduce examples, explore additional analyses, and adapt the techniques to their own financial time-series problems.

## Package highlights

- **Curated financial datasets** spanning equity indexes, interest rates, energy markets, bond curves, and alternative assets (e.g., `DSP`, `SPFUT`, `FRWRD`, `us.bis.yield`, `bitcoin`).
- **Time-series utilities** for filtering, decomposition, forecasting, and volatility modeling (`kalman`, `sstl`, `pred.ar`, `garch`, `sim.garch`).
- **Extremes and heavy tails** with generalized extreme value and generalized Pareto distributions, L-moments, tail plots, and POT estimation (`gev.ml`, `fit.gpd`, `tailplot`).
- **Dependence modeling** through copulas, multivariate normals, kernel density/regression, and correlation diagnostics (`fit.copula`, `rcopula`, `dmvnorm`, `kreg`).
- **Option pricing and fixed-income helpers** including Black–Scholes, implied volatility, Nelson–Siegel forward curves, and coupon bond analytics (`bscall`, `isig`, `fns`, `bns`).

A complete index of exported symbols is available in `INDEX`, while `man/` hosts the Rd documentation files used by R's help system.

## Installation

```r
# Install from the latest GitHub release
if (.Platform$OS.type == "windows") {
  install.packages(
    "https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/Rsafd.zip",
    repos = NULL,
    type = "source"
  )
} else {
  install.packages(
    "https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/Rsafd.tar.gz",
    repos = NULL,
    type = "source"
  )
}

# Or install from a local checkout
install.packages("/path/to/Rsafd", repos = NULL, type = "source")
```

The package requires R 4.2.0 or newer. System requirements match the dependencies declared in `DESCRIPTION`:

- **Imports:** glasso, MASS, plot3D, qgraph, quadprog, quantreg, robustbase, scatterplot3d, splines, timeDate, tseries
- **Suggests:** lattice

## Getting started

```r
library(Rsafd)

# Load one of the bundled datasets
data("DSP")
plot(DSP, type = "l", main = "Daily S&P 500")

# Fit a GARCH(1,1) model to S&P 500 log returns
returns <- diff(log(DSP))[-1]
fit <- garch(returns, order = c(1, 1))
summary(fit)

# Estimate a generalized Pareto tail on extreme returns
library(quantreg)
exceedances <- returns[returns < quantile(returns, 0.95)]
fit_gpd <- fit.gpd(exceedances, threshold = quantile(exceedances, 0.9))
plot(fit_gpd)
```

See `help(package = "Rsafd")` for a searchable list of all datasets and helper functions, and consult the book for guided exercises that pair with each object.

## Repository layout

- `R/` contains the lazily-loaded compiled R code.
- `data/` stores serialized datasets in R's `.rda` format.
- `man/` holds Rd documentation consumed by `?foo` help pages.
- `html/` supplies the pre-built reference index for R CMD help.
- `DESCRIPTION`, `NAMESPACE`, `LICENSE`, and `INDEX` provide standard package metadata.

## Building from source

```sh
R CMD build Rsafd
R CMD check Rsafd_1.0.1.tar.gz
```

Running `R CMD check` is recommended before submitting changes or distributing a new tarball.
