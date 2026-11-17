# Rsafd

[![R-CMD-check](https://github.com/PrincetonUniversity/Rsafd/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/PrincetonUniversity/Rsafd/actions/workflows/R-CMD-check.yml)

Rsafd is an R package that accompanies René Carmona's book *Statistical Analysis of Financial Data in R*. It gathers the datasets and reusable functions referenced throughout the text so readers can reproduce examples, explore additional analyses, and adapt the techniques to their own financial time-series problems.

## Package highlights

- **Curated financial datasets** spanning equity indexes, interest rates, energy markets, bond curves, and alternative assets (e.g., `DSP`, `SPFUT`, `FRWRD`, `us.bis.yield`, `bitcoin`).
- **Time-series utilities** for filtering, decomposition, forecasting, and volatility modeling (`kalman`, `sstl`, `pred.ar`, `garch`, `sim.garch`).
- **Extremes and heavy tails** with generalized extreme value and generalized Pareto distributions, L-moments, tail plots, and POT estimation (`gev.ml`, `fit.gpd`, `tailplot`).
- **Dependence modeling** through copulas, multivariate normals, kernel density/regression, and correlation diagnostics (`fit.copula`, `rcopula`, `dmvnorm`, `kreg`).
- **Option pricing and fixed-income helpers** including Black–Scholes, implied volatility, Nelson–Siegel forward curves, and coupon bond analytics (`bscall`, `isig`, `fns`, `bns`).

Run `help(package = "Rsafd")` for the full index of functions and data sets.

## Installation

Rsafd requires R 4.2.0 or newer. Each release is tested on macOS, Windows and Ubuntu with the current and the previous R release.

```r
# Dependencies (from CRAN)
install.packages(c("MASS", "plot3D", "quantreg", "robustbase", "timeDate", "tseries"))

# Rsafd from the latest GitHub release
if (.Platform$OS.type == "windows") {
  # prebuilt Windows binary
  install.packages(
    "https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/Rsafd.zip",
    repos = NULL
  )
} else {
  # source package (no compiler needed)
  install.packages(
    "https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/Rsafd.tar.gz",
    repos = NULL, type = "source"
  )
}
```

A specific version can be installed from its release page, e.g.
`https://github.com/PrincetonUniversity/Rsafd/releases/download/v1.1.0/Rsafd_1.1.0.tar.gz`.
Each release also attaches a macOS binary (`Rsafd_X.Y.Z.tgz`) and `SHA256SUMS`.

## Getting started

```r
library(Rsafd)

# Bundled data sets are available directly (lazy loaded)
plot(DSP, type = "l", main = "Daily S&P 500")

# Generalized Pareto fit of both tails of the S&P 500 daily log-returns
fit <- fit.gpd(DSPLRet, tail = "two")
tailplot(fit)

# Black-Scholes price and implied volatility
price <- bscall(TAU = 0.25, K = 100, S = 105, R = 0.03, SIG = 0.2)
isig(TAU = 0.25, K = 100, S = 105, R = 0.03, CALL = price)

# Simulate from a Gumbel copula and estimate Kendall's tau
u <- rcopula(gumbel.copula(2), 1000)
cor(u$x, u$y, method = "kendall")
Kendalls.tau(gumbel.copula(2))
```

See the book for guided exercises that pair with each object.

## Jupyter environment (Docker)

A ready-to-use, multi-architecture (Intel and Apple Silicon macOS, Linux, Windows/WSL) Jupyter image with Python and R kernels and Rsafd pre-installed is published with every release:

* R, Rsafd, IRkernel, reticulate, keras, rgl, glmnet
* Python 3.12, TensorFlow, Keras, and a scientific stack (numpy, pandas, scipy, scikit-learn, matplotlib, seaborn), pinned by a hashed lockfile
* reticulate configured for the Python virtualenv; a healthcheck validates TensorFlow, Rsafd and the reticulate bridge

### Quick start

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
2. Run the container from a terminal:
   - macOS/Linux
     ```
     docker run -p 8888:8888 -e JUPYTER_LINK_ONLY=1 -v "$HOME":/workspace/notebooks ghcr.io/princetonuniversity/rsafd:latest
     ```
   - Windows (PowerShell)
     ```
     docker run -p 8888:8888 -e JUPYTER_LINK_ONLY=1 -v "${env:USERPROFILE}:/workspace/notebooks" ghcr.io/princetonuniversity/rsafd:latest
     ```
3. Open the printed URL in a browser.

Replace `latest` with a version (e.g. `1.1.0`) to pin the environment to a specific Rsafd release.

> The image used to be published as `ghcr.io/princetonuniversity/rsafd-docker`. That image is no longer updated; use `ghcr.io/princetonuniversity/rsafd` instead.

| Capability | Details |
|------------|---------|
| UI modes | `-e JUPYTER_UI=notebook` (default) or `-e JUPYTER_UI=lab` |
| Token control | Random by default; set a reusable token with `-e JUPYTER_TOKEN=yourtoken` or disable it with `-e JUPYTER_DISABLE_TOKEN=1` |
| Healthcheck | `/usr/local/bin/healthcheck`, wired into Docker `HEALTHCHECK` |
| Persistent user installs | Packages installed from a notebook persist in the mounted volume (see below) |

### Persistent package installs

When a host directory is mounted at `/workspace/notebooks`, packages installed interactively are kept in a directory specific to the image's commit:

```
/workspace/notebooks/.rsafd-docker-<image-sha>/
    py/   # Python user packages (PIP_TARGET)
    R/    # R user library (R_LIBS_USER)
```

```python
!pip install statsmodels
```

```r
install.packages("xts", repos = "https://cloud.r-project.org")
```

Remove `.rsafd-docker-*/` on the host to reset.

### Troubleshooting

* **Dead kernel:** large data sets may exceed the container's memory. Increase it in [Docker Desktop's settings](https://docs.docker.com/desktop/settings-and-maintenance/settings/) (macOS) or the [WSL configuration](https://learn.microsoft.com/en-us/windows/wsl/wsl-config) (Windows).
* **`install.packages()` hangs:** pass `repos = "https://cloud.r-project.org"`, since Jupyter cannot answer the interactive CRAN-mirror prompt.

## Development

Everything runs in pinned containers (R 4.6.2 with a dated CRAN snapshot); only Docker is needed on the host.

```sh
docker-compose build check
docker-compose run --rm check     # R CMD check --as-cran; fails on any WARNING or ERROR
docker-compose run --rm test      # testthat suite
docker-compose run --rm lint      # lintr
docker-compose up jupyter         # Jupyter image built from this checkout, on 127.0.0.1:8888
```

With R installed locally, the usual `R CMD build .` and `R CMD check --as-cran Rsafd_*.tar.gz` also work.

Repository layout:

- `R/`, `man/`, `data/`, `tests/`: the package (standard R source layout).
- `Dockerfile`, `docker-compose.yml`, `scripts/container-check.sh`: pinned check/test/lint container.
- `docker/jupyter/`: the Jupyter image (`requirements.in` and its hashed `requirements.lock`; regenerate with `scripts/lock-python.sh`).
- `.github/workflows/`: CI (`R-CMD-check.yml`, `jupyter-image.yml`) and releases (`release.yml`).

The code in `R/Financial.R`, `R/Multivariate.R`, `R/Smoothing.R`, `R/timeSeries.R` and `R/Univariate.R` is the author's and keeps his style; `lintr` is applied to new code only.

### Continuous integration

Every push and pull request runs `R CMD check --as-cran` (warnings fail the build) on macOS, Windows and Ubuntu, each with the current (`release`) and previous (`oldrel-1`) R version, plus lint and the containerized check and tests. Changes to the Jupyter image also build it for amd64 and arm64.

### Releases

Versions follow `X.Y.Z` and are tagged `vX.Y.Z`. To release:

1. Add a `# Rsafd X.Y.Z` section to `NEWS.md` and merge to `main`.
2. Run `scripts/release.sh X.Y.Z --dry-run`, then `scripts/release.sh X.Y.Z`. This sets `Version:` in `DESCRIPTION`, commits, tags `vX.Y.Z` and pushes.
3. The `release` workflow verifies that the tag, `DESCRIPTION` and `NEWS.md` agree, reruns the full check matrix, and only then publishes the GitHub release (source package, Windows and macOS binaries, checksums) and the multi-arch Jupyter image (`X.Y.Z`, `X.Y`, `latest`).
