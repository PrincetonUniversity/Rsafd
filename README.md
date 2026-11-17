# Rsafd

[![R-CMD-check](https://github.com/PrincetonUniversity/Rsafd/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/PrincetonUniversity/Rsafd/actions/workflows/R-CMD-check.yml)

Functions and data sets for René Carmona's *Statistical Analysis of Financial Data in R*: heavy tails and extreme values (`fit.gpd`, `gev.ml`, `tailplot`), copulas (`fit.copula`, `rcopula`), smoothing (`kreg`, `kdest`), time series (`garch`, `kalman`, `sstl`), option and bond pricing (`bscall`, `isig`, `bns`), and the book's data (`DSP`, `SPFUT`, `FRWRD`, ...). See `help(package = "Rsafd")`.

## Install

Choose one. Each is tested in CI on macOS (Apple Silicon), Windows and Linux.

### R (4.2 or newer)

```r
install.packages(c("MASS", "plot3D", "quantreg", "robustbase", "timeDate", "tseries"))
url <- "https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/"
if (.Platform$OS.type == "windows") {
  install.packages(paste0(url, "Rsafd.zip"), repos = NULL)
} else {
  install.packages(paste0(url, "Rsafd.tar.gz"), repos = NULL, type = "source")
}
```

On macOS, install [XQuartz](https://www.xquartz.org) first: plot3D loads R's tcltk, which needs it. The Miniconda route below does not.

### Miniconda (Jupyter with R and Python)

Creates an environment named `Renv` with R 4.5, Python 3.12, Jupyter, TensorFlow and Keras. Requires [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install). On Windows, run these commands in [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Ubuntu).

```sh
# macOS (Apple Silicon)
curl -fsSL https://raw.githubusercontent.com/PrincetonUniversity/Rsafd/main/conda/renv-osx-arm64.txt -o /tmp/renv.txt
# Linux or WSL
curl -fsSL https://raw.githubusercontent.com/PrincetonUniversity/Rsafd/main/conda/renv-linux-$(uname -m).txt -o /tmp/renv.txt

conda create -y -n Renv --override-channels -c conda-forge --file /tmp/renv.txt
conda activate Renv
Rscript -e 'install.packages("https://github.com/PrincetonUniversity/Rsafd/releases/latest/download/Rsafd.tar.gz", repos = NULL, type = "source")'
jupyter notebook
```

Later sessions: `conda activate Renv`, then `jupyter notebook`. Remove with `conda env remove -n Renv`.

### Docker (Jupyter with R and Python)

Requires [Docker Desktop](https://www.docker.com/products/docker-desktop/). Open the printed URL in a browser.

```sh
# macOS / Linux
docker run -p 8888:8888 -e JUPYTER_LINK_ONLY=1 -v "$HOME":/workspace/notebooks ghcr.io/princetonuniversity/rsafd:latest
# Windows (PowerShell)
docker run -p 8888:8888 -e JUPYTER_LINK_ONLY=1 -v "${env:USERPROFILE}:/workspace/notebooks" ghcr.io/princetonuniversity/rsafd:latest
```

- **Pin a release:** use a version tag such as `:1.1.0` instead of `:latest`.
- **Options:** `-e JUPYTER_UI=lab` for JupyterLab; `-e JUPYTER_TOKEN=...` sets the token, `-e JUPYTER_DISABLE_TOKEN=1` removes it.
- **Your installs persist:** packages added from a notebook (`!pip install`, `install.packages()`) are kept in `.rsafd-docker-<sha>/` in the mounted folder. Delete it to reset.
- **Dead kernel:** give Docker more memory ([macOS](https://docs.docker.com/desktop/settings-and-maintenance/settings/), [Windows/WSL](https://learn.microsoft.com/en-us/windows/wsl/wsl-config)).
- **`install.packages()` hangs:** add `repos = "https://cloud.r-project.org"`.

The image was formerly `ghcr.io/princetonuniversity/rsafd-docker`, which is no longer updated.

## Example

```r
library(Rsafd)
fit <- fit.gpd(DSPLRet, tail = "two")   # GPD fit to both tails of S&P 500 log-returns
tailplot(fit)
bscall(TAU = 0.25, K = 100, S = 105, R = 0.03, SIG = 0.2)
Kendalls.tau(gumbel.copula(2))
```

## Development

Only Docker is needed. R 4.6.2 and a dated CRAN snapshot are pinned.

```sh
docker-compose run --rm check   # R CMD check --as-cran; warnings fail
docker-compose run --rm test    # testthat
docker-compose run --rm lint    # lintr on new code; the author's R/ files keep his style
docker-compose up jupyter       # Jupyter image from this checkout
```

| Path | Contents |
|------|----------|
| `R/`, `man/`, `data/`, `tests/` | the package |
| `conda/` | `environment.yml` and per-platform lockfiles (`scripts/lock-conda.sh`) |
| `docker/jupyter/` | Jupyter image; Python lockfile via `scripts/lock-python.sh` |
| `Dockerfile`, `docker-compose.yml` | check, test and lint container |
| `.github/workflows/` | CI and releases |

**CI** runs `R CMD check --as-cran` on macOS, Windows and Ubuntu with R release and oldrel-1, plus lint and the container check. Changes to `conda/` or `docker/jupyter/` also build those environments.

**Releases** use `X.Y.Z`. Add a `# Rsafd X.Y.Z` section to `NEWS.md` and merge to `main`. Then run `scripts/release.sh X.Y.Z` (try `--dry-run` first). The tag triggers a full re-check, then publishes the GitHub release (source, Windows and macOS binaries, checksums) and the Jupyter image.
