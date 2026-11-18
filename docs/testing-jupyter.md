# Testing the Jupyter environment

Quick checks that a new Docker container or conda `Renv` environment works. Start Jupyter, open a notebook with the **R** kernel, and run each cell. Plots appear inline. The expected results come from the `1.2.0` image.

Every check should take a few seconds, except the first TensorFlow import (cell 8), which can take up to a minute.

## R kernel

**1. Package loads.** Messages about masking `qqnorm` and `CO2` are normal.

```r
library(Rsafd)
packageVersion("Rsafd")
```

Expect the release you installed, e.g. `'1.2.0'`.

**2. Data sets are available** (lazy loaded, no `data()` needed).

```r
length(DSPLRet)    # 15256
dim(TRGSP)         # 5000 6
head(Charlotte.ts) # daily MinT, MaxT, AvgT from 1961-01-01
```

**3. Option pricing matches the book.** This is the book's vectorized example over 5000 rows.

```r
BSTRG <- bscall(TAU = TRGSP[, 3] / 252, K = TRGSP[, 2], S = TRGSP[, 1],
                R = TRGSP[, 4], SIG = TRGSP[, 5])
range(TRGSP[, 6] - BSTRG)                                   # within about 1e-12 of 0
p <- bscall(TAU = 0.25, K = 100, S = 105, R = 0.03, SIG = 0.2)
p                                                           # 7.574
isig(TAU = 0.25, K = 100, S = 105, R = 0.03, CALL = p)      # about 0.2
```

**4. Heavy tails.** Fits both tails of the S&P 500 log-returns and draws two plots.

```r
fit <- fit.gpd(DSPLRet, tail = "two")
c(fit@upper.par.ests["xi"], fit@lower.par.ests["xi"])      # about 0.404 and 0.105
tailplot(fit)
```

**5. Copulas.** The sample and theoretical Kendall's tau should agree.

```r
set.seed(1)
u <- rcopula(gumbel.copula(2), 2000)
cor(u$x, u$y, method = "kendall")                           # about 0.50
Kendalls.tau(gumbel.copula(2))                              # 0.5

b <- rcopula(bb1.copula(0.5, 1.5), 2000)
c(cor(b$x, b$y, method = "kendall"), Kendalls.tau(bb1.copula(0.5, 1.5)))  # about 0.47 and 0.467
```

**6. Time series plot.**

```r
plot(Charlotte.ts)
```

**7. 3D plot** (uses plot3D, which the Docker image and `Renv` include).

```r
set.seed(2)
X <- rnorm(500); Y <- X + rnorm(500)
k <- kdest(X, Y)                                            # draws a 3D density surface
```

**8. Python from R** (reticulate and TensorFlow).

```r
library(reticulate)
tf <- import("tensorflow")
tf$`__version__`                                            # e.g. "2.21.0"
```

## Python kernel

Open a notebook with the **Python 3** kernel:

```python
import tensorflow as tf, keras, numpy, pandas, sklearn
print(tf.__version__, keras.__version__)
```

## Docker only

- **Container health:** `docker ps` should show the container as `healthy` about 2 minutes after it starts. The healthcheck tests TensorFlow, Rsafd and reticulate.
- **Installs persist:** in an R notebook, run `install.packages("xts", repos = "https://cloud.r-project.org")`. Then restart the container with the same `-v` folder and run `library(xts)`. It should load without reinstalling.

## If something fails

- **`library(Rsafd)` fails:** run `.libPaths()` and check that the environment's library is listed. For conda, make sure you ran `conda activate Renv` before `jupyter notebook`.
- **A 3D plot reports that plot3D is missing:** run `install.packages("plot3D")`. With a plain R install on macOS, also install [XQuartz](https://www.xquartz.org).
- **Dead kernel:** give Docker more memory (see the README's Docker section).
- **Anything else:** open an issue at https://github.com/PrincetonUniversity/Rsafd/issues and include the cell and its output.
