# syntax=docker/dockerfile:1
# Pinned R toolchain for building, checking, testing and linting Rsafd.
#   docker-compose run --rm check   # R CMD check --as-cran (fails on WARNING)
#   docker-compose run --rm test    # testthat suite
#   docker-compose run --rm lint    # lintr
FROM rocker/r-ver:4.6.2 AS deps

# Freeze CRAN at a dated Posit Package Manager snapshot so builds are reproducible.
ARG CRAN_SNAPSHOT=2026-11-16
RUN echo "options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/resolute/${CRAN_SNAPSHOT}'))" \
      >> /usr/local/lib/R/etc/Rprofile.site

RUN apt-get update \
 && apt-get install -y --no-install-recommends qpdf pandoc xvfb xauth libtcl8.6 libtk8.6 \
 && rm -rf /var/lib/apt/lists/*

# Install everything DESCRIPTION asks for (Imports + Suggests) plus the check tooling.
COPY DESCRIPTION /tmp/Rsafd/DESCRIPTION
COPY scripts/install-deps.R /tmp/install-deps.R
RUN Rscript /tmp/install-deps.R /tmp/Rsafd/DESCRIPTION rcmdcheck lintr pkgload

FROM deps AS check
WORKDIR /work
COPY . /work/Rsafd
COPY scripts/container-check.sh /usr/local/bin/rsafd-check
RUN chmod +x /usr/local/bin/rsafd-check
ENTRYPOINT ["/usr/local/bin/rsafd-check"]
CMD ["check"]
