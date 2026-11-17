# syntax=docker/dockerfile:1
## Multi-arch Jupyter (Python + R) environment with Rsafd installed from this repository.
## Targets linux/amd64 and linux/arm64 (Intel and Apple Silicon macOS, Linux, Windows WSL).
## Build from the repository root (the build context must contain the package source):
##   docker build -f docker/jupyter/Dockerfile -t rsafd-jupyter .
## or: docker-compose build jupyter
## Published by .github/workflows/release.yml as ghcr.io/princetonuniversity/rsafd.

FROM ghcr.io/astral-sh/uv:0.9.5 AS uv

FROM rocker/r-ver:4.6.2

ARG GIT_SHA=unknown
ARG CRAN_SNAPSHOT=2026-11-16
ARG PYTHON_VERSION=3.12
ARG USERNAME=developer
ARG USERID=1000
ARG GROUPID=1000
ARG GH_OWNER=PrincetonUniversity

LABEL org.opencontainers.image.source="https://github.com/${GH_OWNER}/Rsafd" \
      org.opencontainers.image.description="Multi-arch Jupyter (Python + R) with Rsafd, tensorflow, keras, reticulate" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.revision="${GIT_SHA}"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

## Freeze CRAN at a dated snapshot so every build resolves the same R packages.
RUN echo "options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/resolute/${CRAN_SNAPSHOT}'))" \
      >> /usr/local/lib/R/etc/Rprofile.site

## System packages: container plumbing plus runtime libraries for the R packages below.
RUN apt-get update -y \
 && apt-get install -y --no-install-recommends \
      ca-certificates curl git openssl sudo tini nodejs \
      libzmq5 libtcl8.6 libtk8.6 libglpk40 libgmp10 libxml2-16 libpng16-16t64 libjpeg-turbo8 \
      libfreetype6 libfontconfig1 libgl1 libglu1-mesa libx11-6 libxext6 libxrender1 libxt6 \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

## Create non-root user (matching host UID/GID optionally for mounted home)
RUN set -euo pipefail; \
	if getent group ${GROUPID} >/dev/null; then \
		EXISTING_GROUP_NAME=$(getent group ${GROUPID} | cut -d: -f1); \
		echo "GID ${GROUPID} already exists as ${EXISTING_GROUP_NAME}"; \
	else \
		groupadd -g ${GROUPID} ${USERNAME}; \
	fi; \
	TARGET_GROUP_NAME=$(getent group ${GROUPID} | cut -d: -f1); \
	if id -u ${USERID} >/dev/null 2>&1; then \
		echo "UID ${USERID} already exists, will re-use and ensure username ${USERNAME}"; \
		if ! id -u ${USERNAME} >/dev/null 2>&1; then \
			usermod -l ${USERNAME} $(getent passwd ${USERID} | cut -d: -f1) || true; \
			usermod -d /home/${USERNAME} -m ${USERNAME} || true; \
		fi; \
	elif id -u ${USERNAME} >/dev/null 2>&1; then \
		echo "User ${USERNAME} exists with different UID; not modifying UID"; \
	else \
		useradd -m -s /bin/bash -u ${USERID} -g ${TARGET_GROUP_NAME} ${USERNAME}; \
	fi; \
	usermod -aG sudo ${USERNAME} && echo "%sudo ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/99-nopasswd

WORKDIR /workspace

## Python virtual environment from the hash-pinned lockfile (uv-managed CPython, PEP 668 friendly)
COPY --from=uv /uv /usr/local/bin/uv
COPY docker/jupyter/requirements.lock /tmp/requirements.lock
ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python UV_NO_CACHE=1
RUN uv python install ${PYTHON_VERSION} \
 && uv venv /opt/venv --python ${PYTHON_VERSION} \
 && uv pip install --python /opt/venv/bin/python --require-hashes -r /tmp/requirements.lock \
 && /opt/venv/bin/python -m ipykernel install --name python3 --display-name "Python 3 (venv)" --sys-prefix

## Environment variables for reticulate to find venv python
ENV PATH="/opt/venv/bin:$PATH" \
    RETICULATE_PYTHON=/opt/venv/bin/python

## R: Jupyter kernel, reticulate bridge, and the packages used alongside Rsafd in the course
COPY scripts/install-deps.R /tmp/install-deps.R
COPY DESCRIPTION /tmp/Rsafd-desc/DESCRIPTION
RUN Rscript /tmp/install-deps.R /tmp/Rsafd-desc/DESCRIPTION \
      IRkernel remotes reticulate \
      quadprog scatterplot3d glasso qgraph keras rgl glmnet \
 && R -q -e "IRkernel::installspec(user = FALSE, name = 'ir', displayname = 'R')"

## Rsafd itself, built and installed from this checkout's source
COPY . /tmp/Rsafd
RUN cd /tmp && R CMD build --no-manual Rsafd \
 && R CMD INSTALL Rsafd_*.tar.gz \
 && rm -rf /tmp/Rsafd /tmp/Rsafd_*.tar.gz \
 && R -q -e "library(Rsafd); cat('Rsafd', format(packageVersion('Rsafd')), 'installed\n')"

## Basic import tests to fail early if something is broken
RUN /opt/venv/bin/python - <<'PYTEST' && \
	R -q -e 'library(Rsafd); cat("Rsafd loaded successfully (smoke test)\n")' && \
	R -q -e 'library(reticulate); tf <- import("tensorflow"); cat(tf[["__version__"]], "\n")'
import tensorflow as tf
import keras
print('TensorFlow:', tf.__version__)
print('Keras:', keras.__version__)
PYTEST

## Expose Jupyter default port
EXPOSE 8888

## Make sure notebooks & data written with correct ownership when running as non-root
RUN mkdir -p /workspace/notebooks && chown -R ${USERNAME}:${GROUPID} /workspace

COPY docker/jupyter/start-jupyter.sh /usr/local/bin/start-jupyter
COPY docker/jupyter/healthcheck.sh /usr/local/bin/healthcheck
RUN chmod +x /usr/local/bin/start-jupyter /usr/local/bin/healthcheck

USER ${USERNAME}
ENV HOME=/home/${USERNAME}

ENV RSAFD_IMAGE_SHA=${GIT_SHA} \
	JUPYTER_PORT=8888 \
	JUPYTER_DISABLE_TOKEN=0 \
	JUPYTER_TOKEN="" \
	JUPYTER_LAB_ARGS="" \
	JUPYTER_PLAIN_URL=0 \
	JUPYTER_DISABLE_LSP=1 \
	JUPYTER_LINK_ONLY=0 \
	JUPYTER_UI=notebook

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["/usr/local/bin/start-jupyter"]

HEALTHCHECK --interval=2m --timeout=20s --start-period=30s --retries=3 CMD /usr/local/bin/healthcheck || exit 1
