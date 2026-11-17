#!/usr/bin/env bash
# Regenerate docker/jupyter/requirements.lock from requirements.in for
# every platform (universal resolution, Python 3.12) using a pinned uv image.
# Safe to re-run; only rewrites the lockfile.
set -euo pipefail

UV_IMAGE="ghcr.io/astral-sh/uv:0.9.5-python3.12-bookworm-slim"
here="$(cd "$(dirname "$0")/.." && pwd)"

docker run --rm --entrypoint uv -v "$here/docker/jupyter":/w -w /w "$UV_IMAGE" \
  pip compile requirements.in \
    --universal --python-version 3.12 \
    --generate-hashes \
    --output-file requirements.lock
echo "wrote docker/jupyter/requirements.lock"
