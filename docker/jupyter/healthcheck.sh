#!/usr/bin/env bash
set -euo pipefail

# Lightweight container health probe.
# Success criteria:
#  1. Python venv present and can import tensorflow symbolically.
#  2. R can load Rsafd package.
#  3. reticulate sees python tensorflow module.

fail() { echo "HEALTH FAIL: $*" >&2; exit 1; }

PYTHON_BIN="/opt/venv/bin/python"
[[ -x "$PYTHON_BIN" ]] || fail "venv python missing"

# Quick python import (suppress verbose TF logs)
TF_VERSION=$($PYTHON_BIN - <<'PY' 2>/dev/null || true
import os; os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
try:
	import tensorflow as tf
	print(tf.__version__)
except Exception as e:
	pass
PY
)
[[ -n "$TF_VERSION" ]] || fail "tensorflow import failed"

# R package checks (quiet)
R -q -e "suppressPackageStartupMessages(library(Rsafd)); cat('Rsafd OK\n')" >/dev/null 2>&1 || fail "Rsafd load failed"
R -q -e "suppressPackageStartupMessages(library(reticulate)); py <- import('tensorflow'); cat('reticulate OK\n')" >/dev/null 2>&1 || fail "reticulate tensorflow import failed"

echo "healthy: tf=${TF_VERSION}"
