#!/usr/bin/env bash
# Entry point of the check image. The package source is at /work/Rsafd
# (baked in at build time, or bind-mounted read-only by docker-compose).
set -euo pipefail

SRC=/work/Rsafd
OUT="${RSAFD_OUT:-/tmp/rsafd-out}"
mode="${1:-check}"

stage() {
  # Never build inside the (possibly read-only) source tree.
  rm -rf "$OUT" && mkdir -p "$OUT"
  cp -R "$SRC" "$OUT/Rsafd"
  rm -rf "$OUT/Rsafd/.git"
}

case "$mode" in
  check)
    stage
    cd "$OUT"
    R CMD build --no-manual Rsafd
    tarball=$(ls Rsafd_*.tar.gz)
    # Xvfb gives tcltk (loaded via plot3D -> misc3d) a display, as on CRAN's Linux
    # check machines. The container runs offline, so skip network-only checks.
    _R_CHECK_CRAN_INCOMING_REMOTE_=false _R_CHECK_SYSTEM_CLOCK_=false \
      xvfb-run -a R CMD check --as-cran --no-manual "$tarball" || {
      cat Rsafd.Rcheck/00check.log; exit 1; }
    # R CMD check exits 0 on WARNINGs; treat them (and ERRORs) as failures.
    if grep -E '^Status: .*(WARNING|ERROR)' Rsafd.Rcheck/00check.log; then
      echo "R CMD check reported warnings or errors" >&2
      exit 1
    fi
    grep '^Status:' Rsafd.Rcheck/00check.log
    ;;
  test)
    stage
    mkdir -p "$OUT/lib"
    R CMD INSTALL --library="$OUT/lib" "$OUT/Rsafd"
    cd "$OUT/Rsafd/tests"
    R_LIBS="$OUT/lib" xvfb-run -a Rscript -e 'res <- testthat::test_dir("testthat", package = "Rsafd", load_package = "installed", reporter = "summary", stop_on_failure = TRUE); df <- as.data.frame(res); if (sum(df$nb) == 0) stop("no expectations ran")'
    ;;
  lint)
    # object_usage_linter resolves symbols against the installed namespace.
    stage
    mkdir -p "$OUT/lib"
    R CMD INSTALL --library="$OUT/lib" "$OUT/Rsafd" >/dev/null
    R_LIBS="$OUT/lib" Rscript -e 'l <- lintr::lint_package("/work/Rsafd"); print(l); if (length(l)) quit(status = 1); cat("no lints\n")'
    ;;
  *)
    echo "usage: rsafd-check [check|test|lint]" >&2
    exit 2
    ;;
esac
