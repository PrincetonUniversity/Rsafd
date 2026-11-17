#!/usr/bin/env bash
# Cut a release: set DESCRIPTION Version to X.Y.Z, commit, tag vX.Y.Z, push.
# Pushing the tag triggers .github/workflows/release.yml, which runs the full
# check matrix before publishing anything.
#
#   scripts/release.sh 1.1.0 --dry-run   # show what would happen
#   scripts/release.sh 1.1.0             # do it
#
# Safe to re-run: steps that are already done are skipped.
# Works under bash and zsh (no bash-4-only features).
set -euo pipefail

usage() { echo "usage: $0 X.Y.Z [--dry-run]" >&2; exit 2; }

version=""
dry_run=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) dry_run=1 ;;
    -h|--help) usage ;;
    *) [ -z "$version" ] && version="$arg" || usage ;;
  esac
done
[ -n "$version" ] || usage

if ! printf '%s' "$version" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
  echo "error: version '$version' is not X.Y.Z" >&2
  exit 1
fi
tag="v$version"

run() {
  if [ "$dry_run" -eq 1 ]; then
    echo "[dry-run] $*"
  else
    echo "+ $*"
    "$@"
  fi
}

cd "$(git rev-parse --show-toplevel)"

branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" != "main" ]; then
  echo "error: releases are cut from main (current branch: $branch)" >&2
  exit 1
fi
if [ -n "$(git status --porcelain)" ]; then
  echo "error: working tree is not clean" >&2
  exit 1
fi
git fetch --quiet origin main --tags
if [ "$(git rev-parse HEAD)" != "$(git rev-parse origin/main)" ]; then
  echo "error: local main is not identical to origin/main; pull or push first" >&2
  exit 1
fi

if ! grep -q "^# Rsafd ${version}\$" NEWS.md; then
  echo "error: NEWS.md needs a '# Rsafd ${version}' section describing the release" >&2
  exit 1
fi

current=$(sed -n 's/^Version:[[:space:]]*//p' DESCRIPTION)
if [ "$current" != "$version" ]; then
  if [ "$dry_run" -eq 1 ]; then
    echo "[dry-run] set DESCRIPTION Version: $current -> $version"
  else
    sed -i.bak "s/^Version:.*/Version: ${version}/" DESCRIPTION && rm -f DESCRIPTION.bak
  fi
  run git commit -m "Release ${version}" -- DESCRIPTION
else
  echo "DESCRIPTION already at ${version}"
fi

if git rev-parse -q --verify "refs/tags/${tag}" >/dev/null; then
  if [ "$(git rev-list -n 1 "$tag")" != "$(git rev-parse HEAD)" ] && [ "$dry_run" -eq 0 ]; then
    echo "error: tag ${tag} already exists on a different commit" >&2
    exit 1
  fi
  echo "tag ${tag} already exists"
else
  run git tag -a "$tag" -m "Rsafd ${version}"
fi

run git push origin main
run git push origin "$tag"
echo "done: the release workflow will publish ${tag} once all checks pass"
