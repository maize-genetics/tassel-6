#!/usr/bin/env bash
# Shallow-clone TASSEL 5 source at a pinned tag into reference/tassel5 (git-ignored)
# so agents can read the code being ported. Read-only: never edit or build it here.
#
# Usage: scripts/fetch-reference.sh [tag]
# Default tag matches the golden-data pin in ADR-0009.
set -euo pipefail

TAG="${1:-${TASSEL5_REF:-v5.2.98}}"
REPO="${TASSEL5_REPO:-https://github.com/maize-genetics/tassel.git}"
DEST="$(cd "$(dirname "$0")/.." && pwd)/reference/tassel5"

if [[ -d "$DEST/.git" ]]; then
  if git -C "$DEST" tag --points-at HEAD | grep -qxF "$TAG"; then
    echo "reference/tassel5 already at $TAG"
    exit 0
  fi
  echo "reference/tassel5 is not at $TAG; replacing it"
  rm -rf "$DEST"
fi

mkdir -p "$(dirname "$DEST")"
git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$TAG" "$REPO" "$DEST"
echo "Fetched TASSEL 5 $TAG into reference/tassel5"
