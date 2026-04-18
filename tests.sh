#!/bin/bash
# tests.sh — dispatcher. ./tests.sh {smoke|short|long|all}
# Backwards compat: ./tests.sh <N> (integer) runs the long test with N CPU.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

LEVEL="${1:-smoke}"
export NCPU="${NCPU:-${2:-2}}"

if [[ "$LEVEL" =~ ^[0-9]+$ ]]; then
  export NCPU="$LEVEL"
  LEVEL="long"
fi

if command -v conda >/dev/null 2>&1 && [ -z "${CONDA_DEFAULT_ENV:-}" ]; then
  eval "$(conda shell.bash hook)"
  if conda env list | grep -q '^tidecluster '; then
    conda activate tidecluster
  fi
fi

case "$LEVEL" in
  smoke) bash "$ROOT/tests/smoke.sh" ;;
  short) bash "$ROOT/tests/short.sh" ;;
  long)  bash "$ROOT/tests/long.sh"  ;;
  all)
    bash "$ROOT/tests/smoke.sh"
    bash "$ROOT/tests/short.sh"
    bash "$ROOT/tests/long.sh"
    ;;
  *)
    echo "usage: $0 {smoke|short|long|all|<NCPU>}" >&2
    exit 2
    ;;
esac
