#!/bin/bash
# tests/unit.sh — fast pure-python unit tests (no external tools / data).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== unit: annotation coverage ==="
python3 "$ROOT/tests/test_annotation_coverage.py"

echo "=== unit: cluster rescue (Pass 5) ==="
python3 "$ROOT/tests/test_cluster_rescue.py"

echo "unit tests OK"
