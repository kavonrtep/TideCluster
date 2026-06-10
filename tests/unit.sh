#!/bin/bash
# tests/unit.sh — fast pure-python unit tests (no external tools / data).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== unit: annotation coverage ==="
python3 "$ROOT/tests/test_annotation_coverage.py"

echo "=== unit: cluster rescue (Pass 5) ==="
python3 "$ROOT/tests/test_cluster_rescue.py"

echo "=== unit: strongest = argmax(id_med) ==="
python3 "$ROOT/tests/test_strongest_by_identity.py"

echo "=== unit: HOR-order confidence tier ==="
python3 "$ROOT/tests/test_hor_order_confidence.py"

echo "=== unit: report TSV bundle + lazy child rows ==="
python3 "$ROOT/tests/test_report_exports.py"

echo "unit tests OK"
