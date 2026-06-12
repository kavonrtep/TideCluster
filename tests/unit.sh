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

echo "=== unit: SSR raw/consensus + founder (real kite fixture) ==="
python3 "$ROOT/tests/test_ssr_raw_fixture.py"

echo "=== unit: harmonic-ladder founder (Pass 7) ==="
python3 "$ROOT/tests/test_harmonic_ladder.py"

echo "=== unit: RepeatMasker empty-.out guard ==="
python3 "$ROOT/tests/test_repeatmasker_empty_out.py"

echo "=== unit: dominant-score ladder founder (Pass 7b) ==="
python3 "$ROOT/tests/test_dominant_ladder.py"

echo "unit tests OK"
