#!/bin/bash
# tests/short.sh — real pipeline on committed 18 MB CEN6 part fasta.
# Runs tidehunter → clustering → tarean (no annotation: no library
# committed). Expected wall time ~ 1-3 min on 2 CPUs.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="$ROOT/test_data/CEN6_ver_220406_part.fasta"
OUT="$ROOT/tmp/tests/short"
NCPU="${NCPU:-2}"

if [ ! -s "$FASTA" ]; then
  echo "FAIL: $FASTA missing"; exit 1
fi

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"

echo "=== run_all on CEN6 part ==="
"$ROOT/TideCluster.py" run_all \
  -c "$NCPU" -pr short -f "$FASTA"

[ -s "$OUT/short_tidehunter.gff3" ] || { echo "FAIL: no tidehunter.gff3"; exit 1; }
[ -s "$OUT/short_clustering.gff3" ] || { echo "FAIL: no clustering.gff3"; exit 1; }
[ -s "$OUT/short_index.html" ]      || { echo "FAIL: no index.html"; exit 1; }

echo
echo "short PASSED (outputs in $OUT)"
