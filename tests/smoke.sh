#!/bin/bash
# tests/smoke.sh — fast sanity check (~15 s).
# Probes every CLI's --help, then runs tidehunter + clustering on a
# 1 MB repeat-rich carve of CEN6 committed to tests/data/smoke/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/tests/data/smoke"
OUT="$ROOT/tmp/tests/smoke"
NCPU="${NCPU:-2}"

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"

echo "=== CLI --version / --help ==="
"$ROOT/TideCluster.py" --version
"$ROOT/TideCluster.py" --help > /dev/null
"$ROOT/tc_update_gff3.py" --help > /dev/null
"$ROOT/tc_reannotate.py" --help > /dev/null
"$ROOT/tc_merge_annotations.py" --help > /dev/null
"$ROOT/tc_comparative_analysis.R" --help > /dev/null
"$ROOT/tc_summarize_comparative_analysis.R" --help > /dev/null

echo
echo "=== tidehunter on smoke fasta ==="
"$ROOT/TideCluster.py" tidehunter \
  -c "$NCPU" -pr smoke -f "$DATA/CEN6_smoke.fasta"

[ -s "$OUT/smoke_tidehunter.gff3" ] || {
  echo "FAIL: smoke_tidehunter.gff3 missing or empty"; exit 1; }

echo
echo "=== clustering (min_length=500 to allow tiny arrays) ==="
"$ROOT/TideCluster.py" clustering \
  -c "$NCPU" -pr smoke -f "$DATA/CEN6_smoke.fasta" -m 500

[ -s "$OUT/smoke_clustering.gff3" ] || {
  echo "FAIL: smoke_clustering.gff3 missing or empty"; exit 1; }

echo
echo "smoke PASSED (outputs in $OUT)"
