#!/bin/bash
# tests/short.sh — real pipeline on a committed 3 MB CEN6 carve.
# Runs tidehunter → clustering → tarean (no annotation: no library
# committed). Local wall time in tidecluster_1.8.0 env: < 30 s on 2 CPU.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="$ROOT/tests/data/short/CEN6_short.fasta"
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
echo "=== per-TRA consensus wrapper ==="
"$ROOT/tc_per_tra_consensus.py" -p short -c "$NCPU"

PTC="$OUT/short_per_tra_consensus"
[ -s "$PTC/per_tra_consensus.fasta" ] || { echo "FAIL: per_tra_consensus.fasta missing/empty"; exit 1; }
[ -s "$PTC/per_tra_metrics.tsv" ]     || { echo "FAIL: per_tra_metrics.tsv missing/empty"; exit 1; }
[ -s "$PTC/summary.log" ]             || { echo "FAIL: summary.log missing"; exit 1; }
[ -s "$PTC/args.json" ]               || { echo "FAIL: args.json missing"; exit 1; }

# Header has the canonical columns.
HDR=$(head -1 "$PTC/per_tra_metrics.tsv")
for col in id source coverage_frac core_coverage quality_grade flags; do
  echo "$HDR" | grep -qw "$col" || { echo "FAIL: column $col missing from per_tra_metrics.tsv"; exit 1; }
done

# At least one TRA emerges.
ROWS=$(($(wc -l < "$PTC/per_tra_metrics.tsv") - 1))
[ "$ROWS" -ge 1 ] || { echo "FAIL: per_tra_metrics.tsv has zero data rows"; exit 1; }

echo
echo "short PASSED (outputs in $OUT, per-TRA rows: $ROWS)"
