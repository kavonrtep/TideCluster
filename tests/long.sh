#!/bin/bash
# tests/long.sh — release gate. Runs the full pipeline on the 180 MB
# CEN6 fasta (committed under test_data/). Override via LONG_FASTA
# env var to point at something larger for deep local runs.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="${LONG_FASTA:-$ROOT/test_data/CEN6_ver_220406.fasta}"
OUT="$ROOT/tmp/tests/long"
NCPU="${NCPU:-2}"

if [ ! -s "$FASTA" ]; then
  echo "FAIL: $FASTA missing (set LONG_FASTA to override)"; exit 1
fi

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"

echo "=== run_all on $FASTA ==="
"$ROOT/TideCluster.py" run_all \
  -c "$NCPU" -pr long -f "$FASTA"

# Loose range assertions: tolerate cross-version drift (playbook §3.7).
[ -s "$OUT/long_tidehunter.gff3" ] || { echo "FAIL: no tidehunter.gff3"; exit 1; }
[ -s "$OUT/long_clustering.gff3" ] || { echo "FAIL: no clustering.gff3"; exit 1; }
[ -s "$OUT/long_index.html" ]      || { echo "FAIL: no index.html"; exit 1; }

NCLUST=$(grep -c 'Name=TRC_' "$OUT/long_clustering.gff3" || true)
[ "${NCLUST:-0}" -ge 1 ] || { echo "FAIL: no TRC_ clusters in long_clustering.gff3"; exit 1; }

echo
echo "long PASSED (clusters found: $NCLUST, outputs in $OUT)"
