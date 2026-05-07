#!/bin/bash
# tests/long.sh — release gate. Runs the full pipeline on a committed
# 10 MB CEN6 carve. Override via LONG_FASTA env var (e.g. point it at
# the full 180 MB test_data/CEN6_ver_220406.fasta) for deep local runs.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="${LONG_FASTA:-$ROOT/tests/data/long/CEN6_long.fasta}"
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
echo "=== per-TRA consensus wrapper ==="
"$ROOT/tc_per_tra_consensus.py" -p long -c "$NCPU"

PTC="$OUT/long_per_tra_consensus"
[ -s "$PTC/per_tra_consensus.fasta" ] || { echo "FAIL: per_tra_consensus.fasta missing/empty"; exit 1; }
[ -s "$PTC/per_tra_metrics.tsv" ]     || { echo "FAIL: per_tra_metrics.tsv missing/empty"; exit 1; }
[ -s "$PTC/summary.log" ]             || { echo "FAIL: summary.log missing"; exit 1; }
[ -s "$PTC/args.json" ]               || { echo "FAIL: args.json missing"; exit 1; }

HDR=$(head -1 "$PTC/per_tra_metrics.tsv")
for col in id source coverage_frac core_coverage quality_grade flags; do
  echo "$HDR" | grep -qw "$col" || { echo "FAIL: column $col missing from per_tra_metrics.tsv"; exit 1; }
done

ROWS=$(($(wc -l < "$PTC/per_tra_metrics.tsv") - 1))
[ "$ROWS" -ge 1 ] || { echo "FAIL: per_tra_metrics.tsv has zero data rows"; exit 1; }

# Grade-A count by header-driven column lookup (robust to column reordering).
NA=$(awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i=="quality_grade") col=i; next }
                 NR>1  && $col=="A"' "$PTC/per_tra_metrics.tsv" | wc -l)

echo
echo "long PASSED (clusters: $NCLUST, per-TRA rows: $ROWS, grade A: $NA, outputs in $OUT)"
