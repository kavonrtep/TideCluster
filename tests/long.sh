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

# Report link-closure: every local ref reachable from the index must resolve,
# both as-built and after a simulated 'maximal' cleanup (kite/tarean/dotplots
# deleted) — the self-contained-report contract downstream tools rely on.
echo
echo "=== report link-closure (intact) ==="
python3 "$ROOT/tests/report_linkcheck.py" "$OUT" long
echo "=== report link-closure (after simulated 'maximal' purge) ==="
python3 "$ROOT/tests/report_linkcheck.py" --purge "$OUT" long
# Optional defense-in-depth: a second, independent opinion from lychee if it is
# installed (handles HTML corners the stdlib walker doesn't). Skipped otherwise.
if command -v lychee >/dev/null 2>&1; then
  echo "=== lychee --offline (defense-in-depth) ==="
  lychee --offline --no-progress --include-fragments "$OUT/long_index.html" \
    || { echo "FAIL: lychee reported broken links"; exit 1; }
else
  echo "(lychee not installed; skipping optional external link check — install to enable)"
fi

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

# ---------------------------------------------------------------------------
# run_all --cleanup: the four guarantees, end to end on a real run.
# A separate run with a low -M so TAREAN actually produces per-TRC directories
# (the default -M leaves this fixture with zero TAREAN output, so the run above
# exercises none of the tarean purge patterns).
# ---------------------------------------------------------------------------
COUT="$ROOT/tmp/tests/long_cleanup"
LIB="$ROOT/test_data/solanum_custom_library_v1_RM_formated.fasta"
rm -rf "$COUT"; mkdir -p "$COUT"
echo
# -M 1000 so TAREAN produces per-TRC directories (at the default -M this fixture
# yields none, exercising no tarean purge pattern), and -l so the ANNOTATION step
# runs: consensus_sequences_all.fasta and the RepeatMasker leftovers that P5
# targets only exist on an annotated run. The library need not match the genome —
# the files are written either way.
echo "=== run_all --cleanup (-M 1000 -l <library> so TAREAN and annotation run) ==="
if [ -s "$LIB" ]; then LIB_ARG="-l $LIB"; else LIB_ARG=""; echo "  (no bundled library; G2 will be partial)"; fi
( cd "$COUT" && "$ROOT/TideCluster.py" run_all -c "$NCPU" -pr clean -M 1000 \
    --cleanup $LIB_ARG -f "$FASTA" > run.log 2>&1 ) || { echo "FAIL: run_all --cleanup exited non-zero"; tail -20 "$COUT/run.log"; exit 1; }

grep -q "^cleanup: removed" "$COUT/run.log" || { echo "FAIL: no cleanup summary line in the run log"; exit 1; }
grep "^cleanup: removed" "$COUT/run.log" | sed 's/^/  /' | cut -c1-100

# The purge actually happened.
[ -d "$COUT/clean_tarean" ] || { echo "FAIL: clean_tarean/ missing entirely (tree-level delete?)"; exit 1; }
for gone in "clean_kite/kitehor.periodogram" "clean_clustering.gff3_1.gff3"; do
  [ -e "$COUT/$gone" ] && { echo "FAIL: $gone survived cleanup"; exit 1; }
done
find "$COUT" -name "*.kmers" -o -name "ggmin.RData" -o -name "monomers.RData" | grep -q . \
  && { echo "FAIL: TAREAN scratch survived cleanup"; exit 1; }

# G1 — the report still resolves, intact and under a further 'maximal' purge.
echo "=== G1: report link-closure after --cleanup ==="
python3 "$ROOT/tests/report_linkcheck.py" "$COUT" clean
python3 "$ROOT/tests/report_linkcheck.py" --purge "$COUT" clean

# G2 — comparative-analysis inputs (the six documented in the README).
echo "=== G2: comparative-analysis inputs after --cleanup ==="
for f in clean_consensus_dimer_library.fasta clean_clustering.gff3; do
  [ -s "$COUT/$f" ] || { echo "FAIL: comparative input $f missing/empty after cleanup"; exit 1; }
done
# SSRS_summary.csv is legitimately EMPTY when a run finds no SSRs (as this
# fixture does), so presence is the contract, not size.
[ -f "$COUT/clean_tarean/SSRS_summary.csv" ] \
  || { echo "FAIL: comparative input clean_tarean/SSRS_summary.csv missing after cleanup"; exit 1; }
if [ -n "$LIB_ARG" ]; then
  # The one comparative input that lives INSIDE a purged directory: P5 operates
  # in clean_consensus/, so this is the file the purge could plausibly eat.
  # (Written by the ANNOTATION step, hence the -l above.) clean_annotation.gff3
  # / .tsv are top-level files no pattern targets, and whether they are written
  # at all depends on the library matching the genome, so they are not asserted.
  [ -s "$COUT/clean_consensus/consensus_sequences_all.fasta" ] \
    || { echo "FAIL: comparative input consensus_sequences_all.fasta missing/empty after cleanup"; exit 1; }
  echo "  ok (incl. consensus_sequences_all.fasta, the one inside a purged directory)"
else
  echo "  ok (consensus_sequences_all.fasta skipped: no library available, so the annotation step did not run)"
fi

# G3 — the report can still be regenerated from what survived.
echo "=== G3: re-render after --cleanup ==="
( cd "$COUT" && python3 "$ROOT/tc_rerender_report.py" --input-dir . --prefix clean > rerender.log 2>&1 ) \
  || { echo "FAIL: re-render after cleanup exited non-zero"; tail -20 "$COUT/rerender.log"; exit 1; }
python3 "$ROOT/tests/report_linkcheck.py" "$COUT" clean

# G4 — per-TRA consensus inputs, and the wrapper itself.
echo "=== G4: per-TRA consensus after --cleanup ==="
for f in clean_kite/monomer_size_top3_estimats.csv clean_tidehunter.gff3 clean_clustering.gff3; do
  [ -s "$COUT/$f" ] || { echo "FAIL: per-TRA input $f missing/empty after cleanup"; exit 1; }
done
[ -d "$COUT/clean_tarean/fasta" ] || { echo "FAIL: clean_tarean/fasta/ deleted — the glob trap"; exit 1; }
( cd "$COUT" && "$ROOT/tc_per_tra_consensus.py" -p clean -c "$NCPU" > ptc.log 2>&1 ) \
  || { echo "FAIL: tc_per_tra_consensus.py failed after cleanup"; tail -20 "$COUT/ptc.log"; exit 1; }
[ -s "$COUT/clean_per_tra_consensus/per_tra_consensus.fasta" ] \
  || { echo "FAIL: per_tra_consensus.fasta missing after cleanup"; exit 1; }

echo
echo "long PASSED (clusters: $NCLUST, per-TRA rows: $ROWS, grade A: $NA, outputs in $OUT)"
echo "  --cleanup guarantees G1-G4 verified in $COUT"
