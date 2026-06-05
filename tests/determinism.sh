#!/bin/bash
# tests/determinism.sh — comparative-analysis reproducibility regression
# (issue #4). Runs tc_comparative_analysis.R several times at different
# thread counts into fresh output dirs and asserts the canonicalised
# satellite-family table is byte-identical across all runs. Guards the
# R1-R4 default-path fixes (--max-seqs, deterministic best-hit dedup,
# sorted graph edges, LC_ALL=C) and the --deterministic (R5) m8 ordering.
#
# Needs a multi-sample comparative input (each sample a TideCluster
# output dir with tc_clustering.gff3). Resolution order:
#   1. $TC_COMPARATIVE_FIXTURE (env override)
#   2. tests/data/comparative   (committed fixture, if added)
#   3. test_data/analysis_1.10.5 (bundled local data, untracked)
# If none is present the test SKIPs (exit 0) so CI without the fixture
# does not fail.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/tmp/tests/determinism"
NCPU="${NCPU:-4}"

FIXTURE="${TC_COMPARATIVE_FIXTURE:-}"
if [ -z "$FIXTURE" ]; then
  for cand in "$ROOT/tests/data/comparative" "$ROOT/test_data/analysis_1.10.5"; do
    [ -d "$cand" ] && { FIXTURE="$cand"; break; }
  done
fi
if [ -z "$FIXTURE" ] || [ ! -d "$FIXTURE" ]; then
  echo "SKIP determinism: no comparative fixture found"
  echo "  (set TC_COMPARATIVE_FIXTURE or add tests/data/comparative)"
  exit 0
fi

rm -rf "$OUT"; mkdir -p "$OUT"
echo "=== comparative determinism test ==="
echo "fixture: $FIXTURE"

# Build the input table from sample dirs that look like TC outputs.
TABLE="$OUT/input_table.tsv"
printf "input_dir\tsample_code\ttidecluster_prefix\n" > "$TABLE"
n_samples=0
for d in "$FIXTURE"/*/; do
  d="${d%/}"
  [ -f "$d/tc_clustering.gff3" ] || continue
  printf "%s\t%s\ttc\n" "$d" "$(basename "$d")" >> "$TABLE"
  n_samples=$((n_samples+1))
done
if [ "$n_samples" -lt 2 ]; then
  echo "SKIP determinism: need >=2 sample dirs with tc_clustering.gff3 (found $n_samples)"
  exit 0
fi
echo "samples: $n_samples"

# Canonicalise trc_satellite_families.tsv: drop cosmetic Satellite_family
# / group_id labels, sort the comma-separated TRC lists within every
# remaining cell, then sort rows. Community numbering is order-stable but
# cosmetic, so two structurally identical results canonicalise equal.
canon() {  # $1 = path to trc_satellite_families.tsv ; prints md5
  python3 - "$1" <<'PY'
import sys, hashlib
path = sys.argv[1]
with open(path) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    drop = {i for i, h in enumerate(header) if h in ("Satellite_family", "group_id")}
    rows = []
    for line in fh:
        cells = line.rstrip("\n").split("\t")
        canon_cells = []
        for i, c in enumerate(cells):
            if i in drop:
                continue
            # sort any comma-separated list within the cell
            parts = [p.strip() for p in c.split(",")]
            canon_cells.append(",".join(sorted(parts)))
        rows.append("\t".join(canon_cells))
    rows.sort()
out = "\n".join(rows).encode()
print(hashlib.md5(out).hexdigest())
PY
}

run_one() {  # $1 = label ; $2 = threads ; $3 = extra args
  local odir="$OUT/run_$1"
  rm -rf "$odir"; mkdir -p "$odir"
  "$ROOT/tc_comparative_analysis.R" \
      -i "$TABLE" -o "$odir" -c "$2" $3 \
      > "$odir/log.txt" 2>&1 || {
        echo "FAIL: comparative run '$1' errored; tail of log:"; tail -20 "$odir/log.txt"; exit 1; }
  local fam="$odir/trc_satellite_families.tsv"
  [ -s "$fam" ] || { echo "FAIL: run '$1' produced no trc_satellite_families.tsv"; exit 1; }
  canon "$fam"
}

echo
echo "=== default path: threads {1, $NCPU} x 2 runs (R1-R4 stability) ==="
md5_t1a=$(run_one "t1a" 1 "");      echo "  t1a (threads 1)   $md5_t1a"
md5_t1b=$(run_one "t1b" 1 "");      echo "  t1b (threads 1)   $md5_t1b"
md5_tna=$(run_one "tna" "$NCPU" ""); echo "  tna (threads $NCPU)   $md5_tna"
md5_tnb=$(run_one "tnb" "$NCPU" ""); echo "  tnb (threads $NCPU)   $md5_tnb"

fail=0
for m in "$md5_t1b" "$md5_tna" "$md5_tnb"; do
  [ "$m" = "$md5_t1a" ] || fail=1
done
if [ "$fail" -ne 0 ]; then
  echo "FAIL: satellite-family table differs across runs/thread counts"
  echo "  diff t1a vs tna:"
  diff <(sort "$OUT/run_t1a/trc_satellite_families.tsv") \
       <(sort "$OUT/run_tna/trc_satellite_families.tsv") | head -20 || true
  exit 1
fi
echo "  OK: canonical satellite-family table identical across all 4 runs"

echo
echo "=== --deterministic: byte-identical mmseqs2_results.tsv (R5) ==="
run_one "det_a" "$NCPU" "--deterministic" > /dev/null
run_one "det_b" "$NCPU" "--deterministic" > /dev/null
da=$(md5sum < "$OUT/run_det_a/mmseqs2_results.tsv" | cut -d' ' -f1)
db=$(md5sum < "$OUT/run_det_b/mmseqs2_results.tsv" | cut -d' ' -f1)
echo "  det_a m8 md5 $da"
echo "  det_b m8 md5 $db"
if [ "$da" != "$db" ]; then
  echo "FAIL: --deterministic mmseqs2_results.tsv not byte-identical"
  exit 1
fi
echo "  OK: --deterministic m8 byte-identical"

echo
echo "determinism PASSED (outputs in $OUT)"
