#!/bin/bash
# tests/rerender.sh — report v2 end-to-end.
# Runs the full short-pipeline output through tc_rerender_report.py,
# then walks the generated HTML asserting every <img>/<a> target
# resolves on disk.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FASTA="$ROOT/tests/data/short/CEN6_short.fasta"
WORK="$ROOT/tmp/tests/rerender"
NCPU="${NCPU:-2}"

if [ ! -s "$FASTA" ]; then
  echo "FAIL: $FASTA missing"; exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"
cd "$WORK"

echo "=== run_all on short fixture ==="
"$ROOT/TideCluster.py" run_all -c "$NCPU" -pr rerender -f "$FASTA"

echo
echo "=== rerender report v2 ==="
"$ROOT/tc_rerender_report.py" --input-dir "$WORK" --prefix rerender

REPORT="$WORK/rerender_report_v2"
[ -s "$REPORT/data/report.json" ]    || { echo "FAIL: report.json missing"; exit 1; }
[ -s "$REPORT/index.html" ]          || { echo "FAIL: index.html missing"; exit 1; }
[ -s "$REPORT/tarean.html" ]         || { echo "FAIL: tarean.html missing"; exit 1; }
[ -s "$REPORT/kite.html" ]           || { echo "FAIL: kite.html missing"; exit 1; }
[ -s "$REPORT/superfamilies.html" ]  || { echo "FAIL: superfamilies.html missing"; exit 1; }
[ -d "$REPORT/assets/datatables" ]   || { echo "FAIL: datatables assets not copied"; exit 1; }
[ -d "$REPORT/trc" ]                 || { echo "FAIL: trc/ dir missing"; exit 1; }

echo
echo "=== check every src/href resolves ==="
python3 - <<PYEOF
import os, re, sys
OUT = "$REPORT"
pages = [os.path.join(OUT, f) for f in os.listdir(OUT) if f.endswith(".html")]
pages += [os.path.join(OUT, "trc", f) for f in os.listdir(os.path.join(OUT, "trc")) if f.endswith(".html")]
bad = 0
for p in pages:
    src = open(p).read()
    refs = re.findall(r'(?:src|href)="([^"]+)"', src)
    for r in refs:
        if r.startswith(("http://","https://","data:","#","mailto:")): continue
        target = r.split("#", 1)[0]
        if not target: continue
        resolved = os.path.normpath(os.path.join(os.path.dirname(p), target))
        if not os.path.exists(resolved):
            print(f"FAIL: {p} -> missing {r}", file=sys.stderr)
            bad += 1
if bad:
    print(f"FAIL: {bad} unresolvable references", file=sys.stderr); sys.exit(1)
print(f"OK: {len(pages)} pages, all src/href resolve")
PYEOF

echo
echo "rerender PASSED"
