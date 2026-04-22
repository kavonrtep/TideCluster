#!/bin/bash
# tests/rerender.sh — report v1 AND v2 end-to-end.
# Runs the full short-pipeline (which now builds v2 automatically at
# the end of tarean), then walks every generated HTML in both the
# v1 root and the v2 subtree asserting every <img>/<a> target
# resolves on disk. Tests both reports to catch regressions in either.
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

echo "=== run_all on short fixture (builds v1 + v2) ==="
"$ROOT/TideCluster.py" run_all -c "$NCPU" -pr rerender -f "$FASTA"

REPORT="$WORK/rerender_report_v2"

# v1 legacy files always present
[ -s "$WORK/rerender_index.html" ]             || { echo "FAIL: v1 index.html missing"; exit 1; }
[ -s "$WORK/rerender_tarean_report.html" ]     || { echo "FAIL: v1 tarean_report.html missing"; exit 1; }
[ -s "$WORK/rerender_trc_superfamilies.html" ] || { echo "FAIL: v1 trc_superfamilies.html missing"; exit 1; }
[ -s "$WORK/rerender_pipeline_stats.json" ]    || { echo "FAIL: pipeline_stats.json missing"; exit 1; }
[ -s "$WORK/rerender_seqid_lengths.tsv" ]      || { echo "FAIL: seqid_lengths.tsv missing"; exit 1; }

# v2 files
[ -s "$REPORT/data/report.json" ]    || { echo "FAIL: v2 report.json missing"; exit 1; }
[ -s "$REPORT/index.html" ]          || { echo "FAIL: v2 index.html missing"; exit 1; }
[ -s "$REPORT/tarean.html" ]         || { echo "FAIL: v2 tarean.html missing"; exit 1; }
[ -s "$REPORT/kite.html" ]           || { echo "FAIL: v2 kite.html missing"; exit 1; }
[ -s "$REPORT/superfamilies.html" ]  || { echo "FAIL: v2 superfamilies.html missing"; exit 1; }
[ -d "$REPORT/assets/datatables" ]   || { echo "FAIL: v2 datatables assets not copied"; exit 1; }
[ -d "$REPORT/trc" ]                 || { echo "FAIL: v2 trc/ dir missing"; exit 1; }

echo
echo "=== check every src/href in v1 + v2 HTML resolves ==="
python3 - <<PYEOF
import os, re, sys
WORK = "$WORK"
REPORT = "$REPORT"

# Collect v1 HTML in WORK root + per-TRC TAREAN + per-TRC KITE subdirs.
v1_pages = []
for f in os.listdir(WORK):
    full = os.path.join(WORK, f)
    if os.path.isfile(full) and f.endswith(".html"):
        v1_pages.append(full)
for sub in ("rerender_tarean", "rerender_kite"):
    p = os.path.join(WORK, sub)
    if not os.path.isdir(p): continue
    for root, _, files in os.walk(p):
        for f in files:
            if f.endswith(".html"):
                v1_pages.append(os.path.join(root, f))

# Collect v2 HTML.
v2_pages = [os.path.join(REPORT, f)
            for f in os.listdir(REPORT) if f.endswith(".html")]
v2_pages += [os.path.join(REPORT, "trc", f)
             for f in os.listdir(os.path.join(REPORT, "trc"))
             if f.endswith(".html")]

def walk_refs(pages, label):
    bad = 0
    for p in pages:
        try: src = open(p).read()
        except Exception: continue
        if not src.strip(): continue            # stub pages (e.g. v1 empty reports)
        for r in re.findall(r'(?:src|href)="([^"]+)"', src):
            if r.startswith(("http://","https://","data:","#","mailto:")): continue
            target = r.split("#", 1)[0]
            if not target: continue
            resolved = os.path.normpath(os.path.join(os.path.dirname(p), target))
            if not os.path.exists(resolved):
                print(f"FAIL [{label}]: {p} -> missing {r}", file=sys.stderr)
                bad += 1
    return bad

bad_v1 = walk_refs(v1_pages, "v1")
bad_v2 = walk_refs(v2_pages, "v2")
print(f"v1: {len(v1_pages)} pages, {bad_v1} unresolvable")
print(f"v2: {len(v2_pages)} pages, {bad_v2} unresolvable")
if bad_v1 or bad_v2: sys.exit(1)

# v2 must link back to v1 index (the new "Legacy report" nav item).
v2_idx = open(os.path.join(REPORT, "index.html")).read()
assert 'rerender_index.html' in v2_idx, "v2 index is missing Legacy report link"
for f in os.listdir(os.path.join(REPORT, "trc")):
    if not f.endswith(".html"): continue
    src = open(os.path.join(REPORT, "trc", f)).read()
    assert 'rerender_index.html' in src, f"v2 {f} missing Legacy report link"
print("OK: Legacy report link present on every v2 page")
PYEOF

echo
echo "rerender PASSED"
