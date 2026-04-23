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

echo "=== run_all on short fixture (builds v1-legacy + v2) ==="
"$ROOT/TideCluster.py" run_all -c "$NCPU" -pr rerender -f "$FASTA"

LEGACY="$WORK/rerender_report_legacy"
REPORT="$WORK/rerender_report"

# v1 legacy files moved into the legacy dir
[ -s "$LEGACY/rerender_index.html" ]             || { echo "FAIL: v1 index.html missing from legacy"; exit 1; }
[ -s "$LEGACY/rerender_tarean_report.html" ]     || { echo "FAIL: v1 tarean_report.html missing from legacy"; exit 1; }
[ -s "$LEGACY/rerender_trc_superfamilies.html" ] || { echo "FAIL: v1 trc_superfamilies.html missing from legacy"; exit 1; }
# pipeline side-cars stay at root (data files, not HTML reports)
[ -s "$WORK/rerender_pipeline_stats.json" ]      || { echo "FAIL: pipeline_stats.json missing"; exit 1; }
[ -s "$WORK/rerender_seqid_lengths.tsv" ]        || { echo "FAIL: seqid_lengths.tsv missing"; exit 1; }
# The three v1 sub-reports must not remain at the root — they moved into
# the legacy dir. (rerender_index.html at the root is the v2 landing;
# its presence is expected.)
for f in "$WORK"/rerender_tarean_report.html \
         "$WORK"/rerender_kite_report.html \
         "$WORK"/rerender_trc_superfamilies.html; do
  [ -f "$f" ] && { echo "FAIL: stale v1 $f at root after legacy move"; exit 1; } || true
done

# v2 files
[ -s "$REPORT/data/report.json" ]    || { echo "FAIL: v2 report.json missing"; exit 1; }
[ -s "$WORK/rerender_index.html" ]   || { echo "FAIL: v2 rerender_index.html missing at root"; exit 1; }
[ -s "$REPORT/tarean.html" ]         || { echo "FAIL: v2 tarean.html missing"; exit 1; }
[ -s "$REPORT/kite.html" ]           || { echo "FAIL: v2 kite.html missing"; exit 1; }
[ -s "$REPORT/superfamilies.html" ]  || { echo "FAIL: v2 superfamilies.html missing"; exit 1; }
[ -d "$REPORT/assets/datatables" ]   || { echo "FAIL: v2 datatables assets not copied"; exit 1; }
[ -d "$REPORT/trc" ]                 || { echo "FAIL: v2 trc/ dir missing"; exit 1; }
# Stale intermediate dir (1.9.0-dev) must be cleaned up
[ -d "$WORK/rerender_report_v2" ] && { echo "FAIL: stale rerender_report_v2/ not cleaned up"; exit 1; } || true

echo
echo "=== check every src/href in v1-legacy + v2 HTML resolves ==="
python3 - <<PYEOF
import os, re, sys
WORK   = "$WORK"
LEGACY = "$LEGACY"
REPORT = "$REPORT"

# v1 HTML now in rerender_report_legacy/ + per-TRC HTML still in data dirs.
v1_pages = [os.path.join(LEGACY, f)
            for f in os.listdir(LEGACY) if f.endswith(".html")]
for sub in ("rerender_tarean", "rerender_kite"):
    p = os.path.join(WORK, sub)
    if not os.path.isdir(p): continue
    for root, _, files in os.walk(p):
        for f in files:
            if f.endswith(".html"):
                v1_pages.append(os.path.join(root, f))

# v2 HTML: root index + subtree + dashboards.
v2_pages = [os.path.join(WORK, "rerender_index.html")]
v2_pages += [os.path.join(REPORT, f)
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

bad_v1 = walk_refs(v1_pages, "v1-legacy")
bad_v2 = walk_refs(v2_pages, "v2")
print(f"v1-legacy: {len(v1_pages)} pages, {bad_v1} unresolvable")
print(f"v2:        {len(v2_pages)} pages, {bad_v2} unresolvable")
if bad_v1 or bad_v2: sys.exit(1)

# v2 must link back to v1 legacy index (the "Legacy report" nav item).
expected_legacy = "rerender_report_legacy/rerender_index.html"
def has_legacy(p, depth_prefix):
    src = open(p).read()
    return f'{depth_prefix}{expected_legacy}' in src
assert has_legacy(os.path.join(WORK, "rerender_index.html"), ""), \
    "root v2 index is missing Legacy report link"
assert has_legacy(os.path.join(REPORT, "tarean.html"), "../"), \
    "subtree v2 tarean.html is missing Legacy report link"
for f in os.listdir(os.path.join(REPORT, "trc")):
    if not f.endswith(".html"): continue
    assert has_legacy(os.path.join(REPORT, "trc", f), "../../"), \
        f"dashboard {f} missing Legacy report link"
print("OK: Legacy report link present and correct depth on every v2 page")
PYEOF

echo
echo "rerender PASSED"
