#!/usr/bin/env python3
"""Link-closure checker for the generated TideCluster v2 HTML report.

Starting from <prefix>_index.html it walks the whole served link-closure — v2
subpages, per-TRC dashboards, the vendored TAREAN drill-down report.html pages,
and the legacy v1 tree — following every local <a>/<img> link and asserting each
referenced local asset (PNG, CSV, HTML, …) exists. External URLs (http/https,
mailto, data:, #-only anchors) are skipped: this is an OFFLINE structural check.

It parses with the stdlib html.parser, so it sees single-quoted, double-quoted,
and unquoted attributes alike — unlike a naive `href="…"` regex, which silently
misses `href='…'` (the class of bug that shipped broken legacy drill-down links).

--purge runs the same check on a throwaway copy with the <prefix>_kite/,
<prefix>_tarean/ and dotplots/ scratch trees deleted, i.e. exactly what a
downstream 'maximal' cleanup (CARP) removes — so it proves the report stays
self-contained after those trees are gone.

Usage:
  python3 tests/report_linkcheck.py [--purge] <run_dir> <prefix>
Exit status is non-zero if any local reference is broken.
"""
import argparse
import os
import shutil
import sys
import tempfile
from html.parser import HTMLParser
from urllib.parse import urldefrag

# Scratch trees a 'maximal' cleanup deletes (relative to the run dir; dotplots
# is unprefixed). <prefix> is substituted at runtime.
PURGE_TREES = ["{prefix}_kite", "{prefix}_tarean", "dotplots"]


class _RefCollector(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.refs = []

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in ("src", "href") and value:
                self.refs.append(value)


def _local_refs(html_path):
    """Local (non-external) src/href targets on a page, fragments/queries stripped."""
    p = _RefCollector()
    with open(html_path, encoding="utf-8", errors="replace") as fh:
        p.feed(fh.read())
    out = []
    for raw in p.refs:
        low = raw.strip().lower()
        if low.startswith(("http://", "https://", "mailto:", "data:", "javascript:")):
            continue
        # drop #fragment and ?query, leaving the file path
        u = urldefrag(raw)[0].split("?", 1)[0]
        if u:
            out.append(u)
    return out


def check_report(run_dir, prefix):
    """Walk the closure from <prefix>_index.html. Returns a list of broken refs
    as (page, ref) tuples (empty when the report is fully self-contained)."""
    start = os.path.normpath(os.path.join(run_dir, f"{prefix}_index.html"))
    broken, seen, pages, assets = [], set(), 0, 0
    stack = [(start, "<entry point>")]
    while stack:
        path, via = stack.pop()
        if path in seen:
            continue
        seen.add(path)
        if not os.path.isfile(path):
            broken.append((via, os.path.relpath(path, run_dir)))
            continue
        if not path.endswith(".html"):
            continue
        pages += 1
        base = os.path.dirname(path)
        for ref in _local_refs(path):
            target = os.path.normpath(os.path.join(base, ref))
            if ref.endswith(".html"):
                stack.append((target, os.path.relpath(path, run_dir)))
            else:
                assets += 1
                if not os.path.isfile(target):
                    broken.append((os.path.relpath(path, run_dir), ref))
    return broken, pages, assets


def _purged_copy(run_dir, prefix):
    tmp = tempfile.mkdtemp(prefix="tc_linkcheck_")
    dst = os.path.join(tmp, "run")
    shutil.copytree(run_dir, dst)
    for tree in PURGE_TREES:
        shutil.rmtree(os.path.join(dst, tree.format(prefix=prefix)), ignore_errors=True)
    return tmp, dst


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", help="run directory holding <prefix>_index.html")
    ap.add_argument("prefix", help="report prefix (e.g. 'long', 'tc')")
    ap.add_argument("--purge", action="store_true",
                    help="check a copy with the kite/tarean/dotplots trees deleted")
    args = ap.parse_args(argv)

    run_dir, cleanup = os.path.abspath(args.run_dir), None
    label = "intact report"
    if args.purge:
        cleanup, run_dir = _purged_copy(run_dir, args.prefix)
        label = "after simulated 'maximal' purge (kite/tarean/dotplots removed)"

    try:
        broken, pages, assets = check_report(run_dir, args.prefix)
    finally:
        if cleanup:
            shutil.rmtree(cleanup, ignore_errors=True)

    print(f"report_linkcheck [{label}]: {pages} pages walked, "
          f"{assets} asset refs checked, {len(broken)} broken")
    for page, ref in broken[:25]:
        print(f"  BROKEN: {ref}   (referenced by {page})")
    if broken:
        print(f"FAIL: {len(broken)} broken local reference(s) in the report closure")
        return 1
    print("OK: report link-closure is fully self-contained")
    return 0


if __name__ == "__main__":
    sys.exit(main())
