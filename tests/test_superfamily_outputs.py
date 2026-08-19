#!/usr/bin/env python3
"""Regression tests for issue #7: superfamily outputs lost to a too-long filename.

`compare_trc_by_blast.R` named each dotplot after every member TRC id joined
together ("TRCS_10_17_74_..._931.png"). At ~61 three-digit members that passes
NAME_MAX (255 bytes), `png()` cannot open the device, and the script dies on the
FIRST loop iteration -- superfamilies are ordered by decreasing member count, so
the largest one is always i = 1. Everything written after the loop is lost: the
CSV, the manifest, `closePage()`. The v2 report then found no CSV and stated
"No TRC superfamilies were identified" for a run that had computed 75 of them,
and `run_all` still exited 0 because no `run_cmd` call site looked at its status.

Covered here:
  1. the R filename helper is bounded regardless of superfamily size (skipped
     if Rscript is unavailable)
  2. `load_superfamilies` finds the new index-based dotplot name, AND still
     finds the pre-fix "TRCS_*" name so old output directories re-render
  3. an ABSENT superfamily CSV is reported as a failed analysis, not as
     "no superfamilies found" -- the wrong-claim half of the bug
  4. `_run_required` turns a failed step into a non-zero exit; `_run_optional`
     warns and continues

Run: python3 tests/test_superfamily_outputs.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

failures = []


def check(label, got, want):
    ok = got == want
    print(F"  {'ok  ' if ok else 'FAIL'} {label}: got {got!r}, want {want!r}")
    if not ok:
        failures.append(label)


def check_true(label, cond, detail=""):
    print(F"  {'ok  ' if cond else 'FAIL'} {label}{'' if cond else ' -- ' + detail}")
    if not cond:
        failures.append(label)


# --- 1. the R helper --------------------------------------------------------
print("R: superfamily_dotplot_filename is bounded for any superfamily size")
rscript = shutil.which("Rscript")
if not rscript:
    print("  skip (Rscript not on PATH)")
else:
    # Source just the helper out of the script (the script itself needs the
    # full Bioconductor stack and CLI args).
    src = Path(ROOT, "tarean", "compare_trc_by_blast.R").read_text().splitlines()
    start = next(i for i, l in enumerate(src)
                 if l.startswith("superfamily_dotplot_filename <- function"))
    end = next(i for i in range(start + 1, len(src)) if src[i].startswith("}"))
    helper = "\n".join(src[start:end + 1])
    prog = helper + """
ns <- c(1, 61, 62, 75, 200, 999, 1000, 5000)
cat(max(sapply(ns, function(i) nchar(superfamily_dotplot_filename(i)))), "\\n")
cat(superfamily_dotplot_filename(1), "\\n")
cat(superfamily_dotplot_filename(75), "\\n")
"""
    out = subprocess.run([rscript, "-e", prog], capture_output=True, text=True)
    if out.returncode != 0:
        check_true("helper runs", False, out.stderr.strip())
    else:
        lines = out.stdout.split()
        check("max filename length over 1..5000 superfamilies stays << 255",
              int(lines[0]) <= 255, True)
        check("name is indexed by rank, not by member ids", lines[1],
              "superfamily_001.png")
        check("75th superfamily", lines[2], "superfamily_075.png")

    # The failing case from the issue: 64 three-digit members produced a
    # 261-byte name. Assert the old scheme really does exceed NAME_MAX, so this
    # test fails loudly if someone reverts to it.
    ids = [10, 17, 74, 142, 148, 150, 163, 172, 186, 188, 194, 199, 203, 212,
           270, 294, 295, 306, 307, 309, 311, 314, 316, 328, 336, 351, 359, 386,
           408, 413, 415, 485, 497, 579, 600, 601, 615, 616, 619, 620, 628, 633,
           634, 640, 648, 650, 654, 698, 718, 729, 737, 741, 746, 747, 780, 812,
           819, 822, 848, 870, 886, 921, 922, 931]
    legacy = "TRCS_" + "_".join(str(i) for i in ids) + ".png"
    check("the reported legacy name (64 members) is 261 bytes", len(legacy), 261)
    check_true("legacy name exceeds NAME_MAX", len(legacy) > 255)


# --- 2/3. the report side ---------------------------------------------------
import tc_rerender_report as rr

CSV = '"Superfamily","TRC","fallback"\n'


def make_run(tmp, csv_rows, dotplot_names=()):
    """Minimal run dir: superfamily CSV (optional) + dotplots/."""
    d = Path(tmp)
    (d / "dotplots").mkdir(parents=True, exist_ok=True)
    for n in dotplot_names:
        (d / "dotplots" / n).write_bytes(b"\x89PNG\r\n\x1a\n")
    if csv_rows is not None:
        (d / "p_trc_superfamilies.csv").write_text(CSV + csv_rows)
    return d


print("report: dotplot lookup (new name, legacy name, absent)")
tmp = tempfile.mkdtemp(prefix="tc_sf_test_")
try:
    rows = ('1,"TRC_10",False\n1,"TRC_17",False\n1,"TRC_74",False\n'
            '2,"TRC_5",False\n2,"TRC_9",True\n')

    d = make_run(os.path.join(tmp, "new"), rows,
                 ["superfamily_001.png", "superfamily_002.png"])
    sf = rr.load_superfamilies(rr.resolve_paths(d, "p"), d)
    check("two superfamilies parsed", [x["id"] for x in sf], [1, 2])
    check("members sorted numerically", sf[0]["trcs"],
          ["TRC_10", "TRC_17", "TRC_74"])
    check("new index-based dotplot found", sf[0]["dotplot"],
          "dotplots/superfamily_001.png")
    check("second superfamily's dotplot", sf[1]["dotplot"],
          "dotplots/superfamily_002.png")
    check("fallback column still read", sf[1]["fallback_trcs"], ["TRC_9"])

    # pre-fix output dir: only the old TRCS_* names exist
    d = make_run(os.path.join(tmp, "legacy"), rows,
                 ["TRCS_10_17_74.png", "TRCS_5_9.png"])
    sf = rr.load_superfamilies(rr.resolve_paths(d, "p"), d)
    check("legacy dotplot name still found (old run dirs re-render)",
          sf[0]["dotplot"], "dotplots/TRCS_10_17_74.png")
    check("legacy name for second superfamily", sf[1]["dotplot"],
          "dotplots/TRCS_5_9.png")

    # new name wins when both are present
    d = make_run(os.path.join(tmp, "both"), rows,
                 ["superfamily_001.png", "TRCS_10_17_74.png", "TRCS_5_9.png"])
    sf = rr.load_superfamilies(rr.resolve_paths(d, "p"), d)
    check("new name preferred over legacy", sf[0]["dotplot"],
          "dotplots/superfamily_001.png")

    d = make_run(os.path.join(tmp, "nopng"), rows)
    sf = rr.load_superfamilies(rr.resolve_paths(d, "p"), d)
    check("missing image -> dotplot None, superfamily still listed",
          (sf[0]["dotplot"], len(sf)), (None, 2))

    print("report: absent CSV is a failure, empty CSV is a result")
    d = make_run(os.path.join(tmp, "empty"), "")          # header only
    check("empty CSV -> csv_present True",
          bool(rr.resolve_paths(d, "p")["superfamilies_csv"]), True)
    d = make_run(os.path.join(tmp, "missing"), None)      # no CSV at all
    check("absent CSV -> csv_present False",
          bool(rr.resolve_paths(d, "p")["superfamilies_csv"]), False)

    # The rendered wording is the part that made the bug invisible.
    def render(csv_present):
        model = {"superfamilies": [], "superfamilies_csv_present": csv_present}
        out = Path(tmp) / F"sf_{csv_present}.html"
        ctx = {"up": "", "root_href": "index.html", "legacy_href": None,
               "assets_href": "assets"}
        run_meta = {"prefix": "p", "version": "test", "generated_at": "now",
                    "stats_line": ""}
        rr.render_superfamilies(model, out, run_meta, ctx)
        return out.read_text()

    html_ok = render(True)
    html_bad = render(False)
    check_true("empty CSV renders 'No TRC superfamilies were identified'",
               "No TRC superfamilies were identified" in html_ok, html_ok[:200])
    check_true("absent CSV does NOT claim none were identified",
               "No TRC superfamilies were identified" not in html_bad,
               html_bad[:200])
    check_true("absent CSV says the analysis did not complete",
               "did not complete" in html_bad, html_bad[:200])
finally:
    shutil.rmtree(tmp, ignore_errors=True)


# --- 4. failures are no longer discarded ------------------------------------
print("pipeline: a failed step aborts instead of exiting 0")
sys.argv = [sys.argv[0]]
import TideCluster as tcl

try:
    tcl._run_required("exit 3", "fake step")
    check_true("_run_required raises on failure", False, "no exception")
except RuntimeError as exc:
    check_true("_run_required raises on failure", True)
    check_true("the error names the step", "fake step" in str(exc), str(exc))

try:
    tcl._run_required("true", "fake step")
    check_true("_run_required is quiet on success", True)
except RuntimeError as exc:
    check_true("_run_required is quiet on success", False, str(exc))

check("_run_optional reports failure but returns False",
      tcl._run_optional("exit 1", "fake optional step"), False)
check("_run_optional returns True on success",
      tcl._run_optional("true", "fake optional step"), True)

check("run_cmd still reports status as [cmd, 'error']",
      tcl.tc.run_cmd("exit 1")[1], "error")
check("run_cmd status on success", tcl.tc.run_cmd("true")[1], "ok")

print()
if failures:
    print(F"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("test_superfamily_outputs.py: all checks passed")
