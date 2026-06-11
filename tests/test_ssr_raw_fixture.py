#!/usr/bin/env python3
"""Real-data regression for the SSR fix + founder pipeline, on a committed fixture.

`tests/data/kite_fixture/` holds a small slice of real kitehor output (the
*S. lycopersicum* TRC_2 = a 9.5 kb satellite that is ATC-rich, and TRC_18 = a
genuine ATC SSR family) + the matching clustering GFF3. Running
`build_monomer_size_csv` over it pins the behaviour that previously regressed:

  • TRC_2 chr2:15535707 — a clean satellite whose kitehor consensus_single SSR
    coverage is an inflated 96 % (raw 6 %) — keeps founder = 9520 (NOT collapsed
    to the 3 bp ATC), and carries both ssr_consensus_* (96 %) and ssr_raw_* (6 %).
  • TRC_18 — repeat_type=SSR — every array gets founder = motif length (3), and
    its raw per-array SSR coverage varies (68/82/18…%), reflecting real
    composition, not a flat consensus 100 %.

Self-contained (no external run dir / tools). Run: python3 tests/test_ssr_raw_fixture.py
"""
import collections
import csv
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402

FX = os.path.join(ROOT, "tests", "data", "kite_fixture")
_fail = 0


def check(cond, msg):
    global _fail
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        _fail += 1


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


out = os.path.join(tempfile.mkdtemp(), "monomer_size_top3_estimats.csv")
tc.build_monomer_size_csv(
    kite_tsv=os.path.join(FX, "kitehor.kite.tsv"),
    ssr_tsv=os.path.join(FX, "kitehor.ssr.tsv"),
    rescored_peaks_tsv=os.path.join(FX, "kitehor.rescored.peaks.tsv"),
    out_csv=out,
    tandem_validate_tsv=os.path.join(FX, "kitehor.tandem_validate.tsv"),
    trc_repeat_type=tc.parse_trc_ssr_motif_len(os.path.join(FX, "clustering.gff3")),
)
rows = list(csv.DictReader(open(out), delimiter="\t"))
by = {(r["TRC_ID"], r["seqid"], r["start"]): r for r in rows}

# Header carries the renamed consensus + new raw SSR columns.
hdr = set(rows[0].keys())
for c in ("ssr_consensus_total_coverage_pct", "ssr_raw_dominant_motif",
          "ssr_raw_total_coverage_pct", "ssr_raw_n_regions"):
    check(c in hdr, f"CSV has column {c}")
check("ssr_total_coverage_pct" not in hdr,
      "old ambiguous ssr_total_coverage_pct column is gone (renamed)")

# TRC_2 satellite: NOT collapsed to the SSR motif; both coverage views present.
p = by.get(("TRC_2", "chr2", "15535707"))
check(p is not None, "TRC_2 chr2:15535707 array present")
if p:
    check(p["founder_period"] == "9520",
          f"TRC_2 satellite keeps founder 9520, not 3 (got {p['founder_period']})")
    check(p["founder_method"] != "ssr" and p["ssr_founder_override"] == "false",
          "TRC_2 satellite not flagged as SSR founder")
    check(abs((_num(p["ssr_consensus_total_coverage_pct"]) or 0) - 96.23) < 0.1,
          "TRC_2 consensus coverage ~96 % (the artifact, kept + labelled)")
    check((_num(p["ssr_raw_total_coverage_pct"]) or 0) < 10,
          f"TRC_2 RAW coverage is low ~6 % (got {p['ssr_raw_total_coverage_pct']})")

# TRC_18 genuine SSR: founder = motif length; raw coverage varies per array.
t18 = [r for r in rows if r["TRC_ID"] == "TRC_18"]
check(t18 and all(r["founder_period"] == "3" for r in t18),
      f"TRC_18 (SSR) every array founder = 3 (n={len(t18)})")
raw18 = sorted({round(_num(r["ssr_raw_total_coverage_pct"]) or 0) for r in t18})
check(len(raw18) >= 3 and max(raw18) - min(raw18) > 20,
      f"TRC_18 raw per-array coverage varies (real composition): {raw18}")

print()
if _fail:
    print(f"{_fail} FAILURE(S)")
    sys.exit(1)
print("ALL PASS")
