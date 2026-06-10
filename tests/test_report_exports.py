#!/usr/bin/env python3
"""Unit tests for the report TSV bundle + lazy (lightweight) Details child rows.

Guards two things added together (docs/report_tsv_and_legend_plan.md):
  • the report/data TSV bundle (trc_table / tra_table / tra_peaks / columns)
    mirrors the model and is fully documented by columns.tsv;
  • the per-array Details child row is NOT inlined as HTML — it is built lazily
    in JS from a compact embedded per-TRC JSON, so satellite-rich dashboards
    stay light. This test fails if the heavy inline-HTML approach ever returns.

Run: python3 tests/test_report_exports.py   (exit 0 = pass)
"""
import csv
import glob
import json
import os
import re
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_rerender_report as rr  # noqa: E402

_fail = 0


def check(cond, msg):
    global _fail
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        _fail += 1


def _peak(rank, period, score, idm, occ=0.5):
    # kitehor.rescored.peaks.tsv raw column names.
    return {
        "rank": str(rank), "period": str(period), "score": str(score),
        "identity_med": str(idm), "identity_iqr": "0.02",
        "scan_occupancy_frac": str(occ), "scan_n_intervals": "100",
        "coverage_frac": "0.9", "spatial_contrast": "0.1",
        "kmer_autocorr_founder": "0.8", "kmer_phase_contrast": "0.7",
        "subrepeat": "false", "phantom": "false", "founder_period": "53",
    }


def _array(seqid, start, end, founder, strongest, mult, tier, peaks):
    return {
        "seqid": seqid, "start": start, "end": end, "length": end - start,
        "founder_period": founder, "founder_id_med": 0.98,
        "delta_id_pp": 1.5, "strongest_period": strongest,
        "strongest_id_med": 0.95, "multiplicity": mult,
        "multiplicity_raw": float(mult), "hor_order_confidence": tier,
        "irregular_multiplicity": False, "founder_method": "strict",
        "founder_fallback": False, "cluster_rescue": False,
        "subrepeat_1_period": None, "subrepeat_1_tier": None,
        "subrepeat_1_occ": None, "subrepeat_2_period": None,
        "subrepeat_2_tier": None, "subrepeat_2_occ": None,
        "weak_short_founder_flag": False, "alt_longer_period": None,
        "alt_cluster_1_period": None, "alt_cluster_2_period": None,
        "ssr_dominant_motif": None, "ssr_total_coverage_pct": None,
        "copy_number": (end - start) / founder, "consensus_period_bp": founder,
        "rescored_peaks": peaks,
    }


def _model():
    a1 = _array("chr1", 100, 5400, 53, 159, 3, "strict",
                [_peak(1, 53, 0.87, 0.98), _peak(2, 159, 0.40, 0.95)])
    a2 = _array("chr1", 9000, 14300, 53, 53, 1, "none",
                [_peak(1, 53, 0.80, 0.97)])
    trc = {
        "id": "TRC_1", "repeat_type": "TR", "n_arrays": 2, "total_size": 10600,
        "min_array": 5300, "median_array": 5300, "max_array": 5300,
        "ssr_motif": None, "superfamily": "SF_1", "annotation": "",
        "kite": {"monomer_primary": 53, "prevalent_founder": 53, "n_hor": 1,
                 "n_hor_weak": 0, "n_subrepeat": 0, "n_ssr": 0},
        "tarean": {"monomer_length": 53, "total_score": 1.0},
        "arrays": [a1, a2],
    }
    return {"trcs": [trc]}


# ---- TSV bundle ----------------------------------------------------------
model = _model()
tmp = tempfile.mkdtemp()
import pathlib
counts = rr.write_table_exports(model, pathlib.Path(tmp))
check(counts == {"trc": 1, "tra": 2, "peaks": 3},
      f"write_table_exports counts == 1/2/3 (got {counts})")

trc_rows = list(csv.DictReader(open(f"{tmp}/trc_table.tsv"), delimiter="\t"))
tra_rows = list(csv.DictReader(open(f"{tmp}/tra_table.tsv"), delimiter="\t"))
pk_rows  = list(csv.DictReader(open(f"{tmp}/tra_peaks.tsv"), delimiter="\t"))
col_rows = list(csv.DictReader(open(f"{tmp}/columns.tsv"), delimiter="\t"))

check(len(trc_rows) == 1 and trc_rows[0]["n_hor_confident"] == "1",
      "trc_table: 1 row, n_hor_confident=1")
check([r["hor_order_confidence"] for r in tra_rows] == ["strict", "none"],
      "tra_table: 2 rows, tiers strict/none")
check(len(pk_rows) == 3, f"tra_peaks: 3 rows (got {len(pk_rows)})")
check(pk_rows[0]["TRC_ID"] == "TRC_1" and pk_rows[0]["period"] == "53"
      and "tier" in pk_rows[0], "tra_peaks: keyed + unfold columns present")

# columns.tsv documents every exported column of every table.
documented = {(r["table"], r["column"]) for r in col_rows}
for tbl, path in [("trc_table", "trc_table.tsv"), ("tra_table", "tra_table.tsv"),
                  ("tra_peaks", "tra_peaks.tsv")]:
    hdr = open(f"{tmp}/{path}").readline().rstrip("\n").split("\t")
    miss = [h for h in hdr
            if (tbl, h) not in documented and h not in ("TRC_ID", "seqid", "start", "end")]
    check(not miss, f"columns.tsv documents all {tbl} columns (missing: {miss})")

# ---- lazy child row: payload is DATA, not inline HTML (the perf guard) ----
script = rr._peaks_payload_script(model["trcs"][0]["arrays"])
check("<td" not in script and "<table" not in script and "tc-tier" not in script,
      "embedded peaks payload carries NO inline child HTML (data only)")
m = re.search(r'id="tc-peaks">(.*?)</script>', script, re.S)
payload = json.loads(m.group(1).replace("<\\/", "</"))
check(len(payload) == 2, f"payload has 2 arrays (got {len(payload)})")
any_arr = next(iter(payload.values()))
check(len(any_arr["pk"][0]) == len(rr._PEAK_COLUMNS),
      "each payload peak row is positional over _PEAK_COLUMNS")
cm = re.search(r'id="tc-coldict">(.*?)</script>', script, re.S)
cd = json.loads(cm.group(1).replace("<\\/", "</"))
check(len(cd["cols"]) == len(rr._PEAK_COLUMNS) and "tiers" in cd,
      "coldict carries column dictionary + tier map for the JS builder")

# ---- size guardrail on real output (skip if the fixture isn't present) ----
FIX = os.path.join(ROOT, "test_data/Solanum_lycopersicum/run_e2e/tc_report/trc")
if os.path.isdir(FIX):
    dashboards = glob.glob(f"{FIX}/*.html")
    inline = [f for f in dashboards if "data-details=" in open(f).read()]
    check(not inline, f"fixture: no dashboard uses inline data-details ({len(inline)} do)")
    sizes = {os.path.basename(f): os.path.getsize(f) for f in dashboards}
    biggest = max(sizes.values())
    # Pre-change TRC_4.html was ~5.3 MB; the lazy design must keep it well under.
    check(biggest < 3_000_000,
          f"fixture: largest dashboard < 3 MB (got {biggest/1e6:.1f} MB)")
else:
    print("SKIP fixture size guardrail (no run_e2e output present)")

# ---- R1: single-source guards (legend ≡ COLUMN_DICT ≡ TIER_DEFS) ----
check(set(rr.COLUMN_DICT) == {"trc_table", "tra_table", "tra_peaks"},
      "COLUMN_DICT has the three report tables")
for tbl, specs in rr.COLUMN_DICT.items():
    check(all(c.header and c.desc for c in specs),
          f"COLUMN_DICT[{tbl}]: every column has a header + description")
# columns.tsv documents exactly COLUMN_DICT (same set, same order per table).
documented_by_table = {}
for r in col_rows:
    documented_by_table.setdefault(r["table"], []).append(r["column"])
for tbl, specs in rr.COLUMN_DICT.items():
    check([c.header for c in specs] == documented_by_table.get(tbl),
          f"columns.tsv {tbl} rows == COLUMN_DICT[{tbl}] (order + set)")
# The HOR-order tier text is single-sourced: every TIER_DEFS meaning must
# appear verbatim in the rendered legend (fails if the legend re-hardcodes it).
legend = rr.arrays_legend()
for t, (label, meaning) in rr.TIER_DEFS.items():
    if t == "none":
        continue
    check(meaning in legend, f"legend embeds TIER_DEFS['{t}'] meaning (single source)")
# The peak-column legend derives from COLUMN_DICT["tra_peaks"].
check([h for h, _ in rr._DETAILS_COL_DESC] == [c.header for c in rr.COLUMN_DICT["tra_peaks"]],
      "_DETAILS_COL_DESC (legend) derives from COLUMN_DICT['tra_peaks']")

print()
if _fail:
    print(f"{_fail} FAILURE(S)")
    sys.exit(1)
print("ALL PASS")
