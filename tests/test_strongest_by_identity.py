#!/usr/bin/env python3
"""Unit tests for Pass-1 `strongest = argmax(identity_med)` semantics.

Regression guard for the S. pimpinellifolium TRC_8 case where kitehor's
rescore picks the rank-1-by-score peak (a 178 bp basic monomer with
id_med 0.736) over a higher-identity long-period peak (~1786 bp with
id_med 0.818). The latter is the real 10x HOR period; trusting kitehor's
column blindly suppresses the HOR call (multiplicity=1 instead of 10).

Drives build_monomer_size_csv() end-to-end through a synthetic
rescored-peaks TSV; no kitehor / blast / external tools required.

Run: python3 tests/test_strongest_by_identity.py   (exit 0 = pass)
"""
import csv
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


# Schema mirrors kitehor.rescored.peaks.tsv. Only the columns Pass 1
# actually reads need real values; the rest can stay empty strings.
_COLS = [
    "case_id", "array_length", "rank", "period", "peak_height",
    "score", "score2", "score2_norm", "background",
    "identity_med", "identity_iqr", "identity_p25", "identity_n",
    "shift_med", "shift_consistency", "phantom", "subrepeat",
    "coverage_frac", "spatial_contrast", "founder_period",
    "kmer_autocorr_founder", "kmer_phase_contrast",
    "scan_n_intervals", "scan_occupancy_frac",
]


def _row(case_id, rank, period, score, id_med, kitehor_founder,
         cov_frac=1.0, spatial=0.0, id_iqr=0.0):
    """Build one rescored-peaks row as a dict (string values)."""
    out = {c: "" for c in _COLS}
    out["case_id"]          = case_id
    out["array_length"]     = "10000"
    out["rank"]             = str(rank)
    out["period"]           = str(period)
    out["score"]            = f"{score:.10f}"
    out["identity_med"]     = f"{id_med:.4f}"
    out["identity_iqr"]     = f"{id_iqr:.4f}"
    out["coverage_frac"]    = f"{cov_frac:.4f}"
    out["spatial_contrast"] = f"{spatial:.4f}"
    out["founder_period"]   = ("" if kitehor_founder is None
                                else str(kitehor_founder))
    out["phantom"]          = "false"
    out["subrepeat"]        = "false"
    out["peak_height"]      = "1000"
    return out


def _run(rows):
    """Write a synthetic rescored-peaks TSV, run build_monomer_size_csv,
    return the per-row output dict keyed by record_id."""
    d = tempfile.mkdtemp()
    rescored = os.path.join(d, "rescored.peaks.tsv")
    out_csv  = os.path.join(d, "out.csv")
    with open(rescored, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_COLS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    tc.build_monomer_size_csv(kite_tsv=None, ssr_tsv=None,
                              rescored_peaks_tsv=rescored,
                              out_csv=out_csv)
    with open(out_csv) as fh:
        return {f"{r['TRC_ID']}:{r['seqid']}_{r['start']}_{r['end']}": r
                for r in csv.DictReader(fh, delimiter="\t")}


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    # (a) S. pimpinellifolium TRC_8 chr12:1314680 mirror. Kitehor picks
    # the rank-1-by-score peak 178 (id_med 0.736) as founder_period;
    # the real high-identity peak is 1786 (id_med 0.818). Strongest
    # should override to 1786 → Pass-1 divisor search finds 178 as
    # the smallest passing divisor (1786/178 ≈ 10.03), giving founder
    # 178 with multiplicity 10.
    case = "TRC_8:chr12_1314680_1327267"
    rows_a = [
        _row(case, 1,  178, 0.0382, 0.736,  kitehor_founder=178),
        _row(case, 2,  356, 0.0146, 0.730,  kitehor_founder=178),
        _row(case, 3, 1962, 0.0143, 0.805,  kitehor_founder=178),
        _row(case, 4, 1788, 0.0117, 0.818,  kitehor_founder=178, id_iqr=0.11),
        _row(case, 5, 1786, 0.0044, 0.820,  kitehor_founder=178, id_iqr=0.11),
    ]
    out = _run(rows_a)
    r = out["TRC_8:chr12_1314680_1327267"]
    check("(a) strongest overridden to ~1786",
          1700 <= int(r["strongest_period"]) <= 1800)
    check("(a) founder is 178",
          int(r["founder_period"]) == 178)
    check("(a) multiplicity is 10",
          int(r["multiplicity"]) == 10)
    check("(a) kitehor_founder_period preserved as 178",
          int(r["kitehor_founder_period"]) == 178)

    # (b) Well-behaved array: kitehor's pick IS the argmax(id_med).
    # Strongest must equal kitehor_founder_period (no override fires)
    # and the existing 6× HOR call must be preserved.
    case2 = "TRC_8:chr5_581360_593243"
    rows_b = [
        _row(case2, 1, 1076, 0.0722, 0.885, kitehor_founder=1076),
        _row(case2, 2,  179, 0.0476, 0.721, kitehor_founder=1076),
        _row(case2, 3,  718, 0.0315, 0.700, kitehor_founder=1076),
        _row(case2, 4,  359, 0.0196, 0.700, kitehor_founder=1076),
    ]
    out = _run(rows_b)
    r = out["TRC_8:chr5_581360_593243"]
    check("(b) strongest unchanged when kitehor already argmax(id_med)",
          int(r["strongest_period"]) == 1076)
    check("(b) kitehor_founder_period also 1076",
          int(r["kitehor_founder_period"]) == 1076)
    check("(b) founder 179, multiplicity 6 (unchanged)",
          int(r["founder_period"]) == 179 and int(r["multiplicity"]) == 6)

    # (c) Kitehor NA → existing fallback path is preserved.
    # All peaks have id_med < 0.7 too, so the in-branch defensive
    # fallback is irrelevant; kitehor NA is what triggers Pass-5
    # eligibility downstream.
    case3 = "TRC_99:chrX_1_10000"
    rows_c = [
        _row(case3, 1, 310, 0.0300, 0.50, kitehor_founder=None),
        _row(case3, 2, 200, 0.0100, 0.50, kitehor_founder=None),
    ]
    out = _run(rows_c)
    r = out["TRC_99:chrX_1_10000"]
    check("(c) kitehor NA → founder_fallback=true",
          r["founder_fallback"] == "true")
    check("(c) kitehor NA → strongest = rank-1 period (310)",
          int(r["strongest_period"]) == 310)
    check("(c) kitehor NA → kitehor_founder_period blank",
          r["kitehor_founder_period"] == "")

    # (d) Defensive: kitehor reports a founder, but NO peak passes the
    # 0.7 id gate. We fall back to kitehor's pick (no behaviour change
    # vs 1.12.0 for this edge case).
    case4 = "TRC_50:chr2_500_1500"
    rows_d = [
        _row(case4, 1, 250, 0.0200, 0.60, kitehor_founder=250),
        _row(case4, 2, 125, 0.0100, 0.55, kitehor_founder=250),
    ]
    out = _run(rows_d)
    r = out["TRC_50:chr2_500_1500"]
    check("(d) all-below-gate → strongest falls back to kitehor pick (250)",
          int(r["strongest_period"]) == 250)
    check("(d) kitehor_founder_period reflects kitehor's pick (250)",
          int(r["kitehor_founder_period"]) == 250)

    # (e) Founder tiebreaker — drapa TRC_2 chr3:56557 case. Strongest is
    # 2176 (id_med 0.996). At kr=2 there are three valid divisors:
    # 1066 (k=2.041, id_med 0.972) — fuzzy, lower id
    # 1088 (k=2.000, id_med 0.995) — cleanest, highest id
    # 1091 (k=1.995, id_med 0.995) — clean, highest id (tie on id_med)
    # The old "smallest period wins" sort picked 1066 because numerically
    # smallest. The fix groups by kr and picks the best per group
    # (highest id_med → cleanest |k-kr| → smallest P), giving 1088.
    case5 = "TRC_2:chr3_56557_126534"
    rows_e = [
        _row(case5, 1, 1088, 0.5914, 0.9954, kitehor_founder=1088),
        _row(case5, 2, 1017, 0.0313, 0.4730, kitehor_founder=1088),
        _row(case5, 3, 1091, 0.0219, 0.9954, kitehor_founder=1088),
        _row(case5, 4, 1066, 0.0097, 0.9719, kitehor_founder=1088),
        _row(case5, 5, 2176, 0.0078, 0.9963, kitehor_founder=1088),
    ]
    out = _run(rows_e)
    r = out["TRC_2:chr3_56557_126534"]
    check("(e) tiebreaker picks 1088 over fuzzy 1066 at same kr=2",
          int(r["founder_period"]) == 1088)
    check("(e) strongest is 2176 (highest id_med)",
          int(r["strongest_period"]) == 2176)
    check("(e) multiplicity is 2",
          int(r["multiplicity"]) == 2)

    # (f) Across-kr semantic preserved — a kr=10 divisor must still win
    # over a kr=2 divisor when both qualify (basic-monomer rule). Set up
    # strongest=1786 with the highest id_med so it's picked as strongest;
    # candidates 178 (kr=10, k=10.034) and 892 (kr=2, k=2.0022) both pass
    # the gate. The new per-kr-winner sort across kr groups picks the
    # smallest period → 178 (deeper decomposition wins, as before).
    case6 = "TRC_99:chrX_1_50000"
    rows_f = [
        _row(case6, 1, 1786, 0.0500, 0.9900, kitehor_founder=1786),  # strongest
        _row(case6, 2,  178, 0.0400, 0.7500, kitehor_founder=1786),  # kr=10
        _row(case6, 3,  892, 0.0300, 0.7500, kitehor_founder=1786),  # kr=2
    ]
    out = _run(rows_f)
    r = out["TRC_99:chrX_1_50000"]
    check("(f) across-kr: smallest P wins (178 over 892) — basic monomer",
          int(r["founder_period"]) == 178)
    check("(f) strongest stays at 1786 (highest id_med)",
          int(r["strongest_period"]) == 1786)
    check("(f) multiplicity is 10",
          int(r["multiplicity"]) == 10)

    # (g) Lever 2 — drapa TRC_26 chr9:1816989 mirror. Strongest 1880,
    # nearby basic-monomer peaks 154/158/161/165 each fail strict
    # |k-kr|≤0.05 individually (gaps 0.21, 0.10, 0.32, 0.39), but
    # cluster as one group (within ±5%) with score-weighted mean ≈
    # 156.8 → k=11.99 → integer-clean. Lever 2's cluster-mean Pass-1
    # picks founder ≈ 157 with multiplicity = 12. Without Lever 2,
    # Pass 1 falls to founder=315 (=2×154 clean divisor) with mult=6,
    # collapsing the HOR call.
    case7 = "TRC_26:chr9_1816989_1825401"
    rows_g = [
        _row(case7, 1,  154, 0.3201, 0.9481, kitehor_founder=154),
        _row(case7, 2,  158, 0.1712, 0.9494, kitehor_founder=154),
        _row(case7, 3,  161, 0.1476, 0.9503, kitehor_founder=154),
        _row(case7, 4,  315, 0.0488, 0.9651, kitehor_founder=154),
        _row(case7, 5, 1880, 0.0170, 0.9915, kitehor_founder=154),  # strongest
        _row(case7, 6,  165, 0.0115, 0.9455, kitehor_founder=154),
        _row(case7, 7,  312, 0.0096, 0.9647, kitehor_founder=154),
    ]
    out = _run(rows_g)
    r = out["TRC_26:chr9_1816989_1825401"]
    check("(g) cluster-mean founder lands in basic-monomer band (150..170)",
          150 <= int(r["founder_period"]) <= 170)
    check("(g) multiplicity = 12 (was 6 with single-peak Pass 1)",
          int(r["multiplicity"]) == 12)
    check("(g) strongest stays at 1880 (highest id_med)",
          int(r["strongest_period"]) == 1880)

    # (h) Lever 3 — TRC with 3 arrays at the basic monomer (131) plus
    # one array where only a single peak in the basic band exists
    # (262 = 2×131 is the only clean divisor of strongest 3395). The
    # Lever 3 expansion lets Pass 2 re-evaluate that mult>1 row using
    # the TRC consensus founder (131) and rescue founder back to 131.
    # The 3 consensus-building arrays have a clean k=2.0 relationship
    # so they keep founder=131 / mult=2 untouched.
    trc_h = "TRC_500"
    rows_h = [
        # Three arrays at the basic monomer (mult=2 each → consensus 131).
        # The 262 peak's id_med is just-higher than 131's so it wins the
        # argmax and becomes strongest; Pass 1's divisor search then
        # finds 131 (k=2.00, clean) as founder.
        _row(f"{trc_h}:chr1_1_10000",     1, 131, 0.50, 0.95, kitehor_founder=131),
        _row(f"{trc_h}:chr1_1_10000",     2, 262, 0.10, 0.96, kitehor_founder=131),  # strongest
        _row(f"{trc_h}:chr1_20000_30000", 1, 131, 0.50, 0.95, kitehor_founder=131),
        _row(f"{trc_h}:chr1_20000_30000", 2, 262, 0.10, 0.96, kitehor_founder=131),
        _row(f"{trc_h}:chr2_1_10000",     1, 131, 0.50, 0.95, kitehor_founder=131),
        _row(f"{trc_h}:chr2_1_10000",     2, 262, 0.10, 0.96, kitehor_founder=131),
        # The "stuck" array — strongest 3395, only one peak in the basic
        # band (131, k=25.92, fails strict gate), peak at 262 happens to
        # be clean (k=12.96; not integer either, but 3395/262=12.96 is
        # clean enough — let's make it clean: strongest 3406, k=26.0/13.0).
    ]
    # Build the stuck array carefully: strongest must be a clean multiple
    # of 262 so Pass 1 picks 262 as founder, but it should ALSO be off
    # 131 by enough that 131 fails the strict ratio gate.
    # 262 × 13 = 3406. 3406 / 131 = 26.00 — actually that IS integer.
    # Need a strongest where 262 is clean divisor but 131 is NOT.
    # 262 × 13 = 3406 → 3406/131 = 26.0 (clean) — bad for the test.
    # 262 × 9 = 2358 → 2358/131 = 18.00 (clean) — bad.
    # Use a non-integer relation: make peak labels 131 and 262 but pick
    # strongest = 3395, peak 262 at k=3395/262=12.96 (NOT clean either).
    # Need 262 to be a clean divisor of strongest. Use strongest=3406,
    # which is 13×262 = 3406 exactly. But 3406/131=26 is also clean.
    # The biology is: the basic monomer is 131 and strongest is N×131,
    # so 262=2×131 is also a divisor. The "stuck" case is when 131 has
    # |k-kr|>0.05 but 262 has |k-kr|≤0.05 — that requires the peak at
    # 131 to be slightly off (e.g. 131 but strongest is at 3404).
    # 3404/131 = 25.985 → |k-kr|=0.015 (clean). Hmm.
    # 3404/262 = 12.992 → |k-kr|=0.008 (clean). Both clean.
    # To make 131 NOT clean while 262 IS clean, the strongest must be
    # such that strongest/131 has gap > 0.05. E.g. strongest=2620 / 131
    # = 20.00 (clean) — bad.
    # Try strongest=2624 / 131 = 20.031 (clean, gap 0.031).
    # strongest=2640 / 131 = 20.153 → gap 0.153 (FAILS strict).
    # strongest=2640 / 262 = 10.076 → gap 0.076 (also FAILS strict).
    # strongest=2620 / 262 = 10.00 clean; 2620/131 = 20 clean — bad.
    # strongest=2624 / 262 = 10.015 clean; 2624/131 = 20.031 clean.
    # The problem: if 262 = 2×131 EXACTLY, the ratios scale together.
    # Solution: make the peak labels 131 and 263 (slightly off 2×131).
    # peak 131, peak 263. strongest=3408 → 3408/131=26.015 clean;
    # 3408/263=12.96 → gap 0.04 clean. Both still pass.
    # Try strongest with 131 OFF but 263 clean: strongest 3419 → 131:
    # 3419/131=26.099 gap 0.099 FAILS; 3419/263=13.00 clean ✓.
    # Use that.
    rows_h.append(_row(f"{trc_h}:chr3_1_10000", 1,  131, 0.20, 0.90, kitehor_founder=131))
    rows_h.append(_row(f"{trc_h}:chr3_1_10000", 2,  263, 0.10, 0.90, kitehor_founder=131))
    rows_h.append(_row(f"{trc_h}:chr3_1_10000", 3, 3419, 0.05, 0.99, kitehor_founder=131))  # strongest
    out = _run(rows_h)
    # The three consensus-building arrays should report founder=131, mult=2
    for k in (f"{trc_h}:chr1_1_10000", f"{trc_h}:chr1_20000_30000",
              f"{trc_h}:chr2_1_10000"):
        ri = out[k]
        check(f"(h) consensus-array {k.split(':')[1]} founder=131 mult=2",
              int(ri["founder_period"]) == 131 and int(ri["multiplicity"]) == 2)
    # The "stuck" array: without Lever 3, Pass 1 would pick founder=263
    # (closest clean divisor of 3419 with kr=13). With Lever 3,
    # Pass 2 expansion swaps to founder=131 (TRC consensus) with
    # mult=26 (irregular).
    rs = out[f"{trc_h}:chr3_1_10000"]
    check("(h) Lever 3 rescues mult>1 stuck row to consensus 131",
          int(rs["founder_period"]) == 131)
    check("(h) Lever 3 sets irregular_multiplicity=true",
          rs["irregular_multiplicity"] == "true")

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
