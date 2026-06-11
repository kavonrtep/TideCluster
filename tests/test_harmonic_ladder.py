#!/usr/bin/env python3
"""Regression for the harmonic-ladder founder (Pass 7), on the committed fixture.

`tests/data/kite_fixture/` includes TRC_47 (S. lycopersicum chr5:68075763) — a
divergent satellite whose basic 179 bp monomer has id_med 0.687 (below the 0.70
gate) but is corroborated by a clean harmonic ladder (rungs at 1×,2×,3×,6×,12×
of 179, with the 6× ≈ 1073 HOR unit conserved at 0.91). The ladder pass must
recover founder = 179 with ×6; and a satellite WITHOUT a ladder (TRC_2, ~9.5 kb)
must be left untouched.

Run: python3 tests/test_harmonic_ladder.py   (exit 0 = pass)
"""
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


out = os.path.join(tempfile.mkdtemp(), "m.tsv")
tc.build_monomer_size_csv(
    kite_tsv=os.path.join(FX, "kitehor.kite.tsv"),
    ssr_tsv=os.path.join(FX, "kitehor.ssr.tsv"),
    rescored_peaks_tsv=os.path.join(FX, "kitehor.rescored.peaks.tsv"),
    out_csv=out,
    tandem_validate_tsv=os.path.join(FX, "kitehor.tandem_validate.tsv"),
    trc_repeat_type=tc.parse_trc_ssr_motif_len(os.path.join(FX, "clustering.gff3")),
)
rows = list(csv.DictReader(open(out), delimiter="\t"))

# Positive: TRC_47 divergent-HOR satellite -> ladder rescue to 179 ×6.
t47 = next(r for r in rows if r["TRC_ID"] == "TRC_47")
check(t47["founder_period"] == "179", f"TRC_47 founder = 179 (got {t47['founder_period']})")
check(t47["strongest_period"] == "1073", f"TRC_47 strongest = 1073 (got {t47['strongest_period']})")
check(t47["multiplicity"] == "6", f"TRC_47 multiplicity = 6 (got {t47['multiplicity']})")
check(t47["founder_method"] == "ladder", f"TRC_47 founder_method = ladder (got {t47['founder_method']})")
check(t47["hor_order_confidence"] == "supported",
      f"TRC_47 tier = supported (got {t47['hor_order_confidence']})")
check(t47["irregular_multiplicity"] == "false", "TRC_47 ×6 is a clean integer (not irregular)")
check(abs(float(t47["founder_id_med"]) - 0.687) < 0.01,
      f"TRC_47 founder_id_med honestly ~0.69 (the diverged base; got {t47['founder_id_med']})")

# Negative: a satellite with NO harmonic ladder must NOT be touched by Pass 7.
for r in rows:
    if r["TRC_ID"] == "TRC_2":
        check(r["founder_method"] != "ladder",
              f"TRC_2 satellite not ladder-rescued ({r['seqid']}:{r['start']} method={r['founder_method']})")

# The helper itself: a clean ×6 ladder returns the fundamental; no ladder -> None.
peaks = [{"period": p, "identity_med": idm, "identity_iqr": iqr, "scan_occupancy_frac": occ}
         for p, idm, iqr, occ in [
             (179, 0.687, 0.045, 0.91), (358, 0.686, 0.027, 0.91),
             (537, 0.674, 0.017, 0.99), (1073, 0.913, 0.119, 1.0)]]
res = tc._harmonic_ladder_founder(peaks, 1073, 0.913)
check(res is not None and int(res[1]) == 179 and res[2] == 6,
      f"_harmonic_ladder_founder clean ladder -> (179, ×6) (got {res and (res[1], res[2])})")
# Only two rungs, double NOT more conserved -> rejected (no exceptional ×2).
flat = [{"period": 179, "identity_med": 0.687, "identity_iqr": 0.045, "scan_occupancy_frac": 0.91},
        {"period": 358, "identity_med": 0.686, "identity_iqr": 0.027, "scan_occupancy_frac": 0.91}]
check(tc._harmonic_ladder_founder(flat, 358, 0.686) is None,
      "×2 with no HOR-conservation gap is rejected (358 not more conserved than 179)")
# Exceptionally-clean ×2 (double genuinely more conserved) -> accepted.
clean2 = [{"period": 111, "identity_med": 0.65, "identity_iqr": 0.04, "scan_occupancy_frac": 0.9},
          {"period": 222, "identity_med": 0.86, "identity_iqr": 0.05, "scan_occupancy_frac": 1.0}]
r2 = tc._harmonic_ladder_founder(clean2, 222, 0.86)
check(r2 is not None and int(r2[1]) == 111 and r2[2] == 2,
      f"exceptionally-clean ×2 (conserved dimer) accepted (got {r2 and (r2[1], r2[2])})")

print()
if _fail:
    print(f"{_fail} FAILURE(S)")
    sys.exit(1)
print("ALL PASS")
