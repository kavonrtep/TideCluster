#!/usr/bin/env python3
"""Unit tests for the HOR-order confidence tier (`hor_order_confidence`).

The tier distinguishes "founder recovered" from "confident HOR order"
(docs/hor_order_confidence_design.md). It is a pure function of columns the
monomer-size row already carries — it never changes founder selection. This
guards the bucketing and the high-k vacuity threshold against regressions, and
confirms the spectral classifier stays retired (no compute_hor / hor_status).

Run: python3 tests/test_hor_order_confidence.py   (exit 0 = pass)
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402

_fail = 0


def check(cond, msg):
    global _fail
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        _fail += 1


f = tc._hor_order_confidence

# (multiplicity, multiplicity_raw, irregular, founder_method, fallback) -> tier
CASES = [
    # mult == 1 -> no higher-order structure, regardless of method
    ((1, 1.0, False, "none", False),   "none"),
    ((1, 1.0, False, "ssr",  False),   "none"),
    # clean low-k strict divisor -> confident HOR order
    ((5, 5.0, False, "strict", False), "strict"),
    ((2, 2.01, False, "strict", False), "strict"),
    # family/anchor-rescued, clean, low-k -> supported
    ((8, 8.02, False, "pass2", False),     "supported"),
    ((6, 6.0,  False, "kh_deeper", False), "supported"),
    # irregular multiplicity -> weak even on the strict path (harmonic-basis)
    ((12, 12.4, True, "strict", False), "weak"),
    ((9,  9.3,  True, "pass2",  False), "weak"),
    # high k: integer test vacuous (k >= _HOR_ORDER_HIGH_K) -> weak
    ((280, 280.0, False, "strict", False), "weak"),
    # relaxed rescue paths -> weak
    ((3, 3.0, False, "pass3",   False), "weak"),
    ((4, 4.0, False, "cluster", False), "weak"),
    ((4, 4.0, False, "strict",  True),  "weak"),   # fallback overrides
]

for args, expected in CASES:
    got = f(*args)
    check(got == expected,
          f"tier{args} == {expected} (got {got})")

# The high-k cutoff is tied to the consensus-anchor tolerance constant.
check(tc._HOR_ORDER_HIGH_K == round(0.5 / tc._RATIO_TOL_FRAC),
      f"_HOR_ORDER_HIGH_K == round(0.5/_RATIO_TOL_FRAC) (= {tc._HOR_ORDER_HIGH_K})")
# Just below the cutoff a clean strict divisor still certifies the order.
below = tc._HOR_ORDER_HIGH_K - 1
check(f(below, float(below), False, "strict", False) == "strict",
      f"k={below} (just below cutoff) strict -> strict")
check(f(tc._HOR_ORDER_HIGH_K, float(tc._HOR_ORDER_HIGH_K), False, "strict", False) == "weak",
      f"k={tc._HOR_ORDER_HIGH_K} (at cutoff) strict -> weak")

# The retired spectral classifier must stay gone from the report module.
import tc_rerender_report as rr  # noqa: E402
check(not hasattr(rr, "compute_hor"),
      "tc_rerender_report has no compute_hor (spectral classifier retired)")
check(rr.__version__ == "3", f"report schema __version__ == 3 (got {rr.__version__})")
check(rr.HOR_REPORTED == ("strict", "supported"),
      "HOR (reported) == strict ∪ supported")

print()
if _fail:
    print(f"{_fail} FAILURE(S)")
    sys.exit(1)
print("ALL PASS")
