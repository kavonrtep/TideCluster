#!/usr/bin/env python3
"""Unit tests for Pass 7b — dominant-score harmonic-ladder founder.

Guards the FabTR-53 divergent-satellite fix: when argmax(id_med) `strongest` is a
long OFF-ladder period (e.g. 3150 = 16.76*188), the founder must be recovered as
the rank-1-by-score monomer P0 (188) via its self-standing harmonic ladder — and
the guards must leave commensurate / unladdered / kitehor-disagreeing arrays
alone. Self-contained: synthetic peak dicts, no external data or kitehor run.

Run: python3 tests/test_dominant_ladder.py   (exit 0 = pass)
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def peak(period, score=0.05, idm=0.90, iqr=0.04, occ=0.9):
    return {"period": period, "score": score, "identity_med": idm,
            "identity_iqr": iqr, "scan_occupancy_frac": occ}


def fabtr53_peaks():
    """Model ptg000030l_rc:2-44995 — P0=188 dominant by score with a 188/376/564
    ladder; strongest (argmax id_med) = 3150, off-ladder (3150/188=16.76)."""
    return [
        peak(188, score=0.179, idm=0.904, iqr=0.048, occ=0.95),  # rank-1 monomer
        peak(376, score=0.025, idm=0.899),                       # 2x rung
        peak(3150, score=0.015, idm=0.9835),                     # off-ladder strongest
        peak(563, score=0.005, idm=0.866),                       # 3x rung
    ]


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    # --- _commensurate ---
    check("3150 not commensurate with 188", not tc._commensurate(3150, 188))
    check("1786 commensurate with 178 (10x)", tc._commensurate(1786, 178))
    check("365 not commensurate with 189 (1.93x)", not tc._commensurate(365, 189))
    check("equal periods commensurate (k=1)", tc._commensurate(188, 188))
    check("376 commensurate with 188 (dimer)", tc._commensurate(376, 188))

    # --- _dominant_ladder_founder: recover P0=188 from off-ladder strongest ---
    pk = fabtr53_peaks()
    res = tc._dominant_ladder_founder(pk, strongest_period=3150, dominant_peak=pk[0])
    check("FabTR-53: founder recovered as 188", res is not None and res[1] == 188)
    check("FabTR-53: multiplicity = round(3150/188) = 17",
          res is not None and res[2] == 17)

    # No ladder (only P0 present, no rungs) -> abstain.
    lone = [peak(188, score=0.179), peak(3150, score=0.015, idm=0.98)]
    check("no ladder rungs -> None",
          tc._dominant_ladder_founder(lone, 3150, lone[0]) is None)

    # P0 fails the IQR gate -> abstain.
    noisy = [peak(188, score=0.179, iqr=0.30), peak(376), peak(563), peak(3150, idm=0.98)]
    check("P0 high IQR -> None",
          tc._dominant_ladder_founder(noisy, 3150, noisy[0]) is None)

    # dominant_peak above strongest -> abstain (P0 must be the fundamental).
    check("P0 >= strongest -> None",
          tc._dominant_ladder_founder(pk, 100, pk[0]) is None)

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
