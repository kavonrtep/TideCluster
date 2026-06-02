#!/usr/bin/env python3
"""Unit tests for Pass 5 (peak-cluster fallback rescue) in tc_utils.

Drives _cluster_peaks_by_period() and _cluster_rescue_founder() with
synthetic peak dicts that mimic the kitehor.rescored.peaks.tsv schema.
No external data, no kitehor run required.

Cases mirror the analysis_kite15k findings (TRC_2 / TRC_4):
  (a) fallback row where cluster gate passes → rescue picks ~10 kb
  (b) fallback row where coverage/spatial gate fails → no rescue
  (c) singleton-cluster fallback → no rescue (n_peaks < 2)
  (d) score-margin too small → no rescue
  (e) clustering correctly groups adjacent periods via single link
  (f) chosen founder's cluster is excluded from alt-cluster list

Run: python3 tests/test_cluster_rescue.py   (exit 0 = pass)
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def peak(rank, period, score, id_med=0.50, cov_frac=0.40, spatial=1.0,
         subrepeat="false"):
    """Build one rescored-peaks row as a dict (string values like the
    real TSV reader produces via csv.DictReader)."""
    return {
        "rank":              str(rank),
        "period":            str(period),
        "score":             f"{score:.10f}",
        "identity_med":      f"{id_med:.4f}",
        "coverage_frac":     f"{cov_frac:.4f}",
        "spatial_contrast":  f"{spatial:.4f}",
        "subrepeat":         subrepeat,
        # diagnostic columns Pass 5 doesn't read but the row format keeps:
        "peak_height":       "1000",
        "phantom":           "false",
        "shift_med":         "0",
        "shift_consistency": "0.5",
        "kmer_autocorr_founder": "",
        "kmer_phase_contrast":   "",
        "scan_n_intervals":      "0",
        "scan_occupancy_frac":   "0",
    }


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    # ---- (a) TRC_2/TRC_4-style fallback: ~10 kb cluster dominates ----
    # Mirrors TRC_4:CP116281.2_2084208_2175313. Rank 1 (period 310) has
    # low coverage; ~10 kb peaks (9384, 9387, 10460-10914) collectively
    # outscore it with high coverage_frac / spatial_contrast.
    peaks_a = [
        peak(1,    310, 0.0263, id_med=0.49, cov_frac=0.07, spatial=0.14),
        peak(2,   9384, 0.0224, id_med=0.51, cov_frac=0.41, spatial=1.0),
        peak(3,     82, 0.0095, id_med=0.50, cov_frac=0.05, spatial=0.18),
        peak(4,    185, 0.0091, id_med=0.50, cov_frac=0.04, spatial=0.14),
        peak(5,   9387, 0.0071, id_med=0.51, cov_frac=0.41, spatial=1.0),
        peak(6,  10687, 0.0069, id_med=0.54, cov_frac=0.27, spatial=1.0),
        peak(7,  10465, 0.0068, id_med=0.53, cov_frac=0.30, spatial=0.95),
        peak(8,  10050, 0.0065, id_med=0.52, cov_frac=0.22, spatial=1.0),
        peak(9,  10914, 0.0062, id_med=0.51, cov_frac=0.24, spatial=1.0),
        peak(10, 10462, 0.0061, id_med=0.53, cov_frac=0.30, spatial=0.95),
        peak(11, 10460, 0.0055, id_med=0.53, cov_frac=0.30, spatial=0.95),
        peak(12, 10690, 0.0020, id_med=0.52, cov_frac=0.22, spatial=1.0),
        peak(13, 10375, 0.0015, id_med=0.53, cov_frac=0.28, spatial=0.95),
    ]
    res = tc._cluster_rescue_founder(peaks_a)
    check("(a) cluster rescue fires on TRC_4-style fallback row",
          res is not None)
    if res is not None:
        cluster, founder = res
        check("(a) founder is in the ~10 kb cluster, not 310",
              8000 <= founder <= 11500)
        check("(a) chosen cluster has >= 2 peaks (single-link)",
              cluster["n_peaks"] >= 2)

    # ---- (b) Gate fail: ~10 kb cluster has low cov_frac ----
    # Same period layout but every ~10 kb peak has cov_frac=0.05 (below
    # _CLUSTER_COV_FRAC_MIN). The rescue must abstain rather than pick a
    # noisy long-period winner.
    peaks_b = []
    for p in peaks_a:
        q = dict(p)
        if int(q["period"]) >= 5000:
            q["coverage_frac"] = "0.05"
            q["spatial_contrast"] = "0.20"
        peaks_b.append(q)
    res = tc._cluster_rescue_founder(peaks_b)
    check("(b) cluster rescue abstains when cov/spatial gate fails",
          res is None)

    # ---- (c) Singleton cluster: lone ~10 kb peak with no neighbour ----
    # Even with score >= rank1, n_peaks=1 should fail _CLUSTER_MIN_PEAKS.
    peaks_c = [
        peak(1, 310,   0.030, cov_frac=0.07, spatial=0.14),
        peak(2, 10000, 0.040, cov_frac=0.40, spatial=1.0),
        peak(3, 200,   0.005, cov_frac=0.05, spatial=0.10),
    ]
    res = tc._cluster_rescue_founder(peaks_c)
    check("(c) singleton cluster (n=1) cannot win even with score > rank1",
          res is None)

    # ---- (d) Score-margin too small: cluster sum < rank1 score ----
    # Two ~10 kb peaks each scoring tiny — their sum stays below 310's score.
    peaks_d = [
        peak(1, 310,   0.030, cov_frac=0.07, spatial=0.14),
        peak(2, 10000, 0.010, cov_frac=0.40, spatial=1.0),
        peak(3, 10050, 0.005, cov_frac=0.40, spatial=1.0),
    ]
    res = tc._cluster_rescue_founder(peaks_d)
    check("(d) cluster rescue abstains when score_sum < rank1 score",
          res is None)

    # ---- (e) Single-link clustering: chain of adjacent periods ----
    # 9384 → 9387 (Δ=3) joins; 9387 → 10000 (Δ=613) does NOT join under
    # the ±5 % / ±100 bp window (window for 9387 is max(469, 100)=469);
    # 10000 → 10050 (Δ=50) joins separately.
    peaks_e = [
        peak(1,  9384, 0.020, cov_frac=0.40, spatial=1.0),
        peak(2,  9387, 0.018, cov_frac=0.40, spatial=1.0),
        peak(3, 10000, 0.010, cov_frac=0.40, spatial=1.0),
        peak(4, 10050, 0.009, cov_frac=0.40, spatial=1.0),
    ]
    clusters = tc._cluster_peaks_by_period(peaks_e)
    # Sorted by score_sum desc; expect 2 clusters of size 2 each.
    check("(e) two clusters formed from two pairs", len(clusters) == 2)
    check("(e) each cluster has 2 members",
          all(c["n_peaks"] == 2 for c in clusters))

    # ---- (f) Alt-cluster exclusion of the founder's own cluster ----
    # Reuse case (a). Recompute clusters_all and filter out the founder.
    res = tc._cluster_rescue_founder(peaks_a)
    if res is not None:
        _, founder = res
        clusters_all = tc._cluster_peaks_by_period(peaks_a)
        window = max(founder * tc._CLUSTER_WINDOW_PCT,
                     tc._CLUSTER_WINDOW_BP_MIN)
        def _is_founder(c):
            for m in c["members"]:
                mp = float(m["period"])
                if abs(mp - founder) <= window:
                    return True
            return False
        alts = [c for c in clusters_all if not _is_founder(c)]
        # 310 should appear as the top alt cluster (rejected founder).
        alt_periods = [int(round(c["median_period"])) for c in alts]
        check("(f) 310 surfaces as an alt cluster after rescue",
              any(abs(p - 310) <= 5 for p in alt_periods))

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
