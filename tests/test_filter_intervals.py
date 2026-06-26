#!/usr/bin/env python3
"""Equivalence guard: the sweep-line tc.filter_intervals() must return exactly
the same objects, in the same order, as the original O(N*M) nested loop.

tc.filter_intervals keeps a GRange g from `gr` iff some interval f in `gr_filter`
on the same seqid has f.name == g.name and g.within(f). The original
implementation was a quadratic nested loop; it was replaced by a per-(seqid,name)
binary-search sweep for speed. This test pins output-identity: it embeds the old
reference implementation and compares (by object identity and order) against the
shipped one on many randomized inputs, plus a few hand-built edge cases.

Self-contained, no external tools/data. Run: python3 tests/test_filter_intervals.py
"""
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def reference_filter_intervals(gr_filter, gr):
    """The original O(N*M) implementation, verbatim, as the ground truth."""
    dict_filter = {}
    dict_gr = {}
    for grange in gr_filter:
        dict_filter.setdefault(grange.seqid, []).append(grange)
    for grange in gr:
        dict_gr.setdefault(grange.seqid, []).append(grange)

    for seqid in dict_filter:
        dict_filter[seqid].sort(key=lambda x: x.start)
    for seqid in dict_gr:
        dict_gr[seqid].sort(key=lambda x: x.start)

    filtered = []
    for seqid in dict_filter:
        if seqid not in dict_gr:
            continue
        for grange_out in dict_gr[seqid]:
            for grange_filter in dict_filter[seqid]:
                if grange_out.name == grange_filter.name:
                    if grange_out.within(grange_filter):
                        filtered.append(grange_out)
                        break
    return filtered


def random_granges(rng, n, seqids, names, coord_max):
    out = []
    for _ in range(n):
        start = rng.randint(1, coord_max)
        end = start + rng.randint(0, coord_max // 4)
        out.append(tc.GRange(rng.choice(seqids), start, end,
                             rng.choice(names), rng.choice("+-")))
    return out


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    # --- randomized equivalence: identical object list (identity + order) ----
    rng = random.Random(12345)
    seqids = ["chr1", "chr2", "chr3"]
    names = ["TRC_1", "TRC_2", "TRC_3"]
    all_ok = True
    nonempty_seen = False
    for trial in range(400):
        # small coord_max => dense overlaps, the quadratic regime this fixes
        coord_max = rng.choice([20, 60, 200])
        gr_filter = random_granges(rng, rng.randint(0, 40), seqids, names, coord_max)
        gr = random_granges(rng, rng.randint(0, 40), seqids, names, coord_max)
        expected = reference_filter_intervals(gr_filter, gr)
        got = tc.filter_intervals(gr_filter, gr)
        if expected:
            nonempty_seen = True
        # identity-and-order equality: same objects, same sequence
        if len(got) != len(expected) or any(a is not b for a, b in zip(got, expected)):
            all_ok = False
            print(f"  trial {trial}: MISMATCH "
                  f"(expected {len(expected)} got {len(got)})")
            break
    check("randomized equivalence (object identity + order)", all_ok)
    check("randomized trials exercised non-empty results", nonempty_seen)

    # --- explicit edge cases ------------------------------------------------
    # exact-boundary containment (within is inclusive on both ends)
    f = [tc.GRange("c", 10, 20, "A", "+")]
    g = [tc.GRange("c", 10, 20, "A", "+")]
    check("inclusive boundary containment kept",
          [x.start for x in tc.filter_intervals(f, g)] == [10])

    # name mismatch => not kept even if geometrically contained
    f = [tc.GRange("c", 1, 100, "A", "+")]
    g = [tc.GRange("c", 10, 20, "B", "+")]
    check("name mismatch excluded", tc.filter_intervals(f, g) == [])

    # containment needs ONE interval to cover both ends (not a union of two)
    f = [tc.GRange("c", 1, 10, "A", "+"), tc.GRange("c", 9, 30, "A", "+")]
    g = [tc.GRange("c", 5, 25, "A", "+")]   # spans the gap; covered by 2nd only
    ref = reference_filter_intervals(f, g)
    check("no spurious union-containment",
          [x is y for x, y in zip(tc.filter_intervals(f, g), ref)] ==
          [True] * len(ref) and len(tc.filter_intervals(f, g)) == len(ref))

    # seqid present in gr but absent in filter => dropped
    f = [tc.GRange("c", 1, 100, "A", "+")]
    g = [tc.GRange("other", 10, 20, "A", "+")]
    check("seqid absent from filter dropped", tc.filter_intervals(f, g) == [])

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
