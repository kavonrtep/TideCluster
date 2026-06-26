#!/usr/bin/env python3
"""Unit test for cross-TRC overlap resolution (tc_utils.resolve_trc_overlaps).

The clustering GFF3 may have different TRCs covering the same genomic bp
(variant arrays clustered separately, overlapping at boundaries). The resolver
must make it non-overlapping by giving each contested span to the dominant TRC
(largest total array length), without losing any base of the union and carrying
per-TRC attributes through. Pure-python, no external tools.

Run: python3 tests/test_overlap_resolution.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def feat(seqid, s, e, name, extra=""):
    a = f"Name={name};repeat_type=TR" + (";" + extra if extra else "")
    return f"{seqid}\tTideCluster\ttandem_repeat\t{s}\t{e}\t1\t+\t.\t{a}\n"


def parse(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = tc.Gff3Feature(line)
        rows.append((f.seqid, f.start, f.end, f.attributes_dict["Name"]))
    return rows


def union_and_overlap(rows):
    """(union_bp, multiply_covered_bp) over (seqid,start,end,name) rows."""
    by_seq = {}
    for seqid, s, e, _n in rows:
        by_seq.setdefault(seqid, []).append((s, e))
    covered = multi = 0
    for seqid, ivs in by_seq.items():
        ev = {}
        for s, e in ivs:
            ev[s] = ev.get(s, 0) + 1
            ev[e + 1] = ev.get(e + 1, 0) - 1
        cur = 0
        last = None
        for p in sorted(ev):
            if last is not None and cur:
                covered += p - last
                if cur >= 2:
                    multi += p - last
            cur += ev[p]
            last = p
    return covered, multi


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    d = tempfile.mkdtemp()
    gin = os.path.join(d, "clustering.gff3")
    gout = os.path.join(d, "resolved.gff3")
    # TRC_A dominant (total 99+200=299); TRC_B total 100; TRC_C total 40 nested in A
    with open(gin, "w") as f:
        f.write("##gff-version 3\n")
        f.write(feat("chr1", 1, 100, "TRC_A", "rDNA_type=45S"))
        f.write(feat("chr1", 1000, 1200, "TRC_A", "rDNA_type=45S"))
        f.write(feat("chr1", 50, 150, "TRC_B"))          # overlaps A on 50-100
        f.write(feat("chr1", 1010, 1050, "TRC_C"))        # nested in A -> absorbed

    before = parse(gin)
    tc.resolve_trc_overlaps(gin, gout)
    after = parse(gout)

    c0, m0 = union_and_overlap(before)
    c1, m1 = union_and_overlap(after)
    check("union (covered bp) preserved", c0 == c1)
    check("no residual cross-TRC overlap", m1 == 0)
    check("there was overlap to resolve", m0 > 0)

    by_name = {}
    for seqid, s, e, n in after:
        by_name.setdefault(n, []).append((s, e))
    # A keeps its contested span (50-100) -> A still ends at 100; B starts at 101
    check("dominant TRC_A wins contested span (A: 1-100)",
          (1, 100) in by_name.get("TRC_A", []))
    check("loser TRC_B trimmed to 101-150", by_name.get("TRC_B") == [(101, 150)])
    check("fully-nested TRC_C absorbed (dropped)", "TRC_C" not in by_name)
    check("A's second region kept whole after absorbing C",
          (1000, 1200) in by_name.get("TRC_A", []))

    # attributes carried through
    a_lines = [l for l in open(gout) if "Name=TRC_A" in l]
    check("per-TRC attributes preserved (rDNA_type on TRC_A)",
          all("rDNA_type=45S" in l for l in a_lines))
    check("repeat_type preserved", all("repeat_type=TR" in l for l in a_lines))

    # deterministic, sorted output
    starts = [s for _seq, s, _e, _n in after]
    check("output sorted by start", starts == sorted(starts))

    # idempotent: resolving an already-resolved file is a no-op on geometry
    gout2 = os.path.join(d, "resolved2.gff3")
    tc.resolve_trc_overlaps(gout, gout2)
    check("idempotent (second pass identical geometry)",
          parse(gout2) == after)

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
