#!/usr/bin/env python3
"""Unit tests for get_repeatmasker_annotation() coverage accounting.

Regression guard for the first-hit-only bug: the per-annotation coverage must be
the UNION of all RepeatMasker hit intervals for a query (bounded to [0, 1]), not
just the first hit. Self-contained: builds synthetic RepeatMasker .out files, no
external data or RepeatMasker run required.

Run: python3 tests/test_annotation_coverage.py   (exit 0 = pass)
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def rm_line(seqid, start, end, cls, name="ref1"):
    # RepeatMasker .out column layout: idx4=query, 5=qstart, 6=qend,
    # 8=strand, 9=matching repeat, 10=class/family. >=12 cols required.
    return f"300 5.0 0.0 0.0 {seqid} {start} {end} (0) + {name} {cls} (0) 1 {end} *\n"


def annotate(rm_lines, seq_lengths, parse_id=False):
    d = tempfile.mkdtemp()
    rm_path = os.path.join(d, "x.out")
    with open(rm_path, "w") as fh:
        fh.write("h1\nh2\nh3\n")  # 3 header lines are skipped by the parser
        fh.writelines(rm_lines)
    tc.get_repeatmasker_annotation(rm_path, seq_lengths, os.path.join(d, "tc"),
                                   parse_id)
    out = {}
    with open(os.path.join(d, "tc_annotation.tsv")) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and c[0]:
                out[c[0]] = (c[1], float(c[2]))
    return out


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    # Three tandem, non-overlapping hits of one family tile a 600 bp consensus.
    # The first-hit bug reported 200/600 = 0.333; the union must be ~1.0.
    r = annotate([rm_line("s0", 1, 200, "CEN"),
                  rm_line("s0", 201, 400, "CEN"),
                  rm_line("s0", 401, 600, "CEN")], {"s0": 600})
    check("three tandem hits -> 1.0", abs(r["s0"][1] - 1.0) < 1e-6)

    # Overlapping hits of one family are union-ed, not summed.
    # 1-200 and 150-400 -> union 1-400 = 400/600; sum would be 451/600.
    r = annotate([rm_line("s0", 1, 200, "CEN"),
                  rm_line("s0", 150, 400, "CEN")], {"s0": 600})
    check("overlapping hits union -> 0.667", abs(r["s0"][1] - 400 / 600) < 1e-6)

    # Two families on one sequence are tracked independently.
    r = annotate([rm_line("s0", 1, 300, "CEN"),
                  rm_line("s0", 301, 600, "rDNA")], {"s0": 600})
    check("two families each 0.5",
          abs(r["s0"][1] - 0.5) < 1e-6 and len(r) == 1)

    # parse_id=True collapses per-array consensus IDs to TRC_<n> and averages.
    # Two consensus of TRC_1, each fully covered -> mean 1.0.
    r = annotate([rm_line("TRC_1_rep0_x", 1, 600, "CEN"),
                  rm_line("TRC_1_rep1_x", 1, 600, "CEN")],
                 {"TRC_1_rep0_x": 600, "TRC_1_rep1_x": 600}, parse_id=True)
    check("parse_id mean over TRC consensus -> 1.0",
          abs(r["TRC_1"][1] - 1.0) < 1e-6)

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
