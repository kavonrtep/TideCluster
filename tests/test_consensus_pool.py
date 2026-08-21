#!/usr/bin/env python3
"""Unit tests for tc_utils.write_consensus_sequences_all.

`<prefix>_consensus/consensus_sequences_all.fasta` is a required input of
tc_comparative_analysis.R. It used to be assembled inside `annotation()`,
because its first consumer was RepeatMasker — so a run made without
`-l/--library` silently lacked it and could not be compared against anything.
It is now written by the clustering step, where the per-TRC dimer files it
concatenates are produced.

Also pins the ordering: `glob` order is filesystem-dependent (verified: on a
real 73-TRC directory it returns TRC_1, TRC_7, TRC_11, TRC_42 …), which made
the file's bytes differ between machines for identical input.

Run: python3 tests/test_consensus_pool.py
"""
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import tc_utils as tc

failures = []


def check(label, got, want):
    ok = got == want
    print(F"  {'ok  ' if ok else 'FAIL'} {label}: got {got!r}, want {want!r}")
    if not ok:
        failures.append(label)


def headers(path):
    return [l[1:].strip() for l in open(path) if l.startswith(">")]


tmp = tempfile.mkdtemp(prefix="tc_pool_test_")
try:
    d = os.path.join(tmp, "tc_consensus")
    os.makedirs(d)
    # deliberately created out of order, and spanning one/two digits so a
    # lexicographic sort would put TRC_10 before TRC_2
    for i in (10, 2, 1, 3):
        with open(os.path.join(d, F"TRC_{i}_dimers.fasta"), "w") as fh:
            fh.write(F">TRC_{i}_rep0\nACGT\n>TRC_{i}_rep1\nTTTT\n")
    # a non-dimer file that must not be swept in
    with open(os.path.join(d, "TRC_1.fasta"), "w") as fh:
        fh.write(">TRC_1_monomer\nAC\n")

    print("writes the pool in natural TRC order")
    out = tc.write_consensus_sequences_all(d)
    check("returns the canonical path", os.path.basename(out),
          "consensus_sequences_all.fasta")
    check("natural TRC order (not lexicographic, not glob)", headers(out),
          ["TRC_1_rep0", "TRC_1_rep1", "TRC_2_rep0", "TRC_2_rep1",
           "TRC_3_rep0", "TRC_3_rep1", "TRC_10_rep0", "TRC_10_rep1"])
    check("only *_dimers.fasta is concatenated",
          any("monomer" in h for h in headers(out)), False)

    print("force=False leaves an existing file alone")
    with open(out, "w") as fh:
        fh.write(">hand_written\nAAAA\n")
    tc.write_consensus_sequences_all(d)
    check("existing file untouched", headers(out), ["hand_written"])
    tc.write_consensus_sequences_all(d, force=True)
    check("force=True rebuilds it", len(headers(out)), 8)

    print("tolerates a directory with nothing to concatenate")
    empty = os.path.join(tmp, "empty_consensus")
    os.makedirs(empty)
    check("no dimer files -> None", tc.write_consensus_sequences_all(empty), None)
    check("and writes nothing",
          os.path.exists(os.path.join(empty, "consensus_sequences_all.fasta")), False)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(F"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("test_consensus_pool.py: all checks passed")
