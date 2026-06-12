#!/usr/bin/env python3
"""Regression guard: run_repeatmasker_with_renaming() must not crash when
RepeatMasker produces no .out file.

RepeatMasker writes no `<input>.out` when the input FASTA is empty or no repeats
are detected (e.g. annotating omitted short regions on a contig with no library
hits). The old code called restore_sequence_names_in_repeatmasker_output() on the
missing file and died with FileNotFoundError, aborting run_all before TAREAN
(observed on Pisum contig ptg000115l_rc). The guard must instead yield an empty
.out that get_repeatmasker_annotation() parses to {} (zero annotations).

Self-contained: monkeypatches subprocess.run so no real RepeatMasker is invoked
and no .out is produced. Run: python3 tests/test_repeatmasker_empty_out.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    d = tempfile.mkdtemp()
    in_fa = os.path.join(d, "consensus.fasta")
    with open(in_fa, "w") as fh:
        fh.write(">TRC_1_rep0\nACGTACGTACGT\n")
    lib = os.path.join(d, "lib.fasta")
    with open(lib, "w") as fh:
        fh.write(">ref#Satellite/Fake\nTTTTAAAACCCC\n")

    # Simulate RepeatMasker finding nothing: a no-op that creates no .out file.
    orig_run = tc.subprocess.run
    tc.subprocess.run = lambda *a, **k: None
    try:
        rm_file = tc.run_repeatmasker_with_renaming(in_fa, lib, cpu=1)
    except Exception as e:  # noqa: BLE001
        rm_file = None
        check("no crash when RepeatMasker writes no .out", False)
        print("  raised:", repr(e))
    finally:
        tc.subprocess.run = orig_run

    if rm_file is not None:
        check("no crash when RepeatMasker writes no .out", True)
        check("returns an existing .out path", os.path.isfile(rm_file))
        # downstream annotation parse must yield zero annotations, not error
        ann = tc.get_repeatmasker_annotation(rm_file, {"TRC_1_rep0": 12},
                                             os.path.join(d, "tc"), parse_id=False)
        check("empty .out -> zero annotations", ann == {})
        # the renamed temp FASTA must be cleaned up
        check("renamed temp FASTA cleaned up",
              not os.path.exists(in_fa + "_renamed.fasta"))

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
