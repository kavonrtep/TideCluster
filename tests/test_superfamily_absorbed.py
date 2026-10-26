#!/usr/bin/env python3
"""Superfamily fallback must not attach TRCs missing from the clustering GFF3.

Issue #12: overlap resolution dropped fully absorbed TRCs from
<prefix>_clustering.gff3 but (up to 1.22.2) left their consensus files, and
the below-threshold fallback in compare_trc_by_blast.R then listed them as
superfamily members with no arrays and no report page.

Fixture tests/data/absorbed_trc/ is cut from the Scutellaria galericulata run
in the issue: the TAREAN dimer library of TRC_5 and TRC_15, the leftover dimer
files of their absorbed variants TRC_25/32/41 (not in the GFF3), and TRC_99, a
copy of TRC_25's dimers that IS in the GFF3 (positive control: the filter must
not stop a real below-threshold TRC from being attached).

Unfixed, the CSV lists SF {TRC_5, TRC_25, TRC_41} and {TRC_15, TRC_32}.

Run: python3 tests/test_superfamily_absorbed.py   (needs Rscript + blastn)
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURE = os.path.join(ROOT, "tests", "data", "absorbed_trc")
SCRIPT = os.path.join(ROOT, "tarean", "compare_trc_by_blast.R")


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="tc_sf_absorbed_")
    try:
        run = os.path.join(work, "run")
        shutil.copytree(FIXTURE, run)
        proc = subprocess.run(
            [SCRIPT, "-i", "tc_consensus_dimer_library.fasta", "-p", "tc",
             "-t", "2", "--consensus_dir", "tc_consensus",
             "--score_threshold", "20"],
            cwd=run, capture_output=True, text=True)
        check("compare_trc_by_blast.R exits 0", proc.returncode == 0)
        if proc.returncode != 0:
            print(proc.stdout[-2000:], proc.stderr[-2000:])
        with open(os.path.join(run, "tc_trc_superfamilies.csv")) as f:
            rows = list(csv.DictReader(f))
        members = {r["TRC"]: r for r in rows}
        ghosts = {"TRC_25", "TRC_32", "TRC_41"}
        check("absorbed TRCs not listed as members",
              not ghosts & set(members))
        check("in-GFF3 below-threshold TRC_99 still attached by the fallback",
              members.get("TRC_99", {}).get("fallback") == "TRUE")
        check("TRC_99 joins TRC_5 (same sequence as TRC_25 did)",
              "TRC_5" in members and members.get("TRC_99", {}).get(
                  "Superfamily") == members["TRC_5"]["Superfamily"])
        skip_line = [l for l in proc.stderr.splitlines()
                     if l.startswith("Fallback: skipping")]
        check("skipped TRCs are named in the log",
              len(skip_line) == 1 and all(t in skip_line[0] for t in ghosts))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
