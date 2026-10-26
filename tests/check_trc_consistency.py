#!/usr/bin/env python3
"""Check that a run's per-TRC outputs name only TRCs in its clustering GFF3.

<prefix>_clustering.gff3 is the authority on which TRCs exist. Up to 1.22.2,
overlap resolution dropped fully absorbed TRCs from it but their consensus
files survived, and the superfamily fallback and the annotation table brought
them back (issue #12, #8). This asserts, for whichever of these exist:

  <prefix>_consensus/TRC_*_dimers.fasta   same TRC set as the GFF3
  <prefix>_trc_superfamilies.csv          subset of the GFF3
  <prefix>_annotation.tsv                 subset of the GFF3

Run: python3 tests/check_trc_consistency.py <prefix>
Exit status 1 on any mismatch.
"""
import csv
import glob
import os
import re
import sys


def gff3_trcs(path):
    names = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            m = re.search(r"(?:^|;)Name=(TRC_\d+)", line.rstrip("\n").split("\t")[-1])
            if m:
                names.add(m.group(1))
    return names


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    prefix = sys.argv[1]
    gff = prefix + "_clustering.gff3"
    if not os.path.exists(gff):
        sys.exit(f"FAIL: {gff} not found")
    in_gff = gff3_trcs(gff)
    failures = []

    def report(label, extra, missing=()):
        if extra or missing:
            msg = f"FAIL {label}:"
            if extra:
                msg += f" not in GFF3: {', '.join(sorted(extra))}"
            if missing:
                msg += f" missing (in GFF3): {', '.join(sorted(missing))}"
            print(msg)
            failures.append(label)
        else:
            print(f"ok   {label}")

    cdir = prefix + "_consensus"
    if os.path.isdir(cdir):
        cons = {os.path.basename(p)[:-len("_dimers.fasta")]
                for p in glob.glob(os.path.join(cdir, "TRC_*_dimers.fasta"))}
        report(f"{cdir}/ == GFF3 TRCs ({len(cons)} vs {len(in_gff)})",
               cons - in_gff, in_gff - cons)

    sf = prefix + "_trc_superfamilies.csv"
    if os.path.exists(sf):
        with open(sf) as f:
            members = {r["TRC"] for r in csv.DictReader(f)}
        report(f"{os.path.basename(sf)} members ⊆ GFF3", members - in_gff)

    ann = prefix + "_annotation.tsv"
    if os.path.exists(ann):
        with open(ann) as f:
            rows = {r["TRC"] for r in csv.DictReader(f, delimiter="\t")}
        report(f"{os.path.basename(ann)} rows ⊆ GFF3", rows - in_gff)

    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
