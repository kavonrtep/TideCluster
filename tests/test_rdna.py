#!/usr/bin/env python3
"""Unit tests for rDNA (45S/5S) identification.

Covers the pure logic without invoking real blast: subprocess.run is
monkeypatched to (a) record the makeblastdb subject and (b) synthesise a blastn
outfmt-6 file from that subject's seqids. Exercises:
  - blastn_rdna_reference_coverage: reference-coverage parsing + identity gate
  - assign_rdna_to_trcs: per-TRC aggregation, threshold, type selection
  - _write_rdna_attributes: rDNA_type/rDNA_coverage added only to matched TRCs
  - identify_rdna: consensus-first search + genomic fallback orchestration

Run: python3 tests/test_rdna.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


# rDNA library used by all tests: a 100 bp 18S (45S) and a 50 bp 5S.
RDNA_LIB_TEXT = (">r45#rDNA_45S/18S\n" + "A" * 100 + "\n"
                 ">r5#rDNA_5S/5S\n" + "C" * 50 + "\n")


def make_fake_run():
    """Fake subprocess.run: makeblastdb records the subject; blastn writes a
    deterministic outfmt-6 derived from the subject seqids.

    Rule: a TRC_1* subject gets a full-length 18S hit (cov 1.0 @ 90% id); a
    TRC_2* subject gets a half-length 5S hit (cov 0.5 @ 90%) plus a weak 18S hit
    (@ 80%, below the default identity gate); anything else gets nothing.
    """
    state = {"subject": None}

    def fake_run(cmd, *a, **k):
        if cmd.startswith("makeblastdb"):
            toks = cmd.split()
            state["subject"] = toks[toks.index("-in") + 1]
            return None
        if cmd.startswith("blastn"):
            out = cmd.rsplit(">", 1)[1].strip()
            rows = []
            for name, _seq in tc.fasta_to_list(state["subject"]):
                trc_n = name.split("_")[1]
                if trc_n == "1":
                    rows.append(("r45#rDNA_45S/18S", name, "90.0", "1", "100"))
                elif trc_n == "2":
                    rows.append(("r5#rDNA_5S/5S", name, "90.0", "1", "25"))
                    rows.append(("r45#rDNA_45S/18S", name, "80.0", "1", "100"))
            with open(out, "w") as f:
                for r in rows:
                    f.write("\t".join(r) + "\n")
            return None
        return None

    return fake_run


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    d = tempfile.mkdtemp()
    lib = os.path.join(d, "rdna_library.fasta")
    write(lib, RDNA_LIB_TEXT)

    orig_run = tc.subprocess.run
    tc.subprocess.run = make_fake_run()
    try:
        # --- blastn_rdna_reference_coverage -----------------------------------
        subj = os.path.join(d, "subject.fasta")
        write(subj, ">TRC_1_c0\nACGT\n>TRC_2_c0\nACGT\n")
        best = tc.blastn_rdna_reference_coverage(lib, subj, cpu=1, min_identity=85.0)
        check("TRC_1 consensus -> full 45S coverage",
              best.get("TRC_1_c0", {}).get("rDNA_45S", (0,))[0] == 1.0)
        check("TRC_2 consensus -> half 5S coverage",
              abs(best.get("TRC_2_c0", {}).get("rDNA_5S", (0,))[0] - 0.5) < 1e-9)
        check("weak (80%) 45S hit excluded by identity gate",
              "rDNA_45S" not in best.get("TRC_2_c0", {}))

        # --- assign_rdna_to_trcs ---------------------------------------------
        calls = tc.assign_rdna_to_trcs(best, min_coverage=0.7)
        check("threshold 0.7: only TRC_1 called 45S",
              calls == {"TRC_1": ("45S", 1.0)})
        calls_lo = tc.assign_rdna_to_trcs(best, min_coverage=0.4)
        check("threshold 0.4: TRC_2 called 5S at 0.5",
              calls_lo.get("TRC_2") == ("5S", 0.5))
    finally:
        tc.subprocess.run = orig_run

    # --- _write_rdna_attributes ----------------------------------------------
    gff = os.path.join(d, "clustering.gff3")
    write(gff, "##gff-version 3\n" + "".join(
        f"chr1\tTideCluster\ttandem_repeat\t{s}\t{e}\t1\t.\t.\tName={n};repeat_type=TR\n"
        for n, s, e in [("TRC_1", 1, 100), ("TRC_2", 200, 300), ("TRC_3", 400, 500)]))
    tc._write_rdna_attributes(gff, {"TRC_1": ("45S", 1.0)})
    body = open(gff).read()
    check("TRC_1 gets rDNA_type + rDNA_coverage",
          "Name=TRC_1;repeat_type=TR;rDNA_type=45S;rDNA_coverage=1.0" in body)
    check("TRC_3 (unmatched) left unchanged",
          "Name=TRC_3;repeat_type=TR\n" in body and "TRC_3;repeat_type=TR;rDNA" not in body)

    # --- identify_rdna end-to-end (consensus + genomic fallback) -------------
    prefix = os.path.join(d, "run")
    write(prefix + "_clustering.gff3", "##gff-version 3\n" + "".join(
        f"chr1\tTideCluster\ttandem_repeat\t{s}\t{e}\t1\t.\t.\tName={n};repeat_type=TR\n"
        for n, s, e in [("TRC_1", 1, 100), ("TRC_2", 200, 300), ("TRC_3", 400, 500)]))
    cons = prefix + "_consensus"
    os.makedirs(cons, exist_ok=True)
    # TRC_1, TRC_2 have consensus; TRC_3 does not -> genomic fallback path
    write(os.path.join(cons, "TRC_1_dimers.fasta"), ">TRC_1_rep0_chr1_1\nACGT\n")
    write(os.path.join(cons, "TRC_2_dimers.fasta"), ">TRC_2_rep0_chr1_1\nACGT\n")
    genome = os.path.join(d, "genome.fasta")
    write(genome, ">chr1\n" + "ACGT" * 250 + "\n")  # 1000 bp

    tc.subprocess.run = make_fake_run()
    try:
        calls = tc.identify_rdna(prefix, genome, lib, cpu=1,
                                 min_coverage=0.7, min_identity=85.0)
    finally:
        tc.subprocess.run = orig_run
    check("identify_rdna labels TRC_1 as 45S", calls.get("TRC_1") == ("45S", 1.0))
    check("identify_rdna leaves TRC_2 below-threshold unlabelled",
          "TRC_2" not in calls)
    gbody = open(prefix + "_clustering.gff3").read()
    check("clustering.gff3 updated in place with rDNA_type",
          "Name=TRC_1;repeat_type=TR;rDNA_type=45S" in gbody)
    check("_rdna.tsv side table written",
          os.path.exists(prefix + "_rdna.tsv")
          and "TRC_1\t45S\t1.0" in open(prefix + "_rdna.tsv").read())
    # TRC_3 had no consensus -> went through the genomic-fallback subject path
    check("genomic fallback ran without error for consensus-less TRC_3",
          "TRC_3" not in calls)  # fake yields no rDNA hit for TRC_3 region

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
