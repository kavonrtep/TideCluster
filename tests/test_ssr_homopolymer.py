#!/usr/bin/env python3
"""Homopolymer SSR arrays: no csv field-limit crash, monomer_size = motif (#13).

kitehor ssr-scan writes `consensus_single` over KITE's period. KITE cannot
report a 1 bp period, so on a long poly-A array its rank-1 peak is a spurious
50-150 kb period and the consensus cell runs to hundreds of kb. Up to 1.22.3
this crashed `build_monomer_size_csv` (csv's 131,072-char field limit), and
once read, `monomer_size` kept the spurious period and the ⚠short flag fired.

Fixture tests/data/homopolymer_kite/:
  TRC_1    3 synthetic 355 kb poly-A arrays (tmp/issue13 var_tbias), kitehor
           0.13.2. The consensus cells are stored as @POLYA@ and expanded here
           to 3 x 70,000 chars.
  TRC_170  drapa ACCCACCCC SSR: chr4 has a lone 2580 bp kite peak (changes to
           9, the 2580 stays as alt_longer_period); chr8 already ranks the
           9 bp motif first (must be unchanged).

Run: python3 tests/test_ssr_homopolymer.py
"""
import csv
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402

FX = os.path.join(ROOT, "tests", "data", "homopolymer_kite")
POLYA = ";".join(["a" * 70000] * 3)


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="tc_ssr_homopolymer_")
    ssr_tsv = os.path.join(work, "kitehor.ssr.tsv")
    with open(os.path.join(FX, "kitehor.ssr.tsv")) as f, open(ssr_tsv, "w") as out:
        out.write(f.read().replace("@POLYA@", POLYA))
    check("fixture has a cell above csv's default field limit",
          len(POLYA) > 131072)

    out_csv = os.path.join(work, "monomer_size_top3_estimats.csv")
    try:
        tc.build_monomer_size_csv(
            kite_tsv=os.path.join(FX, "kitehor.kite.tsv"),
            ssr_tsv=ssr_tsv,
            rescored_peaks_tsv=os.path.join(FX, "kitehor.rescored.peaks.tsv"),
            tandem_validate_tsv=os.path.join(FX, "kitehor.tandem_validate.tsv"),
            out_csv=out_csv,
            trc_repeat_type=tc.parse_trc_ssr_motif_len(
                os.path.join(FX, "clustering.gff3")))
        ok = True
    except csv.Error as err:
        print(err)
        ok = False
    check("build_monomer_size_csv reads a 210 k-char ssr.tsv cell", ok)
    if not ok:
        sys.exit(1)

    rows = list(csv.DictReader(open(out_csv), delimiter="\t"))
    by = {(r["TRC_ID"], r["seqid"]): r for r in rows}
    polya = [r for r in rows if r["TRC_ID"] == "TRC_1"]
    check("3 poly-A rows", len(polya) == 3)
    check("poly-A founder_period = 1",
          all(r["founder_period"] == "1" for r in polya))
    check("poly-A monomer_size = 1 (was the 81-105 kb kite peak)",
          all(r["monomer_size"] == "1" for r in polya))
    check("poly-A score blank (rank-1 peak is not the motif)",
          all(r["score"] == "" for r in polya))
    check("poly-A not flagged weak-short",
          all(r["weak_short_founder_flag"] == "false" for r in polya))

    chr4 = by.get(("TRC_170", "chr4"), {})
    check("TRC_170 chr4 monomer_size = motif 9", chr4.get("monomer_size") == "9")
    check("TRC_170 chr4 keeps 2580 as alt_longer_period",
          chr4.get("alt_longer_period") == "2580")
    check("TRC_170 chr4 not flagged weak-short",
          chr4.get("weak_short_founder_flag") == "false")
    chr8 = by.get(("TRC_170", "chr8"), {})
    check("TRC_170 chr8 (motif already rank 1) keeps its score",
          chr8.get("monomer_size") == "9"
          and chr8.get("score") == "0.4165879875"
          and chr8.get("monomer_size_2") == "4097")

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
