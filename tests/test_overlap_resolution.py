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
    absorbed = tc.resolve_trc_overlaps(gin, gout)
    after = parse(gout)
    check("returns the absorbed TRCs", absorbed == ["TRC_C"])

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
    absorbed2 = tc.resolve_trc_overlaps(gout, gout2)
    check("idempotent (second pass identical geometry)",
          parse(gout2) == after)
    check("second pass absorbs nothing", absorbed2 == [])

    # Consensus files (issue #12). clustering() drops the absorbed TRCs from
    # its consensus set before saving, so <prefix>_consensus/ has the same TRCs
    # as the resolved GFF3. A stale TRC file from an earlier run into the same
    # directory must go too, while non-TRC files there are left alone.
    cdir = os.path.join(d, "tc_consensus")
    os.makedirs(cdir)
    for stale in ("TRC_9.fasta", "TRC_9_dimers.fasta"):
        open(os.path.join(cdir, stale), "w").write(">old\nACGT\n")
    keep = os.path.join(cdir, "consensus_sequences_all.fasta_renamed.fasta.out")
    open(keep, "w").write("rm by-product\n")
    cons = {n: {f"{n}_rep": "ACGT"} for n in ("TRC_A", "TRC_B", "TRC_C")}
    dim = {n: {f"{n}_rep": "ACGTACGT"} for n in cons}
    for trc in absorbed:
        cons.pop(trc, None)
        dim.pop(trc, None)
    tc.save_consensus_files(cdir, cons, dim)
    on_disk = {f[:-len("_dimers.fasta")] for f in os.listdir(cdir)
               if f.endswith("_dimers.fasta")}
    check("consensus dir == TRCs in resolved GFF3", on_disk == set(by_name))
    check("stale TRC_9 files removed",
          not any(f.startswith("TRC_9") for f in os.listdir(cdir)))
    check("non-TRC files in consensus dir kept", os.path.exists(keep))

    # clustering() must use the return value (a guard on the wiring, which
    # needs a full clustering run to exercise end to end).
    src = open(os.path.join(ROOT, "TideCluster.py")).read()
    check("clustering() prunes cons_cls by the absorbed TRCs",
          "absorbed = tc.resolve_trc_overlaps(" in src
          and "cons_cls.pop(trc, None)" in src
          and "cons_cls_dimer.pop(trc, None)" in src)

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
