#!/usr/bin/env python3
"""Unit test for the chunked, pooled RepeatMasker path (FR-1 / P1).

Validates the coordinate machinery WITHOUT invoking real RepeatMasker:

  A) split_fasta_to_chunk_files: long sequences are cut into overlapping
     pieces of ~chunk_size with the right overlap, short sequences emitted
     whole, each piece gets a short numeric header, the matching_table offsets
     map chunk-local coordinates back to genome coordinates, and each chunk
     file holds exactly the planned subsequence.

  B) run_repeatmasker_genome_chunked: with a fake RepeatMasker (monkeypatched
     subprocess.run, inherited by the fork-based Pool on Linux) every chunk hit
     is remapped to the correct genome coordinate and original sequence name,
     and the GFF3 is written sorted (pool-order independent).

  C) hit-less chunks (RepeatMasker writes no .out) contribute nothing instead
     of crashing.

Run: python3 tests/test_chunked_repeatmasker.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import tc_utils as tc  # noqa: E402


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for name, seq in seqs:
            fh.write(F">{name}\n{seq}\n")


# --- fake RepeatMasker -------------------------------------------------------
# Parses the RepeatMasker argv built by _repeatmasker_chunk_worker, reads the
# chunk FASTA, and writes a .out with one hit per sequence at local [2, 4].
def _fake_rm_run(cmds, *a, **k):
    out_dir = cmds[cmds.index("-dir") + 1]
    chunk_fasta = cmds[-1]
    tokens = [h for h in tc.read_fasta_sequence_size(chunk_fasta)]
    out_file = os.path.join(out_dir, os.path.basename(chunk_fasta) + ".out")
    with open(out_file, "w") as fh:
        fh.write("h1\nh2\nh3\n")  # 3 RepeatMasker header lines
        for tok in tokens:
            # cols: score div del ins QUERY begin end (left) strand match NAME ...
            fh.write(F"300 1.0 0.0 0.0 {tok} 2 4 (0) + REP TRCsat 1 3 (0) 1\n")
    return None


def _fake_rm_run_no_out(cmds, *a, **k):
    return None  # writes nothing -> simulates a hit-less / empty chunk


def main():
    failures = []

    def check(name, cond):
        print(("PASS " if cond else "FAIL ") + name)
        if not cond:
            failures.append(name)

    d = tempfile.mkdtemp()
    genome_seq = "".join("ACGT"[i % 4] for i in range(100))
    genome = os.path.join(d, "genome.fasta")
    write_fasta(genome, [("seqA", genome_seq), ("seqB", "ACGTACGTAC")])  # 100 + 10

    # --- A) splitter ---------------------------------------------------------
    out_dir = os.path.join(d, "chunks")
    os.makedirs(out_dir, exist_ok=True)
    file_paths, mt = tc.split_fasta_to_chunk_files(
        genome, out_dir, chunk_size=20, overlap=5
    )
    expected_mt = [
        ["seqA", 0, 0, 25, "0"],
        ["seqA", 1, 20, 45, "1"],
        ["seqA", 2, 40, 65, "2"],
        ["seqA", 3, 60, 85, "3"],
        ["seqA", 4, 80, 100, "4"],
        ["seqB", 0, 0, 10, "5"],
    ]
    check("matching_table tiles genome with overlap", mt == expected_mt)
    check("all chunk files exist", all(os.path.exists(p) for p in file_paths))

    # each token's chunk subsequence == genome[start:end]
    token_seq = {}
    for p in file_paths:
        with open(p) as fh:
            token_seq.update(tc.read_single_fasta_to_dictionary(fh))
    seqmap = {"seqA": genome_seq, "seqB": "ACGTACGTAC"}
    subseq_ok = all(
        token_seq[row[4]] == seqmap[row[0]][row[2]:row[3]] for row in mt
    )
    check("chunk subsequences match genome slices", subseq_ok)

    # --- B) chunked RepeatMasker remap ---------------------------------------
    orig_run = tc.subprocess.run
    tc.subprocess.run = _fake_rm_run
    out_gff3 = os.path.join(d, "rm.gff3")
    try:
        tc.run_repeatmasker_genome_chunked(
            genome, os.path.join(d, "lib.fasta"), 2, "", out_gff3,
            chunk_size=20, overlap=5,
        )
    finally:
        tc.subprocess.run = orig_run

    rows = []
    with open(out_gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = tc.Gff3Feature(line)
            rows.append((f.seqid, f.start, f.end, f.strand,
                         f.attributes_dict["Name"]))
    # local hit [2,4] + offset(token) -> genome coords; offsets 0,20,40,60,80 / 0
    expected = [
        ("seqA", 2, 4, "-", "TRCsat"),
        ("seqA", 22, 24, "-", "TRCsat"),
        ("seqA", 42, 44, "-", "TRCsat"),
        ("seqA", 62, 64, "-", "TRCsat"),
        ("seqA", 82, 84, "-", "TRCsat"),
        ("seqB", 2, 4, "-", "TRCsat"),
    ]
    check("hits remapped to genome coordinates + names", set(rows) == set(expected))
    check("GFF3 written sorted (pool-order independent)",
          rows == sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3], r[4])))
    check("Cluster_ID and Name both present",
          all(";Name=" in l and "Cluster_ID=" in l
              for l in open(out_gff3) if not l.startswith("#")))

    # --- C) hit-less chunks (no .out) ----------------------------------------
    tc.subprocess.run = _fake_rm_run_no_out
    out_gff3_empty = os.path.join(d, "rm_empty.gff3")
    try:
        tc.run_repeatmasker_genome_chunked(
            genome, os.path.join(d, "lib.fasta"), 2, "", out_gff3_empty,
            chunk_size=20, overlap=5,
        )
    except Exception as e:  # noqa: BLE001
        check("hit-less chunks do not crash", False)
        print("  raised:", repr(e))
    else:
        check("hit-less chunks do not crash", True)
        data_lines = [l for l in open(out_gff3_empty) if not l.startswith("#")]
        check("no .out -> empty (header-only) GFF3", data_lines == [])
    finally:
        tc.subprocess.run = orig_run

    print("\nALL PASS" if not failures else f"\nFAILED: {failures}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
