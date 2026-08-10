#!/usr/bin/env python3
"""Regression test for the split_fasta_to_chunk_files "Too many open files" bug.

`split_fasta_to_chunk_files` used to open one handle per chunk file at once
(`{p: open(p, "w") for p in file_paths}`). The number of chunk files grows with
genome size (~genome_bp / chunk_size, ~1800 for a 90 Gbp genome at the default
chunk_size), so on large genomes it exceeds the open-file limit and the split —
and therefore the whole tc_reannotate / chunked-RepeatMasker path — aborts with
`OSError: [Errno 24] Too many open files`.

The fix opens chunk files lazily through a bounded LRU cache (<= max_open_handles
descriptors). This test forces MORE chunk files than the process is allowed to
have open at once, and checks that:
  1. the split completes (no OSError) with the fd limit below the file count, and
  2. every planned piece is written intact (LRU eviction + append-reopen does not
     corrupt, truncate, or misplace any chunk).

Run: python3 tests/test_chunk_fd_limit.py
"""
import os
import resource
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def read_chunk_records(paths):
    """Return {token: sequence} parsed from all chunk FASTA files."""
    records = {}
    for p in paths:
        with open(p) as fh:
            token, buf = None, []
            for line in fh:
                if line.startswith(">"):
                    if token is not None:
                        records[token] = "".join(buf)
                    token, buf = line[1:].strip(), []
                else:
                    buf.append(line.strip())
            if token is not None:
                records[token] = "".join(buf)
    return records


def run_constrained(limit):
    resource.setrlimit(resource.RLIMIT_NOFILE, (limit, limit))
    import tc_utils as tc

    tmpdir = tempfile.mkdtemp(prefix="tc_fd_test_")
    genome = os.path.join(tmpdir, "genome.fasta")
    # A base long enough to be cut into far more pieces than `limit`.
    seq = "".join("ACGT"[(i * 7) % 4] for i in range(50000))  # 50 kb, non-trivial
    with open(genome, "w") as fh:
        fh.write(">chr1\n%s\n" % seq)

    out_dir = os.path.join(tmpdir, "chunks")
    os.makedirs(out_dir)
    # chunk_size=100 -> ~500 pieces -> ~500 chunk files, well above `limit`.
    file_paths, matching_table, token_to_file = tc.split_fasta_to_chunk_files(
        genome, out_dir, chunk_size=100, overlap=10
    )

    assert len(file_paths) > limit, (
        "test misconfigured: need more chunk files (%d) than the fd limit (%d) "
        "to prove bounded handles" % (len(file_paths), limit)
    )
    # Content integrity: each planned piece must be present exactly.
    records = read_chunk_records(file_paths)
    for header, _i, start, end, token in matching_table:
        assert token in records, "missing chunk record for token %s" % token
        assert records[token] == seq[start:end], (
            "chunk %s content mismatch (LRU eviction corrupted a write)" % token
        )
    assert len(records) == len(matching_table), "unexpected number of chunk records"
    print(
        "  constrained(limit=%d): %d chunk files, %d pieces, all intact"
        % (limit, len(file_paths), len(matching_table))
    )


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "--child":
        run_constrained(int(sys.argv[2]))
        sys.exit(0)
    # Run under a lowered hard limit in a child so it cannot be raised away.
    import subprocess

    subprocess.run(
        [sys.executable, os.path.abspath(__file__), "--child", "350"], check=True
    )
    print("test_chunk_fd_limit: PASSED")
