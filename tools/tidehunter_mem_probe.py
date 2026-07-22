#!/usr/bin/env python3
"""
tidehunter_mem_probe.py — measure TideHunter peak RSS + wall time per "part",
per monomer-period round.

Purpose
-------
Fix 5 in docs/large_genome_scaling_notes.md proposes replacing the serial
``run_tidehunter`` part loop (~1,800 serial invocations at 90 Gbp) with a
memory-gated pool of single-threaded workers. The governing constraint is
memory: N concurrent single-threaded TideHunter processes cost ~N × the
per-part working set, so the pool size must be
``floor($AGENT_MEMORY / part_peak)`` — and ``part_peak`` must be MEASURED,
not guessed.

This probe measures ``part_peak`` (peak resident set, from the kernel's
monotonic ``/proc/<pid>/status:VmHWM``) and wall time for a real "part"
(a pack of 500 kb chunks, exactly as ``split_fasta_to_chunks`` produces),
across:
  * the three ``tidehunter_long`` rounds — round 1 (-p 40 -P 3000),
    round 2 (-p 3001 -P 10000), round 3 (-p 10001 -P 25000) — which have
    very different period regimes and are expected to peak very
    differently; the pool must be sized on the WORST round, and
  * several part sizes, so we can see whether peak RSS scales with input
    size (⇒ shrinking parts helps) or is dominated by a fixed
    per-array/base cost (⇒ shrinking parts won't help and dynamic
    admission is needed).

It does NOT touch the pipeline. It reuses ``tc_utils.split_fasta_to_chunks``
so the part content is faithful to what the real pipeline feeds TideHunter.

Usage
-----
    conda activate tidecluster
    ./tools/tidehunter_mem_probe.py \
        --fasta test_data/CEN6_ver_220406.fasta \
        --part-sizes 5,15,30 --threads 1

Notes
-----
* Peak RSS is read from ``/proc/<pid>/status:VmHWM`` (KiB, kernel
  high-water mark, monotonic), polled every ~10 ms. Because VmHWM only
  increases, polling captures the true peak for any run longer than a few
  poll intervals; sub-100 ms runs (which are small anyway) may undersample.
* ``--threads`` defaults to 1 — that is the per-worker footprint the pool
  formula needs. Pass e.g. ``--threads $AGENT_CPUS`` to instead measure
  today's serial-path peak (one process, all threads) for contrast.
* The sandbox has no ``/usr/bin/time``; timing uses
  ``time.perf_counter_ns`` and arithmetic is done in Python (no ``bc``).
"""

import argparse
import os
import re
import shlex
import subprocess
import sys
import tempfile
import threading
import time

# import the pipeline's own chunker so the "part" is faithful
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tc_utils as tc  # noqa: E402

# The three tidehunter_long rounds (TideCluster.py:947-949) + the default
# single-round args (TideCluster.py:1107) which is identical to round 1.
ROUND_ARGS = {
    "round1_p40-3000":     "-p 40 -P 3000 -c 5 -e 0.25",
    "round2_p3001-10000":  "-p 3001 -P 10000 -c 5 -e 0.25",
    "round3_p10001-25000": "-p 10001 -P 25000 -c 5 -e 0.25",
}

CHUNK_SIZE = 500000   # TideCluster.py:792 / :940
OVERLAP = 50000       # TideCluster.py:793 / :941


def check_tidehunter():
    try:
        v = subprocess.check_output(["TideHunter", "-v"],
                                    stderr=subprocess.STDOUT).decode().strip()
    except (OSError, subprocess.CalledProcessError) as e:
        sys.exit(F"TideHunter not runnable on PATH: {e}\n"
                 F"Did you `conda activate tidecluster`?")
    print(F"TideHunter version: {v}")
    return v


def build_part(chunked_fasta, target_bytes, out_path):
    """Pack whole 500 kb chunk-sequences from ``chunked_fasta`` into
    ``out_path`` until at least ``target_bytes`` of sequence is written.
    Returns the number of sequence bases written (input size of the part)."""
    written = 0
    n_seqs = 0
    with open(chunked_fasta) as fh, open(out_path, "w") as out:
        for header, seq in tc.read_single_fasta_as_generator(fh):
            out.write(F">{header}\n{seq}\n")
            written += len(seq)
            n_seqs += 1
            if written >= target_bytes:
                break
    return written, n_seqs


def run_and_measure(part_fasta, th_args, threads, out_file, timeout):
    """Run one TideHunter invocation and measure its peak RSS.

    The authoritative peak is ``os.wait4``'s ``rusage.ru_maxrss`` — the
    kernel's exact high-water RSS for the child, with NO sampling risk (an
    earlier ``/proc/<pid>/status:VmHWM`` poller under/over-counted by up to
    ~1.8x because a late allocation spike fell between polls / after the
    /proc entry vanished). VmHWM sampling and TideHunter's own
    ``Peak RSS: X GB`` stderr line are kept only as cross-checks.

    Returns dict with ru_maxrss_kib, vmhwm_kib, th_self_mb, wall_s,
    returncode, n_features, stderr_tail.
    """
    cmd = (["TideHunter", "-f", "2", "-o", out_file, "-t", str(threads)]
           + shlex.split(th_args) + [part_fasta])
    devnull = os.open(os.devnull, os.O_WRONLY)
    err_r, err_w = os.pipe()

    t0 = time.perf_counter_ns()
    pid = os.fork()
    if pid == 0:  # child
        try:
            os.dup2(devnull, 1)
            os.dup2(err_w, 2)
            os.close(err_r)
            os.execvp("TideHunter", cmd)
        except Exception:  # pragma: no cover
            pass
        os._exit(127)

    # parent
    os.close(err_w)
    os.close(devnull)
    err_chunks = []

    def drain():
        while True:
            b = os.read(err_r, 65536)
            if not b:
                break
            err_chunks.append(b)

    tdr = threading.Thread(target=drain)
    tdr.start()

    vmhwm_kib = 0
    status_path = F"/proc/{pid}/status"
    deadline = t0 + int(timeout * 1e9)
    ru = None
    wstatus = 0
    timed_out = False
    while True:
        wpid, wstatus, rusage = os.wait4(pid, os.WNOHANG)
        if wpid == pid:
            ru = rusage
            break
        try:  # VmHWM cross-check only
            with open(status_path) as sf:
                for line in sf:
                    if line.startswith("VmHWM:"):
                        vmhwm_kib = max(vmhwm_kib, int(line.split()[1]))
                        break
        except (OSError, ValueError):
            pass
        if time.perf_counter_ns() > deadline:
            os.kill(pid, 9)
            _wp, wstatus, ru = os.wait4(pid, 0)
            timed_out = True
            break
        time.sleep(0.01)
    wall_s = (time.perf_counter_ns() - t0) / 1e9

    tdr.join()
    os.close(err_r)
    err = b"".join(err_chunks).decode(errors="replace")

    ru_maxrss_kib = ru.ru_maxrss if ru else -1
    ret = -9 if timed_out else os.waitstatus_to_exitcode(wstatus)
    m = re.search(r"Peak RSS:\s*([\d.]+)\s*GB", err)
    th_self_mb = float(m.group(1)) * 1024 if m else -1.0
    err_tail = "TIMEOUT" if timed_out else (
        " | ".join(err.strip().splitlines()[-2:]) if err.strip() else "")
    n_feat = -1
    if os.path.exists(out_file):
        with open(out_file) as f:
            n_feat = sum(1 for _ in f)
    return dict(ru_maxrss_kib=ru_maxrss_kib, vmhwm_kib=vmhwm_kib,
                th_self_mb=th_self_mb, wall_s=wall_s, ret=ret,
                n_features=n_feat, err=err_tail)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", default="test_data/CEN6_ver_220406.fasta",
                    help="input FASTA to chunk into parts (default: bundled CEN6)")
    ap.add_argument("--part-sizes", default="5,15,30",
                    help="comma-separated part sizes in MB (default: 5,15,30)")
    ap.add_argument("--rounds", default="round1_p40-3000,round2_p3001-10000,"
                    "round3_p10001-25000",
                    help="comma-separated round keys to run (default: all 3)")
    ap.add_argument("--threads", type=int, default=1,
                    help="TideHunter -t (default 1 = per-worker footprint)")
    ap.add_argument("--timeout", type=float, default=1800,
                    help="per-run timeout seconds (default 1800)")
    ap.add_argument("--out", default=None, help="TSV output path")
    ap.add_argument("--workdir", default=None,
                    help="scratch dir for parts + TH output (default: mkdtemp)")
    ap.add_argument("--keep", action="store_true", help="keep part/output files")
    args = ap.parse_args()

    check_tidehunter()
    if not os.path.exists(args.fasta):
        sys.exit(F"input FASTA not found: {args.fasta}")

    part_sizes_mb = [float(x) for x in args.part_sizes.split(",")]
    round_keys = [r for r in args.rounds.split(",") if r]
    for r in round_keys:
        if r not in ROUND_ARGS:
            sys.exit(F"unknown round '{r}'; choose from {list(ROUND_ARGS)}")

    work = args.workdir or tempfile.mkdtemp(prefix="th_mem_probe_")
    os.makedirs(work, exist_ok=True)
    out_tsv = args.out or os.path.join(work, "tidehunter_mem_probe.tsv")

    print(F"Input FASTA : {args.fasta} ({os.path.getsize(args.fasta)/1e6:.1f} MB)")
    print(F"Chunking to {CHUNK_SIZE} bp pieces (overlap {OVERLAP}) ...")
    chunked, _mt = tc.split_fasta_to_chunks(args.fasta, CHUNK_SIZE, OVERLAP)
    print(F"Chunked file: {os.path.getsize(chunked)/1e6:.1f} MB")
    print(F"Threads     : {args.threads}   Work dir: {work}\n")

    rows = []
    # build each part once, reuse across rounds
    for size_mb in part_sizes_mb:
        part_path = os.path.join(work, F"part_{size_mb:g}MB.fasta")
        in_bp, n_seqs = build_part(chunked, int(size_mb * 1e6), part_path)
        for rk in round_keys:
            th_args = ROUND_ARGS[rk]
            out_file = os.path.join(work, F"part_{size_mb:g}MB_{rk}.out")
            m = run_and_measure(part_path, th_args, args.threads, out_file,
                                args.timeout)
            peak_mb = m["ru_maxrss_kib"] / 1024.0    # authoritative
            vmhwm_mb = m["vmhwm_kib"] / 1024.0
            per_mbp = peak_mb / (in_bp / 1e6) if in_bp else 0.0
            rows.append(dict(round=rk, args=th_args, part_mb=size_mb,
                             input_bp=in_bp, n_seqs=n_seqs, threads=args.threads,
                             peak_mb=peak_mb, vmhwm_mb=vmhwm_mb,
                             th_self_mb=m["th_self_mb"], per_mbp=per_mbp,
                             wall_s=m["wall_s"], n_features=m["n_features"],
                             ret=m["ret"], err=m["err"]))
            status = "ok" if m["ret"] == 0 else F"FAIL(ret={m['ret']})"
            print(F"  {rk:22s} part={size_mb:>5g}MB in={in_bp/1e6:6.1f}Mbp "
                  F"peak={peak_mb:8.1f}MB ({per_mbp:6.1f} MB/Mbp) "
                  F"[vmhwm {vmhwm_mb:.0f} / th {m['th_self_mb']:.0f}] "
                  F"t={m['wall_s']:7.1f}s feats={m['n_features']:>8d} {status}"
                  + (F"  {m['err']}" if m["err"] else ""))
            if not args.keep and os.path.exists(out_file):
                os.remove(out_file)
        if not args.keep and os.path.exists(part_path):
            os.remove(part_path)

    if not args.keep and os.path.exists(chunked):
        os.remove(chunked)

    with open(out_tsv, "w") as f:
        f.write("round\targs\tpart_mb\tinput_bp\tn_seqs\tthreads\t"
                "peak_rss_mb\tvmhwm_mb\tth_selfreport_mb\tpeak_mb_per_mbp\t"
                "wall_s\tn_features\treturncode\n")
        for r in rows:
            f.write(F"{r['round']}\t{r['args']}\t{r['part_mb']:g}\t{r['input_bp']}\t"
                    F"{r['n_seqs']}\t{r['threads']}\t{r['peak_mb']:.1f}\t"
                    F"{r['vmhwm_mb']:.1f}\t{r['th_self_mb']:.1f}\t{r['per_mbp']:.2f}\t"
                    F"{r['wall_s']:.2f}\t{r['n_features']}\t{r['ret']}\n")

    # pool-sizing hint from the worst (highest-peak) row per part size
    print("\n=== pool-sizing hint (worst round per part size) ===")
    for size_mb in part_sizes_mb:
        sub = [r for r in rows if r["part_mb"] == size_mb and r["ret"] == 0]
        if not sub:
            continue
        worst = max(sub, key=lambda r: r["peak_mb"])
        print(F"  part={size_mb:g}MB : worst={worst['round']} "
              F"peak={worst['peak_mb']:.0f}MB  ->  at a 64 GB budget, "
              F"pool_size≈{int(64000 // max(1, worst['peak_mb']))}")
    print(F"\nTSV written: {out_tsv}")


if __name__ == "__main__":
    main()
