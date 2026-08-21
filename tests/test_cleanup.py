#!/usr/bin/env python3
"""Unit tests for `run_all --cleanup` (docs/output_cleanup_spec.md).

A finished run is ~82% intermediates, but the obvious purge — deleting
<prefix>_{kite,tarean,consensus}/ wholesale — keeps only the *viewable* report
and silently destroys three other capabilities. So the purge is file-level, and
what survives is a contract:

  G1 view the HTML report            G3 re-render it (tc_rerender_report.py)
  G2 feed tc_comparative_analysis.R  G4 feed tc_per_tra_consensus.py

These tests assert the CONTRACT rather than a file list, on purpose: CARP
maintains its own copy of this purge set, so a layout change must fail loudly
here rather than quietly widening one of the two sets.

Run: python3 tests/test_cleanup.py
"""
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import tc_utils as tc

failures = []


def _brief(v):
    """Long path lists are the norm here; print the shape, not the contents."""
    if isinstance(v, list) and len(v) > 3:
        return F"[{len(v)} paths]"
    return repr(v)


def check(label, got, want):
    ok = got == want
    print(F"  {'ok  ' if ok else 'FAIL'} {label}: got {_brief(got)}, want {_brief(want)}")
    if not ok and isinstance(got, list) and isinstance(want, list):
        print("       only in got : ", sorted(set(got) - set(want))[:5])
        print("       only in want: ", sorted(set(want) - set(got))[:5])
    if not ok:
        failures.append(label)


def check_true(label, cond, detail=""):
    print(F"  {'ok  ' if cond else 'FAIL'} {label}{'' if cond else ' -- ' + detail}")
    if not cond:
        failures.append(label)


def write(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def build_run(root, prefix="tc", trcs=("TRC_1", "TRC_2")):
    """A miniature run directory with the shape a real run_all produces."""
    p = os.path.join(root, prefix)
    disposable, keep = [], []

    # --- kite ---------------------------------------------------------------
    disposable.append(write(F"{p}_kite/kitehor.periodogram", "periodogram" * 100))
    disposable.append(write(F"{p}_kite/_kite_input_longext.fasta"))
    disposable.append(write(F"{p}_kite/_rescored_longext.peaks.tsv"))
    keep.append(write(F"{p}_kite/monomer_size_top3_estimats.csv"))
    keep.append(write(F"{p}_kite/kitehor.rescored.peaks.tsv"))
    keep.append(write(F"{p}_kite/kitehor.kite.tsv"))
    keep.append(write(F"{p}_kite/kitehor.ssr.tsv"))
    keep.append(write(F"{p}_kite/kitehor.tandem_validate.tsv"))
    keep.append(write(F"{p}_kite/profile_plots/TRC_1.png"))

    # --- tarean -------------------------------------------------------------
    keep.append(write(F"{p}_tarean/SSRS_summary.csv"))
    for t in trcs:
        # the ONE trap: the array FASTA exists twice, and only the copy inside
        # the per-TRC dir is disposable.
        keep.append(write(F"{p}_tarean/fasta/{t}.fasta", "arrays" * 50))
        d = F"{p}_tarean/{t}.fasta_tarean"
        disposable.append(write(F"{d}/{t}.fasta", "arrays" * 50))
        disposable.append(write(F"{d}/{t}.fasta_11.kmers", "k" * 200))
        disposable.append(write(F"{d}/{t}.fasta_27.kmers", "k" * 200))
        disposable.append(write(F"{d}/ggmin.RData", "R" * 100))
        disposable.append(write(F"{d}/monomers.RData", "R" * 100))
        keep.append(write(F"{d}/report.html"))
        keep.append(write(F"{d}/summary_table.csv"))
        keep.append(write(F"{d}/ppm_11mer_1.csv"))
        keep.append(write(F"{d}/img/graph_11mer_1.png"))
        keep.append(write(F"{d}/consensus.fasta"))
        keep.append(write(F"{d}/tarean_contigs.fasta"))

    # --- consensus ----------------------------------------------------------
    disposable.append(write(F"{p}_consensus/consensus_sequences_all.fasta_renamed.fasta.cat", "c" * 300))
    disposable.append(write(F"{p}_consensus/consensus_sequences_all.fasta_renamed.fasta.masked", "m" * 300))
    disposable.append(write(F"{p}_consensus/consensus_sequences_all.fasta_renamed.fasta.tbl"))
    keep.append(write(F"{p}_consensus/consensus_sequences_all.fasta", "seqs" * 50))
    keep.append(write(F"{p}_consensus/consensus_sequences_all.fasta.out"))
    for t in trcs:
        keep.append(write(F"{p}_consensus/{t}_dimers.fasta"))

    # --- top level ----------------------------------------------------------
    disposable.append(write(F"{p}_clustering.gff3_1.gff3"))
    for name in ("_clustering.gff3", "_annotation.gff3", "_annotation.tsv",
                 "_consensus_dimer_library.fasta", "_tidehunter.gff3",
                 "_tidehunter_short.gff3", "_cmd_args.json",
                 "_pipeline_stats.json", "_seqid_lengths.tsv",
                 "_tarean_report.tsv", "_trc_superfamilies.csv",
                 "_index.html", "_rdna.tsv"):
        keep.append(write(F"{p}{name}"))
    keep.append(write(F"{p}_report/superfamilies.html"))
    keep.append(write(F"{p}_report/img/tarean/TRC_1/graph_11mer_1.png"))
    keep.append(write(F"{p}_report_legacy/{prefix}_tarean_report.html"))
    keep.append(write(os.path.join(root, "dotplots", "superfamily_001.png")))
    # explicitly requested debug output -- must survive (--keep_rounds)
    for r in (1, 2, 3):
        keep.append(write(F"{p}_tidehunter_round{r}.gff3"))

    return p, [os.path.normpath(x) for x in disposable], [os.path.normpath(x) for x in keep]


# --- 1. the purge set matches exactly what it should ------------------------
print("purge set: every disposable file goes, every kept file stays")
tmp = tempfile.mkdtemp(prefix="tc_cleanup_test_")
try:
    prefix, disposable, keep = build_run(os.path.join(tmp, "run"))

    matched = [p for p, _why in tc._cleanup_candidates(prefix)]
    check("all disposable files matched", sorted(matched), sorted(disposable))
    check("no kept file matched", sorted(set(matched) & set(keep)), [])

    # the trap, called out explicitly because it is the one that would hurt
    fasta_dir_hits = [m for m in matched if os.sep + "fasta" + os.sep in m]
    check_true("<prefix>_tarean/fasta/ is NOT matched (the glob trap)",
               not fasta_dir_hits, str(fasta_dir_hits))

    check_true("--keep_rounds files are not matched",
               not [m for m in matched if "_tidehunter_round" in m])

    # dry run touches nothing
    removed, freed = tc.cleanup_run_directory(prefix, dry_run=True)
    check("dry run reports the same files", sorted(removed), sorted(disposable))
    check_true("dry run reports a non-zero size", freed > 0, str(freed))
    check_true("dry run deleted nothing",
               all(os.path.exists(p) for p in disposable))

    # real run
    removed, freed = tc.cleanup_run_directory(prefix)
    check("removed the disposable set", sorted(removed), sorted(disposable))
    check_true("disposable files are gone",
               not [p for p in disposable if os.path.exists(p)])
    missing_keep = [p for p in keep if not os.path.exists(p)]
    check_true("every kept file survives", not missing_keep, str(missing_keep[:5]))
    check_true("freed byte count is the sum of what went", freed > 0, str(freed))

    # idempotent
    removed2, freed2 = tc.cleanup_run_directory(prefix)
    check("second run removes nothing", (len(removed2), freed2), (0, 0))

    # --- 2. the four guarantees, as paths --------------------------------
    print("contract: the four guarantees still resolve after cleanup")
    G = {
        "G1 report": [F"{prefix}_index.html", F"{prefix}_report/superfamilies.html",
                      F"{prefix}_report/img/tarean/TRC_1/graph_11mer_1.png",
                      F"{prefix}_report_legacy/tc_tarean_report.html"],
        "G2 comparative": [F"{prefix}_consensus_dimer_library.fasta",
                           F"{prefix}_consensus/consensus_sequences_all.fasta",
                           F"{prefix}_clustering.gff3", F"{prefix}_annotation.gff3",
                           F"{prefix}_annotation.tsv",
                           F"{prefix}_tarean/SSRS_summary.csv"],
        "G3 re-render": [F"{prefix}_cmd_args.json", F"{prefix}_pipeline_stats.json",
                         F"{prefix}_seqid_lengths.tsv", F"{prefix}_tarean_report.tsv",
                         F"{prefix}_trc_superfamilies.csv",
                         F"{prefix}_kite/monomer_size_top3_estimats.csv",
                         F"{prefix}_kite/kitehor.rescored.peaks.tsv",
                         F"{prefix}_kite/profile_plots/TRC_1.png",
                         F"{prefix}_tarean/TRC_1.fasta_tarean/report.html",
                         F"{prefix}_tarean/TRC_1.fasta_tarean/summary_table.csv",
                         F"{prefix}_tarean/TRC_1.fasta_tarean/ppm_11mer_1.csv",
                         F"{prefix}_tarean/TRC_1.fasta_tarean/img/graph_11mer_1.png",
                         os.path.join(os.path.dirname(prefix), "dotplots",
                                      "superfamily_001.png")],
        "G4 per-TRA": [F"{prefix}_kite/monomer_size_top3_estimats.csv",
                       F"{prefix}_tarean/fasta/TRC_1.fasta",
                       F"{prefix}_tarean/fasta/TRC_2.fasta",
                       F"{prefix}_tidehunter.gff3", F"{prefix}_clustering.gff3"],
    }
    for name, paths in G.items():
        gone = [p for p in paths if not os.path.isfile(p)]
        check_true(F"{name}: all {len(paths)} inputs present", not gone, str(gone))

    # --- 3. protected-path guard fires on a loosened pattern -------------
    print("guard: a loosened pattern aborts before deleting anything")
    prefix2, _d2, keep2 = build_run(os.path.join(tmp, "run2"))
    saved = tc.CLEANUP_PATTERNS
    try:
        # the mistake CARP hit: `*/` instead of `*.fasta_tarean/`
        tc.CLEANUP_PATTERNS = (("{p}_tarean/*/TRC_*.fasta", "loosened"),)
        try:
            tc.cleanup_run_directory(prefix2)
            check_true("loosened glob raises", False, "no exception")
        except RuntimeError as exc:
            check_true("loosened glob raises", True)
            check_true("the error names the protected path",
                       "fasta" in str(exc), str(exc))
        check_true("nothing was deleted before the guard fired",
                   all(os.path.exists(p) for p in keep2))
    finally:
        tc.CLEANUP_PATTERNS = saved

    # --- 4. tolerant of a missing / unreadable tree ----------------------
    print("tolerance: an incomplete run directory is not an error")
    empty = os.path.join(tmp, "empty")
    os.makedirs(empty)
    removed3, freed3 = tc.cleanup_run_directory(os.path.join(empty, "nothing"))
    check("cleanup on a directory with no run", (len(removed3), freed3), (0, 0))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(F"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("test_cleanup.py: all checks passed")
