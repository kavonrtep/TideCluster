#!/usr/bin/env python3
"""End-to-end per-TRA consensus generation for a TideCluster run.

Builds one consensus monomer per Tandem Repeat Array (TRA) using two
complementary methods (array-MSA and TideHunter consensus), validates
each by self-BLAST against the source array, picks the better one per
TRA, and assigns a quality grade plus diagnostic flags.

The wrapper orchestrates four R scripts that live under
`tarean/consensus_prototype/`. Outputs land in
`<prefix>_per_tra_consensus/` (override with `--out-dir`).
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROTOTYPE_DIR = SCRIPT_DIR / "tarean" / "consensus_prototype"
RSCRIPT_PATHS = {
    "msa":      PROTOTYPE_DIR / "consensus_msa.R",
    "validate": PROTOTYPE_DIR / "tests" / "validate_drapa_blast.R",
    "th":       PROTOTYPE_DIR / "tests" / "validate_th_single.R",
    "selector": PROTOTYPE_DIR / "consensus_ensemble.R",
}


def banner(msg):
    line = "=" * len(msg)
    print(f"\n{line}\n{msg}\n{line}", flush=True)


def fail(msg, exit_code=1):
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(exit_code)


def check_executables(required):
    """Verify external commands are on PATH."""
    missing = [exe for exe in required if shutil.which(exe) is None]
    if missing:
        fail(f"required executables not found on PATH: {', '.join(missing)}")


def resolve_inputs(prefix):
    """Derive the four input paths from the TideCluster prefix."""
    prefix = Path(prefix).resolve()
    inputs = {
        "kite_tsv":         prefix.parent / f"{prefix.name}_kite" / "monomer_size_top3_estimats.csv",
        "tarean_dir":       prefix.parent / f"{prefix.name}_tarean" / "fasta",
        "tidehunter_gff3":  prefix.parent / f"{prefix.name}_tidehunter.gff3",
        "clustering_gff3":  prefix.parent / f"{prefix.name}_clustering.gff3",
    }
    missing = [(k, str(v)) for k, v in inputs.items() if not v.exists()]
    if missing:
        msg = "\n".join(f"  {k}: {p}" for k, p in missing)
        fail(f"missing required inputs derived from prefix '{prefix}':\n{msg}")
    return inputs


def run_step(label, cmd, log_path, force, sentinel, verbose):
    """Run one Rscript step, optionally skipping when its sentinel exists.

    `sentinel` is the file whose presence means "already done"; ignored
    when `force` is True.
    """
    if sentinel and not force and sentinel.exists():
        print(f"[skip] {label}: cached output present "
              f"({sentinel.relative_to(sentinel.anchor)})", flush=True)
        return
    banner(f"[run] {label}")
    if verbose:
        print("  cmd: " + " ".join(str(x) for x in cmd), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if proc.returncode != 0:
        sys.stderr.write(open(log_path).read())
        fail(f"{label} failed (exit {proc.returncode}); see {log_path}")
    print(f"[ok]  {label}  ({dt:.1f}s, log: {log_path.name})", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--prefix", required=True,
                        help="TideCluster run prefix (resolves all input paths)")
    parser.add_argument("-o", "--out-dir", default=None,
                        help="output directory (default: <prefix>_per_tra_consensus)")
    parser.add_argument("-c", "--cpu", type=int, default=4,
                        help="thread count for MSA + BLAST (default: 4)")
    parser.add_argument("--pident", type=int, default=85,
                        help="self-BLAST percent-identity cutoff (default: 85)")
    parser.add_argument("--evalue", default="1e-10",
                        help="self-BLAST e-value cutoff (default: 1e-10)")
    parser.add_argument("--cov-strict",         type=float, default=0.90,
                        help="grade A floor (default: 0.90)")
    parser.add_argument("--cov-qualified",      type=float, default=0.80,
                        help="grade B floor (default: 0.80)")
    parser.add_argument("--cov-low",            type=float, default=0.50,
                        help="grade C floor (default: 0.50)")
    parser.add_argument("--core-cutoff",        type=float, default=0.95,
                        help="boundary_overext flag trigger (default: 0.95)")
    parser.add_argument("--gap-cutoff-bp",      type=int,   default=1000,
                        help="internal_gap flag trigger in bp (default: 1000)")
    parser.add_argument("--pident-sd-cutoff",   type=float, default=1.5,
                        help="heterogeneous flag trigger (default: 1.5)")
    parser.add_argument("--pident-p05-cutoff",  type=float, default=92.0,
                        help="low_pident flag trigger (default: 92.0)")
    parser.add_argument("-f", "--force", action="store_true",
                        help="re-run cached intermediate steps")
    parser.add_argument("--no-keep-intermediate", action="store_true",
                        help="delete msa/ and th/ subdirectories after success "
                             "(saves disk; loses BLAST tsvs needed for re-grading)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="echo full Rscript command lines")
    args = parser.parse_args()

    # ---- pre-flight checks --------------------------------------------
    for key, path in RSCRIPT_PATHS.items():
        if not path.exists():
            fail(f"R script missing: {path}")
    check_executables(["Rscript", "mafft", "blastn", "makeblastdb"])
    inputs = resolve_inputs(args.prefix)

    out_dir = Path(args.out_dir) if args.out_dir else \
        Path(f"{args.prefix}_per_tra_consensus")
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    msa_dir      = out_dir / "msa"
    msa_blast    = msa_dir / "blast_validation"
    th_dir       = out_dir / "th"
    selected_dir = out_dir / "selected"
    log_dir      = out_dir / "logs"

    # Persist resolved arguments (matches TideCluster's args.json convention).
    args_record = {
        "prefix":              str(Path(args.prefix).resolve()),
        "inputs":              {k: str(v) for k, v in inputs.items()},
        "out_dir":              str(out_dir),
        "cpu":                  args.cpu,
        "pident":               args.pident,
        "evalue":               args.evalue,
        "cov_strict":           args.cov_strict,
        "cov_qualified":        args.cov_qualified,
        "cov_low":              args.cov_low,
        "core_cutoff":          args.core_cutoff,
        "gap_cutoff_bp":        args.gap_cutoff_bp,
        "pident_sd_cutoff":     args.pident_sd_cutoff,
        "pident_p05_cutoff":    args.pident_p05_cutoff,
        "force":                args.force,
        "keep_intermediate":    not args.no_keep_intermediate,
    }
    with open(out_dir / "args.json", "w") as fh:
        json.dump(args_record, fh, indent=2)

    t_start = time.time()

    # ---- step 1: array-MSA consensus per TRA --------------------------
    run_step(
        "1/4 array-MSA consensus",
        [
            "Rscript", str(RSCRIPT_PATHS["msa"]),
            "--mode",       "batch",
            "--kite-tsv",   str(inputs["kite_tsv"]),
            "--tarean-dir", str(inputs["tarean_dir"]),
            "--out-dir",    str(msa_dir),
            "--cpu",        str(args.cpu),
        ],
        log_path=log_dir / "01_msa.log",
        force=args.force,
        sentinel=msa_dir / "per_tra_consensus.fasta",
        verbose=args.verbose,
    )

    # ---- step 2: self-BLAST validation of MSA consensus ---------------
    # Builds the array DB used by step 3 too.
    run_step(
        "2/4 self-BLAST of array-MSA consensus",
        [
            "Rscript", str(RSCRIPT_PATHS["validate"]),
            "--consensus-fasta", str(msa_dir / "per_tra_consensus.fasta"),
            "--diagnostics-tsv", str(msa_dir / "per_tra_diagnostics.tsv"),
            "--tarean-dir",      str(inputs["tarean_dir"]),
            "--kite-tsv",        str(inputs["kite_tsv"]),
            "--out-dir",         str(msa_blast),
            "--pident",          str(args.pident),
            "--evalue",          str(args.evalue),
            "--cpu",             str(args.cpu),
        ],
        log_path=log_dir / "02_msa_blast.log",
        force=args.force,
        sentinel=msa_blast / "blast.tsv",
        verbose=args.verbose,
    )

    # ---- step 3: TideHunter consensus + self-BLAST --------------------
    # Re-uses the array DB built in step 2.
    blast_db = msa_blast / "all_arrays"
    run_step(
        "3/4 TideHunter consensus + self-BLAST",
        [
            "Rscript", str(RSCRIPT_PATHS["th"]),
            "--tidehunter-gff3", str(inputs["tidehunter_gff3"]),
            "--clustering-gff3", str(inputs["clustering_gff3"]),
            "--blast-db",        str(blast_db),
            "--baseline-tsv",    str(msa_blast / "per_tra_blast_metrics.tsv"),
            "--out-dir",         str(th_dir),
            "--pident",          str(args.pident),
            "--evalue",          str(args.evalue),
            "--cpu",             str(args.cpu),
        ],
        log_path=log_dir / "03_th_consensus.log",
        force=args.force,
        sentinel=th_dir / "blast.tsv",
        verbose=args.verbose,
    )

    # ---- step 4: selector + grading -----------------------------------
    # Always re-run; cheap and picks up any threshold change.
    run_step(
        "4/4 consensus selector + grading",
        [
            "Rscript", str(RSCRIPT_PATHS["selector"]),
            "--mb-fasta",          str(msa_dir / "per_tra_consensus.fasta"),
            "--mb-blast-tsv",      str(msa_blast / "blast.tsv"),
            "--th-fasta",          str(th_dir / "per_tra_th_consensus.fasta"),
            "--th-blast-tsv",      str(th_dir / "blast.tsv"),
            "--kite-tsv",          str(inputs["kite_tsv"]),
            "--out-dir",           str(selected_dir),
            "--cov-strict",        str(args.cov_strict),
            "--cov-qualified",     str(args.cov_qualified),
            "--cov-low",           str(args.cov_low),
            "--core-cutoff",       str(args.core_cutoff),
            "--gap-cutoff-bp",     str(args.gap_cutoff_bp),
            "--pident-sd-cutoff",  str(args.pident_sd_cutoff),
            "--pident-p05-cutoff", str(args.pident_p05_cutoff),
        ],
        log_path=log_dir / "04_selector.log",
        force=True,                              # always re-run
        sentinel=None,
        verbose=args.verbose,
    )

    # ---- promote final outputs to top of out-dir ----------------------
    for fname in ("per_tra_consensus.fasta",
                  "per_tra_metrics.tsv",
                  "summary.log"):
        src = selected_dir / fname
        if src.exists():
            shutil.copy2(src, out_dir / fname)
        else:
            fail(f"selector did not produce expected output: {src}")

    # ---- optional cleanup ---------------------------------------------
    if args.no_keep_intermediate:
        for sub in (msa_dir, th_dir):
            shutil.rmtree(sub, ignore_errors=True)
        print("[clean] removed intermediate msa/ and th/ directories",
              flush=True)

    dt_total = time.time() - t_start
    banner("done")
    print(f"output:    {out_dir}", flush=True)
    print(f"  per_tra_consensus.fasta  selected per-TRA consensus", flush=True)
    print(f"  per_tra_metrics.tsv      grade + flags + diagnostics", flush=True)
    print(f"  summary.log              totals and threshold values", flush=True)
    print(f"  logs/                    per-step Rscript stdout/stderr",
          flush=True)
    print(f"wall time: {dt_total:.1f}s", flush=True)


if __name__ == "__main__":
    main()
