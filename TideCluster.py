#!/usr/bin/env python
"""Wrapper of TideHunter
Input is a fasta file of DNA sequences
Output is a table of detected tandem repeats by tidehunter
this wrapper split DNA sequence to chunks, run tidehunter on each chunk and merge results
coordinates in the table must be recalculated to the original sequence

"""
import glob
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool
import tc_utils as tc
import argparse
from version import __version__

# minimal python version is 3.6
assert sys.version_info >= (3, 6), "Python 3.6 or newer is required"


# --max_memory is offered by every subcommand that sizes a worker pool
# (tidehunter, tarean, run_all); the resolution chain behind it lives in
# tc_utils.memory_budget_mb.
MAX_MEMORY_HELP = (
    "Memory limit for this run, in GB. Used to size the TideHunter worker pool "
    "and cap TAREAN threads. Set this on a cluster or in a container: without "
    "it the budget falls back to the scheduler environment (PBS_RESC_MEM, "
    "SLURM_MEM_PER_NODE, ...), then the cgroup limit, then /proc/meminfo "
    "MemAvailable -- which reports the whole node's memory, not the job's limit."
)


def _cleanup_outputs(prefix):
    """Delete the run's intermediates and report what went.

    Called only at the very end of a SUCCESSFUL run_all -- never from a
    `finally:`, never for an individual subcommand. What survives is a contract
    (report viewing, comparative-analysis input, report re-render, per-TRA
    consensus); see docs/output_cleanup_spec.md and tc_utils.CLEANUP_PATTERNS.
    """
    removed, freed = tc.cleanup_run_directory(prefix)
    print(F"cleanup: removed {len(removed)} intermediate file(s), freed "
          F"{freed / 1e6:.1f} MB. The run directory is no longer byte-complete; "
          F"report, re-render, comparative-analysis and per-TRA consensus inputs "
          F"are kept.")
    # Record it in the run's own stats so a pruned directory says so months
    # later. The report is built before cleanup runs, so this shows up only on
    # a subsequent re-render -- the JSON is the durable record either way.
    stats_path = F"{prefix}_pipeline_stats.json"
    if os.path.exists(stats_path):
        try:
            with open(stats_path) as f:
                stats = json.load(f)
            stats["cleanup_files_removed"] = len(removed)
            stats["cleanup_bytes_freed"] = freed
            with open(stats_path, "w") as f:
                json.dump(stats, f, indent=2)
        except (OSError, ValueError) as e:
            print(F"WARNING: could not record cleanup in {stats_path}: {e}",
                  file=sys.stderr)


def _run_required(cmd, what):
    """Run a pipeline step whose failure invalidates the run's outputs.

    ``tc.run_cmd`` returns ``[cmd, 'error']`` and prints the child's stderr, but
    every call site used to discard that status, so a dead step left the
    pipeline running to a report that quietly denied its own missing results --
    a killed ``compare_trc_by_blast.R`` produced a report stating "no TRC
    superfamilies were identified" for a genome with 75 of them, and ``run_all``
    still exited 0 (issue #7). Aborting is the safe default: a non-zero exit is
    recoverable, a confidently wrong report is not.
    """
    if tc.run_cmd(cmd)[1] == 'error':
        raise RuntimeError(
            F"{what} failed; its outputs are missing or incomplete. See the "
            F"error above. Command: {cmd}")


def _run_optional(cmd, what):
    """Run a step whose failure degrades the report but does not invalidate it.

    Used for purely presentational artefacts (e.g. profile plots): say so
    loudly, then carry on. The counterpart of :func:`_run_required`; the point
    of having both is that every ``run_cmd`` result is now looked at.
    """
    if tc.run_cmd(cmd)[1] == 'error':
        print(F"WARNING: {what} failed; the run continues without it "
              F"(report will be missing these outputs).", file=sys.stderr)
        return False
    return True


def _tarean_max_threads(cpu, max_memory=None):
    """Cap the per-job TAREAN thread count by a memory budget.

    A TAREAN job's internal ``mclapply`` forks scale with the TRC size, so on a
    memory-constrained host giving a large TRC many threads can exhaust RAM.
    Cap threads at ``budget / per_thread``, where the budget and its origin come
    from :func:`tc_utils.memory_budget_mb` (``--max_memory``, then the scheduler
    environment, then the cgroup limit, then host ``MemAvailable``); if nothing
    is readable, fall back to ``cpu`` (no cap). ``per_thread`` is a deliberately
    conservative constant -- this is a safety cap against runaway forks, not a
    precise controller.

    Returns ``(max_threads, budget_mb, source)`` so the caller can report which
    source won: a cap computed from a host-wide ``MemAvailable`` under a batch
    scheduler is not a cap at all (issue #6).
    """
    per_thread_mb = 4000.0
    budget_mb, source = tc.memory_budget_mb(max_memory)
    tc.warn_if_host_memory_budget(source)
    if budget_mb is None:
        return cpu, None, source
    return max(1, min(cpu, int(budget_mb / per_thread_mb))), budget_mb, source


def tarean(prefix, gff, fasta=None, cpu=4, min_total_length=50000, args=None,
           version=__version__, max_memory=None):
    """
    Run tarean on genomic sequences specified in gff3 file
    from gff record extract sequence end analyse it with tarean algorithm
    :param prefix: prefix for output files, if gff is not specied it is used to find gff3
    file
    :param gff: gff3 file with tidehunter results
    :param fasta: reference fasta
    :param cpu: number of cpu cores to use
    :param min_total_length: minimal total length of sequences to run tarean
    :param max_memory: memory limit in GB (--max_memory); caps per-job threads
    :return:
    """
    script_path = os.path.dirname(os.path.realpath(__file__))
    if gff is None:
        # use preferentially gff3 file with annotation
        if os.path.exists(prefix + "_annotation.gff3"):
            gff = prefix + "_annotation.gff3"
        elif os.path.exists(prefix + "_clustering.gff3"):
            gff = prefix + "_clustering.gff3"
        else:
            print(F"gff3 file annotation not found in {prefix}")
            return

    # directory for tarean results
    tarean_dir = prefix + "_tarean"
    # directory for fasta files extracted from gff3, it could be deleted after tarean run
    tarean_dir_fasta = tarean_dir + "/fasta"
    # create directory for tarean results
    if not os.path.exists(tarean_dir):
        # create directory for tarean results
        os.mkdir(tarean_dir)
    # create directory for fasta files
    if not os.path.exists(tarean_dir_fasta):
        os.mkdir(tarean_dir_fasta)

    fasta_dict = tc.extract_sequences_from_gff3(gff, fasta, tarean_dir_fasta)
    # get per-contig lengths and the total (reused below for the summary
    # block and persisted as a TSV side-car for downstream tools e.g.
    # tc_rerender_report.py genome-distribution visualisation).
    seqid_lengths = tc.read_fasta_sequence_size(fasta)
    input_fasta_length = sum(seqid_lengths.values())
    with open(F"{prefix}_seqid_lengths.tsv", "w") as _f:
        _f.write("seqid\tlength\n")
        for _sid, _len in sorted(seqid_lengths.items(),
                                 key=lambda kv: (-kv[1], kv[0])):
            _f.write(F"{_sid}\t{int(_len)}\n")
    l_debug = 0
    cmd_list = []

    ssr = {}
    with open(gff, "r") as f:
        for i in f:
            if i.startswith("#"):
                continue
            gff_record = tc.Gff3Feature(i)
            if "ssr" in gff_record.attributes:
                ssr[gff_record.attributes_dict['Name']] = gff_record.attributes_dict[
                    "ssr"]

    print("preparaing sequences for TAREAN")
    omitted_clusters = []
    trc_total_length = int(0)
    trc_total_length_omitted = 0
    for k, v in fasta_dict.items():
        l_debug += 1
        v_basename = os.path.basename(v)
        with open(v, "r") as f:
            seqs = tc.read_single_fasta_to_dictionary(f)
        total_length = sum([len(i) for i in seqs.values()]) / 2  # it is dimer!!
        trc_total_length += total_length
        if total_length < min_total_length:
            trc_total_length_omitted += total_length
            print(
                    F"total length of sequences in {k} is less than {min_total_length} "
                    F"nt, "
                    F"skipping"
                    )
            omitted_clusters.append((k, total_length, len(seqs)))
            continue
        # does seqs contain only one sequence?
        if len(seqs) > 1:
            # more than one sequence in fasta file, set orientation
            seqs2 = tc.group_sequences_by_orientation(seqs, k=8)
            for n in seqs:
                if n in seqs2['reverse']:
                    seqs[n] = tc.reverse_complement(seqs[n])
        # save oriented sequences to file
        tc.save_fasta_dict_to_file(seqs, v)
        # get current script path
        # do not run if ssrs
        if k in ssr:
            print(F"{k} is SSR, skipping")
            continue
        tarean_out = F"{tarean_dir}/{v_basename}_tarean"
        # Store (total array length, input fasta, output dir); the per-job -n
        # thread count and dispatch order are decided below so large TRCs can be
        # given more cores. tarean.R output is thread-count-independent, so this
        # changes only wall time, not results.
        cmd_list.append((total_length, v, tarean_out))
    # run cmd tarean in parallel using multiprocessing module
    if len(omitted_clusters) > 0:
        with open(F"{tarean_dir}/omitted_clusters.txt", "w") as f:
            f.write("cluster_id\ttotal_length\tnumber_of_arrays\n")
            # sort by total length
            omitted_clusters.sort(key=lambda x: x[1], reverse=True)
            for i in omitted_clusters:
                f.write(F"{i[0]}\t{i[1]}\t{i[2]}\n")

    # RUN kitehor (Rust k-mer periodicity + rescore + SSR scan).
    #
    # Pipeline (kitehor >=0.12.0): kite-periodicity emits the kite top-3
    # TSV + the long-format peaks TSV + the H[d]/bg periodogram bundle,
    # then `rescore` and `ssr-scan` run **in parallel** against that
    # peaks TSV. The slow `analyze` cascade (rule-classify + tandem-
    # validate + summary-merge) is dropped — TideCluster derives the
    # per-array founder / strongest / subrepeat-candidate structure
    # itself from rescore's identity_med + scan_occupancy_frac etc.
    # See docs/kitehor_integration_plan.md for the full mapping.
    kite_dir = F"{prefix}_kite"
    os.makedirs(kite_dir, exist_ok=True)
    multi_fa = F"{kite_dir}/_kite_input.fasta"
    n_kite_records = tc.build_kite_multifasta(
        F"{prefix}_tarean/fasta", multi_fa)
    if n_kite_records == 0:
        print("No arrays passed to kitehor; skipping KITE step.")
    else:
        rescore_max_period = int(getattr(args, "kite_rescore_max_period", 10000))
        rescore_top_n      = int(getattr(args, "kite_rescore_top_n",      20))
        print(F"Running kitehor kite-periodicity on {n_kite_records} array(s).")
        _run_required(F"kitehor kite-periodicity {multi_fa}"
                      F" --out {kite_dir}/kitehor.kite.tsv"
                      F" --out-peaks {kite_dir}/kitehor.kite.peaks.tsv"
                      F" --periodogram {kite_dir}/kitehor.periodogram"
                      F" --threads {cpu}", "kitehor kite-periodicity")
        # rule-classify is cheap and produces the per-array verdicts that
        # `tandem-validate` (kitehor >= 0.13.0's unified subrepeat detector,
        # spec v5) consumes. Run it first, then fan rescore || ssr-scan ||
        # tandem-validate out concurrently. rayon inside each binary still
        # uses --threads; rescore is the O(period²) bottleneck so it keeps
        # the lion's share of the CPU budget, the other two are quick.
        _run_required(F"kitehor rule-classify {kite_dir}/kitehor.kite.peaks.tsv"
                      F" --out {kite_dir}/kitehor", "kitehor rule-classify")
        rescore_threads = max(1, cpu - 2)
        side_threads    = max(1, cpu // 2)
        rescore_cmd = (
            F"kitehor rescore"
            F" --peaks {kite_dir}/kitehor.kite.peaks.tsv"
            F" --out {kite_dir}/kitehor.rescored"
            F" --max-period {rescore_max_period}"
            F" --top-n {rescore_top_n}"
            F" --threads {rescore_threads} {multi_fa}")
        ssr_cmd = (
            F"kitehor ssr-scan"
            F" --kite-peaks {kite_dir}/kitehor.kite.peaks.tsv"
            F" --out {kite_dir}/kitehor {multi_fa}")
        # tandem-validate: founder/host monomer is inferred from the
        # verdicts; per array it reports the dominant nested-TR candidate
        # with density (occupancy) / spatial_contrast / phase_contrast and
        # a decision_hint (localized_subrepeat | confirms_host | ...).
        # TideCluster gates these against its OWN founder downstream.
        tv_cmd = (
            F"kitehor tandem-validate"
            F" --verdicts {kite_dir}/kitehor.verdicts.tsv"
            F" --peaks {kite_dir}/kitehor.kite.peaks.tsv"
            F" --out {kite_dir}/kitehor"
            F" --threads {side_threads} {multi_fa}")
        print(F"Running kitehor rescore (max-period {rescore_max_period},"
              F" top-n {rescore_top_n}) || ssr-scan || tandem-validate"
              F" in parallel.")
        rescore_proc = subprocess.Popen(rescore_cmd, shell=True)
        ssr_proc     = subprocess.Popen(ssr_cmd,     shell=True)
        tv_proc      = subprocess.Popen(tv_cmd,      shell=True)
        rc1 = rescore_proc.wait()
        rc2 = ssr_proc.wait()
        rc3 = tv_proc.wait()
        if rc1 != 0:
            raise RuntimeError(F"kitehor rescore failed (exit {rc1})")
        if rc2 != 0:
            raise RuntimeError(F"kitehor ssr-scan failed (exit {rc2})")
        if rc3 != 0:
            raise RuntimeError(F"kitehor tandem-validate failed (exit {rc3})")
        # Selective long-period extension: arrays whose dominant monomer is
        # above the rescore cap get NA identity and fall back to a spurious
        # short peak. Re-rescore only those at the higher cap (cheap — rescore
        # is O(period²), confined to the few flagged arrays). See case 4.
        ext_max_period = int(getattr(args, "kite_rescore_max_period_ext", 25000))
        if ext_max_period > rescore_max_period:
            flagged = tc.extend_long_period_rescore(
                kite_dir, multi_fa, rescore_max_period, ext_max_period,
                rescore_top_n, rescore_threads)
            if flagged:
                print(F"Re-rescored {len(flagged)} array(s) with a dominant "
                      F"period > {rescore_max_period} bp at max-period "
                      F"{ext_max_period}.")
        # SSR families are classified at clustering (repeat_type=SSR + ssr motif
        # in the clustering GFF3); their founder is the fundamental motif length
        # for all arrays, set in build_monomer_size_csv.
        trc_repeat_type = tc.parse_trc_ssr_motif_len(F"{prefix}_clustering.gff3")
        tc.build_monomer_size_csv(
            kite_tsv=F"{kite_dir}/kitehor.kite.tsv",
            ssr_tsv=F"{kite_dir}/kitehor.ssr.tsv",
            rescored_peaks_tsv=F"{kite_dir}/kitehor.rescored.peaks.tsv",
            tandem_validate_tsv=F"{kite_dir}/kitehor.tandem_validate.tsv",
            out_csv=F"{kite_dir}/monomer_size_top3_estimats.csv",
            trc_repeat_type=trc_repeat_type)
        print("Rendering per-TRC profile heatmaps.")
        _run_optional(F"{script_path}/tarean/kite_heatmaps.R"
                      F" --periodogram {kite_dir}/kitehor.periodogram"
                      F" --top3-csv {kite_dir}/monomer_size_top3_estimats.csv"
                      F" --out-dir {kite_dir}/profile_plots",
                      "KITE profile heatmaps (kite_heatmaps.R)")
    if os.path.exists(multi_fa):
        os.remove(multi_fa)


    # copy index.html as prefix_index.html
    html_src = script_path + "/tarean/index.html"
    html_dst = F"{prefix}_index.html"
    # Format string to report analysis settings (from args and version)
    saved = load_args_from_file(prefix)
    for key, val in saved.items():
        if not hasattr(args, key):
            setattr(args, key, val)

    # Use original_fasta if available (preserves user-provided path), otherwise use args.fasta
    input_fasta = getattr(args, 'original_fasta', args.fasta)

    # Determine TideHunter mode and format info
    tidehunter_mode = getattr(args, 'long', False)
    if tidehunter_mode:
        th_mode_str = "Three-round (--long)"
        th_args_str = "Round 1: p=40-3000, Round 2: p=3001-10000, Round 3: p=10001-25000"
    else:
        th_mode_str = "Standard"
        th_args_str = getattr(args, 'tidehunter_arguments', 'n/a')

    # The info block summarises run_all-level settings. When tarean is invoked
    # standalone, run_all-only args (tidehunter_arguments, min_length, no_dust,
    # library) are normally back-filled from <prefix>_cmd_args.json above; fall
    # back to 'n/a' so a standalone run without that side-car still prints the
    # block instead of raising AttributeError.
    # Record the run's memory assumption alongside the CPU count: an unset
    # --max_memory means the budget was inferred (scheduler env / cgroup /
    # host MemAvailable), which is exactly the case that can silently be the
    # node's memory rather than the job's (issue #6).
    _mm = getattr(args, 'max_memory', None)
    mem_limit_str = (F"{_mm:g} GB (--max_memory)" if _mm else
                     "auto (scheduler env / cgroup / MemAvailable)")

    settings = (F"Input file                 : {input_fasta}\n"
                F"Prefix                     : {args.prefix}\n"
                F"Minimum TRC total length   : {args.min_total_length}\n"
                F"Minimum array length       : {getattr(args, 'min_length', 'n/a')}\n"
                F"Dust filter                : {'no' if getattr(args, 'no_dust', False) else 'yes'}\n"
                F"TideHunter mode            : {th_mode_str}\n"
                F"TideHunter arguments       : {th_args_str}\n"
                F"CPU                        : {args.cpu}\n"
                F"Memory limit               : {mem_limit_str}\n"
                F"Library                    : {getattr(args, 'library', 'n/a')}\n"
                F"TideCluster version        : {version}\n")

    # create summar with number of clusters, number of omitted clusters, number of SSRs
    # total length all clusters, total length of reference
    # Count number of features in clustering.gff3 file
    gff_feature_count = sum(1 for line in open(gff) if not line.startswith("#"))

    summary = (
        F"Number of TRCs                  : {l_debug}\n"
        F"Number of TRCs above threshold  : {l_debug - len(omitted_clusters)}\n"
        F"Number of SSRs in TRCs          : {len(ssr)}\n"
        F"Total length of TRCs            : {int(trc_total_length)} nt\n"
        F"Number of TRAs                  : {gff_feature_count}\n"
        F"Input sequence length           : {input_fasta_length} nt\n"
    )

    # Persist the same numbers as structured JSON so downstream tools
    # (tc_rerender_report.py, external parsers) don't have to scrape the
    # rendered index.html for them.
    stats_json = {
        "n_trcs_total":           l_debug,
        "n_trcs_above_threshold": l_debug - len(omitted_clusters),
        "n_ssrs":                 len(ssr),
        "total_tr_length":        int(trc_total_length),
        "n_tras":                 gff_feature_count,
        "input_sequence_length":  input_fasta_length,
        "tidecluster_version":    version,
    }
    with open(F"{prefix}_pipeline_stats.json", "w") as f:
        json.dump(stats_json, f, indent=2)

    # replace all PREFIX_PLACEHOLDER with prefix value and save to new file
    # replace all SETTINGS_PLACEHOLDER with settings value and save to new file
    with open(html_src, "r") as f, open(html_dst, "w") as f2:
        for line in f:
            prefix_basename = os.path.basename(prefix)
            new_line = line.replace("PREFIX_PLACEHOLDER", prefix_basename)
            new_line = new_line.replace("SETTINGS_PLACEHOLDER", settings)
            new_line = new_line.replace("SUMMARY_PLACEHOLDER", summary)
            f2.write(new_line)

    if len(cmd_list) == 0:
        print("No remaining sequences for TAREAN analysis; exiting")
        with open(F"{prefix}_tarean_report.html", "w") as f:
            f.write("No TRC passed the minimum total length threshold for TAREAN "
                    "analysis. ")
            # create empty TRC clustering report
            html_file = F"{prefix}_trc_superfamilies.html"
            # print to html file
            print("TRC similarity clustering not performed due to lack of TAREAN "
                  "consensus sequences.", file=open(html_file, "w"))

        # This no-TAREAN path never runs compare_trc_by_blast.R, so emit the
        # canonical empty superfamily CSV + manifest here too (mirrors that
        # script's make_empty_outputs) -- keeps the artefact contract uniform
        # regardless of which path produced the run.
        _write_empty_superfamily_outputs(prefix)

        _maybe_identify_rdna(prefix, fasta, args, cpu)
        _move_v1_to_legacy(prefix)
        _build_report_v2(prefix)
        return

    print("running TAREAN")

    def _tarean_cmd(input_fasta, out_dir, threads):
        return (F"{script_path}/tarean/tarean.R -i {input_fasta} -s 0 "
                F"-n {threads} -o {out_dir}")

    # TAREAN jobs are independent (each writes its own output dir) and tarean.R
    # output is thread-count-independent (verified: -n1 and -nK give identical
    # consensus / summary), so per-job thread counts and dispatch order can be
    # chosen freely to minimise wall time without changing results. TAREAN cost
    # scales steeply with array size and TRC lengths span ~1000x, so one large
    # TRC can otherwise dominate the whole stage on a single core. Sort
    # longest-first (LPT); (-len, in, out) is a deterministic total order.
    jobs = sorted(cmd_list, key=lambda x: (-x[0], x[1], x[2]))
    total_jobs = len(jobs)
    max_threads, budget_mb, budget_src = _tarean_max_threads(cpu, max_memory)
    if budget_mb is None:
        print(F"TAREAN: no memory budget ({budget_src}); threads capped by "
              F"-c {cpu} only")
    else:
        print(F"TAREAN: memory budget ~{budget_mb:.0f} MB from {budget_src} -> "
              F"at most -n {max_threads} per job (cap {cpu})")
    completed = [0]
    failed_jobs = []

    def _progress(_res):
        completed[0] += 1
        # run_cmd returns [cmd, 'ok'|'error']; a per-TRC failure costs that
        # TRC its consensus but leaves the rest of the run meaningful, so it
        # is collected and reported rather than aborting mid-stage. What it
        # must not do is pass unnoticed (issue #7).
        if isinstance(_res, (list, tuple)) and len(_res) == 2 and _res[1] == 'error':
            failed_jobs.append(_res[0])
        print(F"completed {completed[0]} of {total_jobs}")

    W = sum(length for length, _i, _o in jobs)
    if total_jobs <= cpu:
        # Fewer jobs than cores: a -n1 pool would leave cores idle. Give every
        # job an equal share of the cores and run them all at once.
        threads = max(1, min(max_threads, cpu // total_jobs))
        print(F"TAREAN: {total_jobs} job(s) <= {cpu} core(s) -> all concurrent, "
              F"-n {threads} each")
        cmds = [_tarean_cmd(i, o, threads) for _l, i, o in jobs]
        with Pool(min(cpu, total_jobs)) as p:
            for _res in p.imap(tc.run_cmd, cmds):
                _progress(_res)
    else:
        # Saturated pool. The bulk pools efficiently at -n1, but an extreme
        # outlier TRC would run ~alone at the tail on one core (the reported
        # 55h straggler). Pull out the leading jobs that are each larger than
        # the whole remaining pool's per-core load (len > sum(smaller)/cpu):
        # everything else finishes before such a job, so it is on the critical
        # path and worth multi-threading. Run the -n1 pool first (all cores on
        # the bulk), then the stragglers multi-threaded (all cores each) -- both
        # phases keep every core busy.
        k = 0
        suffix = W
        for length, _i, _o in jobs:
            rest_below = suffix - length
            if rest_below > 0 and length > rest_below / cpu:
                k += 1
                suffix -= length
            else:
                break
        stragglers, bulk = jobs[:k], jobs[k:]
        if stragglers:
            print(F"TAREAN: {total_jobs} jobs -> {len(bulk)} pooled at -n1, then "
                  F"{len(stragglers)} straggler(s) at -n {max_threads}")
        else:
            print(F"TAREAN: {total_jobs} jobs pooled at -n1 (no straggler outlier)")
        bulk_cmds = [_tarean_cmd(i, o, 1) for _l, i, o in bulk]
        with Pool(cpu) as p:
            for _res in p.imap(tc.run_cmd, bulk_cmds):
                _progress(_res)
        for _l, i, o in stragglers:
            _progress(tc.run_cmd(_tarean_cmd(i, o, max_threads)))

    if failed_jobs:
        print(F"WARNING: {len(failed_jobs)} of {total_jobs} TAREAN job(s) "
              F"failed; those TRCs have no consensus and are absent from the "
              F"report. Failed commands:", file=sys.stderr)
        for c in failed_jobs:
            print(F"  {c}", file=sys.stderr)
        if len(failed_jobs) == total_jobs:
            raise RuntimeError(
                F"all {total_jobs} TAREAN jobs failed; nothing to report on")
    print("TAREAN finished")
    # get SSR info for tarean report from gff3 file

    # export SSR info to csv file
    with open(F"{tarean_dir}/SSRS_summary.csv", "w") as f:
        for k, v in ssr.items():
            f.write(F"{k}\t{v}\n")

    # final tarean report:
    cmd = (F"{script_path}/tarean/tarean_report.R -i {tarean_dir} -o"
           F" {prefix}_tarean_report -g {gff}")
    print("Making final tarean report.")
    # tarean_report.R builds the consensus dimer library that the superfamily
    # step below consumes, so its failure cascades.
    _run_required(cmd, "tarean_report.R")

    # Compare TRC by blast   - it require consensus dimers library generated by tarean_report.R.
    # --consensus_dir + --annotation_tsv enable the below-TAREAN-threshold
    # fallback: TRCs whose array total length fell below `min_total_length`
    # (default 50 kb) have no TAREAN consensus, so are absent from the dimer
    # library. The fallback BLASTs their raw TideHunter dimer consensus
    # (preserved by the clustering step) against the dimer-library DB and
    # attaches qualifying small TRCs to an existing superfamily, or promotes
    # a (small, previously-singleton-big) pair to a new SF. See the R
    # script for the score gate and annotation-consistency safety check.
    consensus_dir_arg = F"{prefix}_consensus"
    annotation_tsv_arg = F"{prefix}_annotation.tsv"
    superfamily_score = float(getattr(args, "superfamily_score", 20))
    print('running Compare TRC by blast')
    cmd = (
        F"{script_path}/tarean/compare_trc_by_blast.R -i {prefix}_consensus_dimer_library.fasta"
        F" -p {prefix} -t {cpu}"
        F" --consensus_dir {consensus_dir_arg}"
        F" --annotation_tsv {annotation_tsv_arg}"
        F" --score_threshold {superfamily_score}")
    _run_required(cmd, "compare_trc_by_blast.R (superfamily analysis)")

    _maybe_identify_rdna(prefix, fasta, args, cpu)
    _move_v1_to_legacy(prefix)
    _build_report_v2(prefix)


def _write_empty_superfamily_outputs(prefix):
    """Write the canonical empty superfamily CSV + manifest for a run that
    produced no superfamilies.

    Mirrors tarean/compare_trc_by_blast.R:make_empty_outputs so the artefact
    contract is uniform no matter which path ran: the CSV is always present
    under the canonical <prefix>_trc_superfamilies.csv name with the same
    header and zero data rows, and a small manifest declares the outputs.
    The caller writes the accompanying HTML stub."""
    base = os.path.basename(prefix)
    # Header-only CSV, same quoted header + schema as write.csv() in the R side.
    with open(F"{prefix}_trc_superfamilies.csv", "w") as f:
        f.write('"Superfamily","TRC","fallback"\n')
    manifest = {
        "producer": "TideCluster/TideCluster.py",
        "schema_version": 1,
        "superfamilies_found": False,
        "n_superfamilies": 0,
        "outputs": {
            "csv": F"{base}_trc_superfamilies.csv",
            "html": F"{base}_trc_superfamilies.html",
        },
        "csv_columns": ["Superfamily", "TRC", "fallback"],
    }
    with open(F"{prefix}_trc_superfamilies.manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)


def _move_v1_to_legacy(prefix):
    """Move the four top-level v1 HTML reports into <prefix>_report_legacy/
    so the root of the output directory is left clean for the v2 landing
    page (`<prefix>_index.html`, written by _build_report_v2).

    Embedded `src=`/`href=` references that point into the v1 data
    directories (<prefix>_tarean/, <prefix>_kite/, dotplots/) get a
    leading '../' prepended so they still resolve from the legacy dir.

    Idempotent: if a file is already missing (e.g. pipeline ran only
    partial steps) it is simply skipped. Data directories are not
    moved — they contain per-TRC v1 HTML that references its own
    siblings with local paths."""
    prefix_dir  = os.path.dirname(prefix) or "."
    prefix_name = os.path.basename(prefix)
    legacy_dir  = F"{prefix}_report_legacy"
    os.makedirs(legacy_dir, exist_ok=True)

    def _move_and_rewrite(src, dst, path_prefixes):
        if not os.path.exists(src):
            return
        with open(src) as f:
            html = f.read()
        for p in path_prefixes:
            html = html.replace(F'src="{p}',  F'src="../{p}')
            html = html.replace(F'href="{p}', F'href="../{p}')
        with open(dst, "w") as f:
            f.write(html)
        os.remove(src)

    # Index: sibling links (e.g. drapa_tarean_report.html) resolve inside
    # the legacy dir once tarean_report.html is moved there too, so no
    # ref rewriting is needed for the index itself.
    _move_and_rewrite(
        F"{prefix}_index.html",
        F"{legacy_dir}/{prefix_name}_index.html",
        [])
    _move_and_rewrite(
        F"{prefix}_tarean_report.html",
        F"{legacy_dir}/{prefix_name}_tarean_report.html",
        [F"{prefix_name}_tarean/"])
    _move_and_rewrite(
        F"{prefix}_kite_report.html",
        F"{legacy_dir}/{prefix_name}_kite_report.html",
        [F"{prefix_name}_kite/"])
    _move_and_rewrite(
        F"{prefix}_trc_superfamilies.html",
        F"{legacy_dir}/{prefix_name}_trc_superfamilies.html",
        ["dotplots/"])
    print(F"Legacy v1 reports moved to {legacy_dir}/")


def _build_report_v2(prefix):
    """Generate the modern v2 HTML report alongside the legacy output.

    Runs at the end of tarean() (both the normal and the early-return
    paths) so every pipeline run ships a v2 report. Wrapped in
    try/except: a rerender failure never fails the pipeline, since
    the legacy v1 HTML is already written by the time we get here."""
    try:
        import tc_rerender_report
        prefix_dir = os.path.dirname(os.path.abspath(prefix)) or "."
        prefix_name = os.path.basename(prefix)
        print("Building report v2")
        tc_rerender_report.build_report(prefix_dir, prefix=prefix_name, quiet=True)
        print(F"Report v2 written: {prefix_dir}/{prefix_name}_index.html"
              F" + {prefix_dir}/{prefix_name}_report/")
    except Exception as e:
        print(F"WARNING: report v2 generation failed: {e}", file=sys.stderr)


def _default_rdna_library():
    """Path to the bundled rDNA reference library (next to this script)."""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        "data", "rdna_library.fasta")


def _maybe_identify_rdna(prefix, fasta, args, cpu):
    """Run rDNA (45S/5S) identification before the report is built.

    Default-on; disabled by --no_rdna. Uses the bundled rDNA library unless
    --rdna_library overrides it. Wrapped in try/except so a failure never
    fails the pipeline (the labels are an enrichment, like the v2 report)."""
    if getattr(args, "no_rdna", False):
        return
    rdna_library = getattr(args, "rdna_library", None) or _default_rdna_library()
    if not os.path.exists(rdna_library):
        print(F"WARNING: rDNA library not found ({rdna_library}); "
              F"skipping rDNA identification", file=sys.stderr)
        return
    try:
        tc.identify_rdna(
            prefix, fasta, rdna_library, cpu,
            min_coverage=getattr(args, "rdna_min_coverage", 0.7),
            min_identity=getattr(args, "rdna_min_identity", 85.0),
        )
    except Exception as e:
        print(F"WARNING: rDNA identification failed: {e}", file=sys.stderr)


def annotation(prefix, library, gff=None, consensus_dir=None, cpu=1):
    """
    Run annotation on sequences defined in gff3 file based on coresponding
    consensus sequences in stored in consensu directory
    produce gff3 file with updated annotation  information
    :param prefix: prefix - base naame for input and output files
    :param library: library file for RepeatMasker
    :param gff: gff3 file with tidehunter results
    :param consensus_dir: directory with consensus sequences
    :param cpu: number of cpu cores to use
    :return:
    """
    if consensus_dir is None:
        consensus_dir = prefix + "_consensus"
    gff_short = None
    gff_short_annot = None
    if gff is None:
        gff = prefix + "_clustering.gff3"
        gff_short = prefix + "_tidehunter_short.gff3"
        gff_short_annot = prefix + "_tidehunter_short_annotation.gff3"
    gff_out = prefix + "_annotation.gff3"
    gff3_dir_split_files = prefix + "_annotation_split_files"
    # get list consensus sequences from consensus directory
    # naming scheme is TRC_10_dimer.fasta
    # use glob to get all files in directory
    consensus_files = glob.glob(consensus_dir + "/TRC*dimers.fasta")
    # it is possible that consensus_files does not exist
    if len(consensus_files) > 0:
        print(F"Annotating based on consensus sequences in {consensus_dir}")
        # conncatenate all consensus sequences to one file
        # Normally already written by clustering(); this covers `annotation`
        # run standalone over an output directory from an older version.
        consensus_files_concat = tc.write_consensus_sequences_all(consensus_dir)
        seq_lengths = tc.read_fasta_sequence_size(consensus_files_concat)
        
        # run RepeatMasker with automatic sequence name renaming
        rm_file = tc.run_repeatmasker_with_renaming(
            consensus_files_concat, library, cpu, 
            additional_params=F"-dir {consensus_dir}"
        )
        
        # parse RepeatMasker output
        rm_annotation = tc.get_repeatmasker_annotation(rm_file, seq_lengths, prefix)
        # add annotation to gff3 file, only if gff3 file exists
        if os.path.exists(gff):
            tc.add_attribute_to_gff(gff, gff_out, "Name", "annotation", rm_annotation)
            tc.split_gff3_by_cluster_name(gff_out, gff3_dir_split_files)
        else:
            print(F"gff3 file {gff} does not exist, no annotation added to gff3 file")
    else:
        # when consensus is not available, it is expected that gff3 file contains
        # consensus sequences, which can be annotated
        print("No consensus sequences found in {consensus_dir}")
        print(F"Annotating based on consensus sequences in stored in {gff}")
        tc.annotate_gff(gff, gff_out, library, cpu=cpu)
    if gff_short is not None:
        if os.path.exists(gff_short):
            print('Running annotation of omitted short regions from TideHunter')
            tc.annotate_gff(gff_short, gff_short_annot, library, cpu=cpu)


def clustering(fasta, prefix, gff3=None, min_length=None, dust=True, cpu=4,
               cluster_identity=75, cluster_coverage=0.8, resolve_overlaps=True):
    """
    Run clustering on sequences defined in gff3 file and fasta file
    produce gff3 file with cluster information
    :param fasta: fasta file with sequences
    :param prefix: prefix - base naame for input and output files
    :param gff3: gff3 file with tidehunter results
    :param min_length: minimal length of repeat to be included in clustering
    :param dust: use dust filter in blast search
    :param cpu: number of cpu cores to use
    :param cluster_identity: minimum BLASTN percent identity for a clustering edge
    :param cluster_coverage: minimum alignment coverage over the shorter sequence
    :return:

    """
    gff3_out = prefix + "_clustering.gff3"
    gff3_dir_split_files = prefix + "_clustering_split_files"
    fasta = fasta
    if gff3 is None:
        gff3 = prefix + "_tidehunter.gff3"
    if min_length is not None:
        print('running filtering on gff3 file')
        gff3 = tc.filter_gff_by_length(
                gff3,
                gff_short=prefix + "_tidehunter_short.gff3",
                min_length=min_length
                )

    # filtering on duplicates in gff3 file
    gff3 = tc.filter_gff_remove_duplicates(gff3)
    # check if gff3 has more than 1 sequence, it is enough to read
    # just beginning of the file
    with open(gff3, "r") as f:
        count = 0
        for i in f:
            if i.startswith("#"):
                continue
            count += 1
            if count > 1:
                break
    if count == 0:
        print("No tandem repeats found in gff3 file after filtering, exiting")
        exit(0)
    # get consensus sequences for clustering
    consensus_file = tempfile.NamedTemporaryFile(delete=False).name
    consensus_dimers_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(consensus_file, "w") as f, open(consensus_dimers_file, "w") as f2:
        # Only the consensus_sequence attribute (cons) and the ID are used here,
        # not the genomic sequence — load_sequence=False skips loading the whole
        # genome into RAM (fasta_to_dict), which otherwise OOMs on large genomes.
        for seq_id, seq, cons in tc.gff3_to_fasta(
                gff3, fasta, "consensus_sequence", load_sequence=False):
            mult = round(1 + 10000 / len(cons))
            consensus = cons * mult
            consensus_dimers = cons * 4
            if len(consensus) > 10000:
                consensus = consensus[0:10000]
            # write consensus sequence to file
            f.write(F">{seq_id}\n{consensus}\n")
            f2.write(F">{seq_id}\n{consensus_dimers}\n")
    # run dustmasker first, sequences which are completely masked
    # will not be used in clustering.
    mask_prop = tc.get_ssrs_proportions(consensus_dimers_file)
    # if count mask_prop above 0.9
    ssrs_id = [k for k, v in mask_prop.items() if v > 0.9]
    # remove ssrs from consensus sequences, they will be added back later
    # but not used for clustering
    ssrs_description = {}
    ssrs_seq = {}
    ssrs_dimers = {}
    # iterate over consensus_file and consensus_dimers_file
    consensus_file_filtered = tempfile.NamedTemporaryFile(delete=False).name
    consensus_dimers_file_filtered = tempfile.NamedTemporaryFile(delete=False).name
    with open(consensus_file, "r") as f, open(consensus_dimers_file, "r") as f2:
        for id, seq in tc.read_single_fasta_as_generator(f):
            if id not in ssrs_id:
                with open(consensus_file_filtered, "a") as f_out:
                    f_out.write(F">{id}\n{seq}\n")
    with open(consensus_dimers_file, "r") as f, open(consensus_dimers_file_filtered, "a") as f2:
        for id, seq in tc.read_single_fasta_as_generator(f):
            if id not in ssrs_id:
                f2.write(F">{id}\n{seq}\n")
            else:
                # NOTE if there are high proportion is simple
                # repeats, this could use a lot of memory!
                ssrs_dimers[id] = seq
                ssrs_description[id] = tc.get_ssrs_description(seq)
                ssrs_seq[id] = " ".join(
                        [i.split(" ")[0] for i in ssrs_description[id].split(
                                "\n"
                                )]
                        )

    # find unique ssrs seq
    ssrs_clusters = {}
    ssrs_representative = {}
    for k, ssrs in ssrs_seq.items():
        if ssrs not in ssrs_representative:
            ssrs_representative[ssrs] = k
    for k, ssrs in ssrs_seq.items():
        ssrs_clusters[k] = ssrs_representative[ssrs]
    # recalculate description for each ssrs_cluster
    dimers_ssrs_clusters = {}
    for n, repre_id in ssrs_clusters.items():
        if repre_id not in dimers_ssrs_clusters:
            dimers_ssrs_clusters[repre_id] = []
        dimers_ssrs_clusters[repre_id].append(ssrs_dimers[n])
    # recalculating ssrs description
    ssrs_cluster_description = {}
    for k, v in dimers_ssrs_clusters.items():
        ssrs_cluster_description[k] = tc.get_ssrs_description_multiple(v)
    if os.path.getsize(consensus_file_filtered) == 0:
        print("No tandem repeats left after dustmasking, skipping clustering")
    # first round of clustering by mmseqs2
    clusters1 = tc.find_cluster_by_mmseqs2(consensus_file_filtered, cpu=cpu)
    representative_id = set(clusters1.values())
    consensus_fasta_representative = tempfile.NamedTemporaryFile(delete=False).name
    tc.filter_fasta_file(consensus_dimers_file_filtered, consensus_fasta_representative,
                          representative_id
                          )
    # second round of clustering by blastn
    clusters2 = tc.find_clusters_by_blast_connected_component(
            consensus_fasta_representative, dust=dust, cpu=cpu,
            perc_identity=cluster_identity, min_coverage=cluster_coverage
            )
    # combine clusters
    clusters_final = clusters1.copy()

    for k, v in clusters1.items():
        if v in clusters2:
            clusters_final[k] = clusters2[v]
        else:
            clusters_final[k] = v

    # add ssrs_id back to clusters_final
    # and also to clusters1 - these are saved as well
    for k in ssrs_id:
        clusters_final[k] = ssrs_clusters[k]
        clusters1[k] = k
    # get total size of each cluster, store in dict
    cluster_size = tc.get_cluster_size2(gff3, clusters_final)

    # representative id sorted by cluster size
    representative_id = sorted(cluster_size, key=cluster_size.get, reverse=True)
    # rename values in clusters dictionary
    cluster_names = {}
    for i, v in enumerate(representative_id):
        cluster_names[v] = F"TRC_{i + 1}"

    ssrs_info = {}  # store ssrs info for gff3 file
    for k, v in clusters_final.items():
        clusters_final[k] = cluster_names[v]
        if k in ssrs_id:
            ssrs_info[cluster_names[v]] = "SSR"
        else:
            ssrs_info[cluster_names[v]] = "TR"
    #
    ssrs_description_final = {}  # with TRC names
    for k, v in ssrs_cluster_description.items():
        ssrs_description_final[cluster_names[k]] = v

    cons_cls, cons_cls_dimer = tc.add_cluster_info_to_gff3(gff3, gff3_out, clusters_final)

    tc.merge_overlapping_gff3_intervals(gff3_out, gff3_out)

    gff_tmp = gff3_out + "_tmp"

    tc.add_attribute_to_gff(gff3_out, gff_tmp, "Name", "repeat_type", ssrs_info)
    os.rename(gff_tmp, gff3_out)
    tc.add_attribute_to_gff(gff3_out, gff_tmp, "Name", "ssr", ssrs_description_final)
    os.rename(gff_tmp, gff3_out)

    # Make the clustering GFF3 non-overlapping across TRCs: variant arrays of a
    # satellite (e.g. rDNA) can be clustered into separate TRCs that interleave
    # and overlap at their boundaries; TideCluster's aim is to annotate each
    # region once. Each contested span goes to the dominant TRC (largest total
    # array length). Disabled with --keep_overlaps.
    if resolve_overlaps:
        tc.resolve_trc_overlaps(gff3_out, gff3_out)

    # save also first round of clustering for debugging
    cons_cls1, cons_cls_dimer1_ = tc.add_cluster_info_to_gff3(
            gff3, gff3_out + "_1.gff3", clusters1
            )
    tc.merge_overlapping_gff3_intervals(gff3_out + "_1.gff3", gff3_out + "_1.gff3")
    # split gff3 file to parts by cluster
    tc.split_gff3_by_cluster_name(gff3_out, gff3_dir_split_files)

    # for debugging
    # write consensus sequences by clusters to directory
    #  used gff3_out as base name for directory
    consensus_dir = prefix + "_consensus"
    tc.save_consensus_files(consensus_dir, cons_cls, cons_cls_dimer )   # this is used
    # for comparative analysis later
    # The concatenated pool is a clustering artefact -- it is exactly these
    # dimer files joined. It used to be built inside annotation() (its first
    # consumer was RepeatMasker), which left runs made without -l/--library
    # unable to feed tc_comparative_analysis.R at all. force=True because the
    # dimer files were just (re)written above.
    tc.write_consensus_sequences_all(consensus_dir, force=True)

    # consensus_dir = prefix + "_consensus_1"  # this was just for debugging
    # tc.save_consensus_files(consensus_dir, cons_cls1, cons_cls_dimer1_)
    # remove all temporary files
    os.remove(consensus_file)
    os.remove(consensus_dimers_file)
    os.remove(consensus_file_filtered)
    os.remove(consensus_dimers_file_filtered)
    os.remove(consensus_fasta_representative)


def tidehunter(fasta, tidehunter_arguments, prefix, cpu=4, max_memory=None):
    """
    run tidehunter on fasta file
    :param fasta: file with sequences
    :param tidehunter_arguments: tidehunter arguments
    :param prefix: prefix - base name for input and output files
    :param cpu: number of cpu cores to use
    :param max_memory: memory limit in GB (--max_memory); sizes the worker pool
    :return:

    """
    # get size of input file
    chunk_size = 500000
    overlap = 50000
    output = prefix + "_tidehunter.gff3"
    output_chunks = prefix + "_chunks.bed"
    # check is tidehunter _arguments contain specification for
    # number of threads to use - format is -t <number>
    # if not add it from cpu variable
    if " -t" not in tidehunter_arguments:
        tidehunter_arguments += F" -t {cpu}"

    # this fill split sequences to chunk and all is stored in single file
    fasta_file_chunked, matching_table = tc.split_fasta_to_chunks(
            fasta, chunk_size, overlap
            )
    results = tc.run_tidehunter(
            fasta_file_chunked, tidehunter_arguments, max_memory=max_memory
            )
    # O(1) chunk-token -> row index so per-feature coordinate remap is a dict
    # lookup, not a linear scan of the (up to ~genome_bp/chunk_size) table.
    token_index = tc.build_matching_table_token_index(matching_table)
    with open(output, "w") as out:
        # write GFF3 header
        out.write("##gff-version 3\n")

        with open(results) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                feature = tc.TideHunterFeature(line)
                # TideHunter is also returning sequences of NNN - to not include them
                # in the output
                if feature.consensus == "N" * feature.cons_length:
                    continue
                feature.recalculate_coordinates(matching_table, token_index)
                out.write(feature.gff3() + "\n")
    # clean up
    os.remove(fasta_file_chunked)
    os.remove(results)

    with open(output_chunks, "w") as out:
        for m in matching_table:
            out.write(F'{m[0]}\t{m[2]}\t{m[3]}\t{m[4]}\n')


def parse_tidehunter_results_to_gff3(results_file, matching_table, round_num):
    """
    Parse TideHunter results file into GFF3 features with round suffix.

    :param results_file: path to TideHunter output file
    :param matching_table: matching table for coordinate recalculation
    :param round_num: round number for suffix (1, 2, or 3)
    :return: list of GFF3Feature objects
    """
    gff3_list = []
    round_suffix = F"_rnd{round_num}"
    # O(1) chunk-token -> row index (see build_matching_table_token_index)
    token_index = tc.build_matching_table_token_index(matching_table)

    with open(results_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            feature = tc.TideHunterFeature(line)
            if feature.consensus == "N" * feature.cons_length:
                continue
            feature.recalculate_coordinates(matching_table, token_index)
            feature.repeat_ID = feature.repeat_ID + round_suffix
            gff3_list.append(feature)

    return gff3_list


def save_gff3_to_file(gff3_list, filepath):
    """
    Save GFF3 features to file with header.

    :param gff3_list: list of GFF3Feature objects
    :param filepath: output file path
    """
    with open(filepath, "w") as f:
        f.write("##gff-version 3\n")
        for feature in gff3_list:
            f.write(feature.gff3() + "\n")


def concat_gff3_files(gff3_files, out_path):
    """
    Concatenate per-round GFF3 files into ``out_path``: one ``##gff-version 3``
    header followed by every non-comment line, in the given file order.

    Byte-identical to ``save_gff3_to_file`` applied to the concatenation of the
    same features (each input file was written by ``save_gff3_to_file``), so it
    replaces holding all rounds' features in memory
    (``all_features_for_masking``) with a streamed on-disk merge.

    :param gff3_files: ordered list of GFF3 file paths
    :param out_path: output GFF3 path
    """
    with open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for fp in gff3_files:
            with open(fp) as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    out.write(line)


def run_tidehunter_round(fasta_input, tidehunter_args, chunk_size, overlap,
                         round_num, prefix, cpu, keep_rounds=False,
                         max_memory=None):
    """
    Run a single round of TideHunter analysis.

    :param fasta_input: input FASTA file (may be masked)
    :param tidehunter_args: TideHunter command line arguments
    :param chunk_size: size of chunks for splitting
    :param overlap: overlap size between chunks
    :param round_num: round number (1, 2, or 3)
    :param prefix: output prefix for debugging files
    :param cpu: number of CPUs
    :param keep_rounds: if True, save permanent copy of round results
    :param max_memory: memory limit in GB (--max_memory); sizes the worker pool
    :return: tuple of (gff3_features_list, temp_gff3_file_path, matching_table)
    """
    # Add CPU threads if not specified
    if " -t" not in tidehunter_args:
        tidehunter_args += F" -t {cpu}"

    # Split FASTA and run TideHunter
    fasta_file_chunked, matching_table = tc.split_fasta_to_chunks(
        fasta_input, chunk_size, overlap
    )
    results = tc.run_tidehunter(fasta_file_chunked, tidehunter_args,
                                max_memory=max_memory)

    # Parse results into GFF3
    gff3_list = parse_tidehunter_results_to_gff3(results, matching_table, round_num)

    # Write to temporary file
    temp_gff3_file = tempfile.NamedTemporaryFile(delete=False, suffix=".gff3").name
    save_gff3_to_file(gff3_list, temp_gff3_file)

    # Save permanent copy if requested
    if keep_rounds:
        permanent_gff3_file = F"{prefix}_tidehunter_round{round_num}.gff3"
        save_gff3_to_file(gff3_list, permanent_gff3_file)
        print(f"Saved round {round_num} results to: {permanent_gff3_file}")

    print(f"Round {round_num} complete: {len(gff3_list)} features found")

    # Clean up intermediate files
    os.remove(fasta_file_chunked)
    os.remove(results)

    return gff3_list, temp_gff3_file, matching_table


def tidehunter_long(fasta, prefix, cpu=4, keep_rounds=False, max_memory=None):
    """
    Run TideHunter in three rounds with increasing monomer size ranges.
    Results from each round are used to mask the sequence for the next round.

    Round 1: -p 40 -P 3000 -c 5 -e 0.25 (default long monomers)
    Round 2: -p 3001 -P 10000 -c 5 -e 0.25 (medium-long monomers, masked for round 1 results)
    Round 3: -p 10001 -P 25000 -c 5 -e 0.25 (very long monomers, masked for rounds 1-2 results)

    All three outputs are merged into a single GFF3 file.

    :param fasta: file with sequences
    :param prefix: prefix - base name for input and output files
    :param cpu: number of cpu cores to use
    :param keep_rounds: if True, keep intermediate GFF3 files from each round for debugging
    :param max_memory: memory limit in GB (--max_memory); sizes the worker pool.
        Per-part peak RSS grows ~3x from round 1 to round 3, so a pool size that
        is safe early is fatal late -- the gate re-measures per round, but only
        protects the run if the budget is the job's and not the host's.
    :return: None
    """
    print("Starting TideHunter long analysis with 3 rounds")
    if keep_rounds:
        print("Intermediate round files will be kept for debugging")

    chunk_size = 500000
    overlap = 50000
    output = prefix + "_tidehunter.gff3"
    output_chunks = prefix + "_chunks.bed"

    # Define rounds configuration: (description, args, input_fasta, mask_file)
    rounds_config = [
        ("ROUND 1: Long monomers (p=40-3000)", "-p 40 -P 3000 -c 5 -e 0.25", fasta, None),
        ("ROUND 2: Medium-long monomers (p=3001-10000, masked for round 1)", "-p 3001 -P 10000 -c 5 -e 0.25", None, None),
        ("ROUND 3: Very long monomers (p=10001-25000, masked for rounds 1-2)", "-p 10001 -P 25000 -c 5 -e 0.25", None, None),
    ]

    round_counts = []       # per-round feature counts (for the summary only)
    temp_gff3_files = []     # per-round GFF3 files on disk (features live here)

    for round_num, (description, tidehunter_args, input_fasta, mask_file) in enumerate(rounds_config, 1):
        print(f"\n=== {description} ===")

        # Determine input FASTA for this round
        if round_num == 1:
            current_fasta = fasta
        else:
            # Masking input = all previous rounds' features, streamed from their
            # per-round GFF3 files on disk (no in-RAM accumulation of every
            # feature across rounds). Byte-identical to the previous
            # save_gff3_to_file(all_features_for_masking, ...).
            merged_gff3_file = tempfile.NamedTemporaryFile(delete=False, suffix=".gff3").name
            concat_gff3_files(temp_gff3_files, merged_gff3_file)
            current_fasta = tc.mask_fasta_with_gff3(fasta, merged_gff3_file)
            os.remove(merged_gff3_file)

        # Run the round
        gff3_list, temp_gff3_file, matching_table = run_tidehunter_round(
            current_fasta, tidehunter_args, chunk_size, overlap,
            round_num, prefix, cpu, keep_rounds, max_memory=max_memory
        )

        temp_gff3_files.append(temp_gff3_file)
        round_counts.append(len(gff3_list))
        # gff3_list is not retained across rounds; its features are in
        # temp_gff3_file on disk, so peak RAM is one round's parse, not all three.

        # Clean up masked FASTA (except for round 1)
        if round_num > 1:
            os.remove(current_fasta)

    # Merge all three rounds into final output by streamed on-disk concatenation
    print("\n=== Merging all three rounds ===")
    concat_gff3_files(temp_gff3_files, output)

    # Clean up temporary GFF3 files (unless keep_rounds is enabled)
    if not keep_rounds:
        for temp_file in temp_gff3_files:
            os.remove(temp_file)

    # Write chunks BED file (use matching table from first round, as it covers full genome)
    with open(output_chunks, "w") as out:
        for m in matching_table:
            out.write(F'{m[0]}\t{m[2]}\t{m[3]}\t{m[4]}\n')

    total_features = sum(round_counts)
    round_stats = ", ".join([f"Round {i}: {round_counts[i-1]}" for i in range(1, 4)])
    print(f"\nTideHunter long analysis complete: {total_features} total features found")
    print(round_stats)


def save_args_to_file(args):
    args_file = f"{args.prefix}_cmd_args.json"
    # load any previously saved args
    saved = {}
    if os.path.exists(args_file):
        with open(args_file, "r") as f:
            saved = json.load(f)
    # update with current args (only those not None)
    for k, v in vars(args).items():
        if v is not None:
            saved[k] = v
    # write back merged args
    with open(args_file, "w") as f:
        json.dump(saved, f)

def load_args_from_file(prefix):
    args_file = f"{prefix}_cmd_args.json"
    if os.path.exists(args_file):
        with open(args_file, "r") as f:
            return json.load(f)
    return {}


def validate_threshold_args(args):
    """Validate the user-tunable clustering / superfamily thresholds.

    Hard errors abort runs that would silently yield empty or meaningless
    groupings (for example a coverage fraction mistakenly entered as a
    percent); warnings flag biologically unusual but permissible settings.
    Parameters absent from the active subcommand are skipped.
    """
    errors, warnings = [], []

    ci = getattr(args, "cluster_identity", None)
    if ci is not None:
        if not 0 < ci <= 100:
            errors.append(
                F"--cluster_identity must be a percent in (0, 100]; got {ci}")
        elif ci < 50:
            warnings.append(
                F"--cluster_identity={ci} is very low; unrelated arrays may be "
                "merged into one TRC")

    cc = getattr(args, "cluster_coverage", None)
    if cc is not None:
        if not 0 < cc <= 1:
            errors.append(
                F"--cluster_coverage must be a fraction in (0, 1] "
                F"(e.g. 0.8, not 80); got {cc}")
        elif cc < 0.3:
            warnings.append(
                F"--cluster_coverage={cc} is very low; arrays sharing only a "
                "short region may be clustered together")

    ss = getattr(args, "superfamily_score", None)
    if ss is not None:
        if ss < 0:
            errors.append(F"--superfamily_score must be >= 0; got {ss}")
        elif ss == 0:
            warnings.append(
                "--superfamily_score=0 collapses any TRC pair with a BLAST hit "
                "into a single superfamily")
        elif ss > 100:
            warnings.append(
                F"--superfamily_score={ss} exceeds the practical maximum (~100) "
                "and will likely prevent all superfamily grouping")

    mm = getattr(args, "max_memory", None)
    if mm is not None:
        if mm <= 0:
            errors.append(F"--max_memory must be > 0 GB; got {mm}")
        elif mm < 4:
            warnings.append(
                F"--max_memory={mm} GB is below one TideHunter part's typical "
                "peak; pools will run serially")

    for w in warnings:
        print(F"WARNING: {w}", file=sys.stderr)
    if errors:
        for e in errors:
            print(F"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-v", "--version", action="version", version=__version__)

    subparsers = parser.add_subparsers(dest='command', help='TideHunter wrapper')

    # TideHunter
    parser_tidehunter = subparsers.add_parser(
            'tidehunter', help='Run wrapper of TideHunter'
            )

    parser_tidehunter.add_argument(
            '-f', '--fasta', type=str, required=True,
            help='Path to reference sequence in fasta format (gzipped files supported)'
            )
    parser_tidehunter.add_argument(
            '-pr', '--prefix', type=str, required=True, help='Base name for output files'
            )
    parser_tidehunter.add_argument(
            "-c", "--cpu", type=int, default=4,
            help="Number of CPUs to use"
            )

    # Create mutually exclusive group for TideHunter arguments vs long mode
    tidehunter_mode = parser_tidehunter.add_mutually_exclusive_group()
    tidehunter_mode.add_argument(
            '-T', '--tidehunter_arguments', type=str, nargs="?", required=False,
            default="-p 40 -P 3000 -c 5 -e 0.25",
            help=('additional arguments for TideHunter in quotes'
                  ', default value: %(default)s)'),
            )
    tidehunter_mode.add_argument(
            "--long", action="store_true", required=False, default=False,
            help="Run TideHunter in three rounds with increasing monomer sizes (40-25000 nt)"
            )
    parser_tidehunter.add_argument(
            "--keep_rounds", action="store_true", required=False, default=False,
            help="Keep intermediate GFF3 files from each round for debugging (only with --long)"
            )
    parser_tidehunter.add_argument(
            "--max_memory", "--max-memory", type=float, default=None,
            metavar="GB",
            help=MAX_MEMORY_HELP
            )
    # Clustering
    parser_clustering = subparsers.add_parser(
            'clustering', help='Run clustering on TideHunter output'
            )
    parser_clustering.add_argument(
            "-f", "--fasta", help="Reference fasta (gzipped files supported)", required=True
            )

    parser_clustering.add_argument(
            "-m", "--min_length", help="Minimum length of tandem repeat array to be "
                                       "included in clustering step. Shorter arrays are "
                                       "discarded, default (%(default)s)",
            required=False,
            default=5000, type=int
            )

    parser_clustering.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF file from previous tidehunter "
                                     "step"),
            required=True
            )
    parser_clustering.add_argument(
            "-g", "--gff", help=("GFF3 output file from tidehunter step. If not provided "
                                 "the file named 'prefix_tidehunter.gff3' will be used"),
            required=False, default=None
            )
    parser_clustering.add_argument(
            '-nd', '--no_dust', required=False, default=False, action='store_true',
            help='Do not use dust filter in blastn when clustering'
            )
    parser_clustering.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )
    parser_clustering.add_argument(
            "--cluster_identity", type=float, default=75,
            help=("Minimum BLASTN percent identity for an array-level (TRC) "
                  "clustering edge. Lower values give looser clusters. "
                  "Default (%(default)s)")
            )
    parser_clustering.add_argument(
            "--cluster_coverage", type=float, default=0.8,
            help=("Minimum alignment coverage over the shorter sequence for an "
                  "array-level (TRC) clustering edge. Default (%(default)s)")
            )
    parser_clustering.add_argument(
            "--keep_overlaps", action="store_true", default=False,
            help=("Do not resolve overlapping TRC regions. By default the "
                  "clustering GFF3 is made non-overlapping across TRCs (each "
                  "contested span assigned to the TRC with the largest total "
                  "array length).")
            )

    # Annotation
    parser_annotation = subparsers.add_parser(
            'annotation', help=('Run annotation on output from clustering step'
                                ' using reference library of tandem repeats')
            )
    parser_annotation.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF3 file from previous clustering "
                                     "step"),
            required=True
            )

    parser_annotation.add_argument(
            "-g", "--gff", help=("GFF3 output file from clustering step. If not provided "
                                 "the file named 'prefix_clustering.gff3' will be used"),
            required=False, default=None
            )
    parser_annotation.add_argument(
            "-cd", "--consensus_directory",
            help=("Directory with consensus sequences which are to be "
                  "annotated. If not provided the directory named 'prefix_consensus' "
                  "will be used"), required=False, default=None
            )

    parser_annotation.add_argument(
            "-l", "--library", help="Path to library of tandem repeats", required=True, )

    parser_annotation.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )

    # tarean
    parser_tarean = subparsers.add_parser(
            'tarean', help='Run TAREAN on clusters to extract representative sequences'
            )
    parser_tarean.add_argument(
            "-g", "--gff", help=("GFF3 output file from annotation or clustering step"
                                 "If not provided the file named "
                                 "'prefix_annotation.gff3' "
                                 "will be used instead. If 'prefix_annotation.gff3' is "
                                 "not "
                                 "found, 'prefix_clustering.gff3' will be used"
                                 ),
            required=False, default=None
            )
    parser_tarean.add_argument(
            "-f", "--fasta", help="Reference fasta (gzipped files supported)", required=True
            )
    parser_tarean.add_argument(
            "-pr", "--prefix", help=("Prefix is used as a base name for output files."
                                     "If --gff is not provided, prefix will be also used"
                                     "to identify GFF3 files from previous clustering/"
                                     "annotation step"),
            required=True
            )
    parser_tarean.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )
    parser_tarean.add_argument(
            "-M", "--min_total_length", type=int, default=50000,
            help=("Minimum combined length of tandem repeat arrays within a single "
                  "cluster, required for inclusion in TAREAN analysis."
                  "Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--max_memory", "--max-memory", type=float, default=None,
            metavar="GB",
            help=MAX_MEMORY_HELP
            )
    parser_tarean.add_argument(
            "--kite_rescore_max_period", type=int, default=10000,
            help=("kitehor `rescore --max-period` cap (bp). Peaks above this stay "
                  "NA in rescore output; TideCluster falls back to the top-scored "
                  "kite peak for those arrays. Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--kite_rescore_max_period_ext", type=int, default=25000,
            help=("Higher `rescore --max-period` for the selective long-period "
                  "re-search: arrays whose dominant monomer exceeds "
                  "--kite_rescore_max_period and have no confident founder below "
                  "it are re-rescored at this cap. Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--kite_rescore_top_n", type=int, default=20,
            help=("kitehor `rescore --top-n` cap per array (number of kite peaks "
                  "rescored). Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--superfamily_score", type=float, default=20,
            help=("Minimum BLASTN score, (alignment_length * percent_identity - "
                  "gap_openings) / longer_consensus_length, for a superfamily edge "
                  "between two TRC consensus sequences. Lower values give looser "
                  "superfamilies. Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--rdna_library", default=None,
            help=("rDNA reference library (RepeatMasker name#class format, "
                  "classes rDNA_45S/* and rDNA_5S/*) for rDNA identification. "
                  "Defaults to the bundled data/rdna_library.fasta.")
            )
    parser_tarean.add_argument(
            "--no_rdna", action="store_true", default=False,
            help="Disable rDNA (45S/5S) identification of TRCs."
            )
    parser_tarean.add_argument(
            "--rdna_min_coverage", type=float, default=0.7,
            help=("Minimum best-subunit reference coverage to call a TRC rDNA. "
                  "Default (%(default)s)")
            )
    parser_tarean.add_argument(
            "--rdna_min_identity", type=float, default=85.0,
            help=("Minimum percent identity for an rDNA reference hit to count. "
                  "Default (%(default)s)")
            )

    parser_run_all = subparsers.add_parser(
            'run_all', help='Run all steps of TideCluster'
            )

    parser_run_all.add_argument(
            "-f", "--fasta", help="Reference fasta (gzipped files supported)", required=True, type=str
            )
    parser_run_all.add_argument(
            "-pr", "--prefix", help="Base name used for input and output files",
            required=True, type=str
            )
    parser_run_all.add_argument(
            "-l", "--library", help="Path to library of tandem repeats", required=False,
            type=str, default=None
            )
    parser_run_all.add_argument(
            "-m", "--min_length", help=("Minimum length of tandem repeat"
                                        " (%(default)s)"), required=False,
            default=5000, type=int
            )
    parser_run_all.add_argument(
            "-nd", "--no_dust", help="Do not use dust filter in blastn when clustering",
            action="store_true", required=False, default=False
            )
    parser_run_all.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )

    # Create mutually exclusive group for TideHunter arguments vs long mode in run_all
    run_all_tidehunter_mode = parser_run_all.add_mutually_exclusive_group()
    run_all_tidehunter_mode.add_argument(
            '-T', '--tidehunter_arguments', type=str, nargs="?", required=False,
            default="-p 40 -P 3000 -c 5 -e 0.25",
            help=('additional arguments for TideHunter in quotes'
                  ', default value: %(default)s)'),
            )
    run_all_tidehunter_mode.add_argument(
            "--long", action="store_true", required=False, default=False,
            help="Run TideHunter in three rounds with increasing monomer sizes (40-25000 nt)"
            )

    parser_run_all.add_argument(
            "--keep_rounds", action="store_true", required=False, default=False,
            help="Keep intermediate GFF3 files from each round for debugging (only with --long)"
            )

    parser_run_all.add_argument(
            "--max_memory", "--max-memory", type=float, default=None,
            metavar="GB",
            help=MAX_MEMORY_HELP
            )

    parser_run_all.add_argument(
            "--cleanup", action="store_true", default=False,
            help=("Delete intermediate files once the run finishes successfully "
                  "(~80%% of the output directory). The report, its re-render "
                  "inputs, the comparative-analysis inputs and the per-TRA "
                  "consensus inputs are all kept; --keep_rounds files are never "
                  "touched. Nothing is deleted if any step fails.")
            )

    parser_run_all.add_argument(
            "-M", "--min_total_length", type=int, default=50000,
            help=("Minimum combined length of tandem repeat arrays within a single "
                  "cluster, required for inclusion in TAREAN analysis."
                  "Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--kite_rescore_max_period", type=int, default=10000,
            help=("kitehor `rescore --max-period` cap (bp). Peaks above this stay "
                  "NA in rescore output; TideCluster falls back to the top-scored "
                  "kite peak for those arrays. Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--kite_rescore_max_period_ext", type=int, default=25000,
            help=("Higher `rescore --max-period` for the selective long-period "
                  "re-search: arrays whose dominant monomer exceeds "
                  "--kite_rescore_max_period and have no confident founder below "
                  "it are re-rescored at this cap. Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--kite_rescore_top_n", type=int, default=20,
            help=("kitehor `rescore --top-n` cap per array (number of kite peaks "
                  "rescored). Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--cluster_identity", type=float, default=75,
            help=("Minimum BLASTN percent identity for an array-level (TRC) "
                  "clustering edge. Lower values give looser clusters. "
                  "Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--cluster_coverage", type=float, default=0.8,
            help=("Minimum alignment coverage over the shorter sequence for an "
                  "array-level (TRC) clustering edge. Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--keep_overlaps", action="store_true", default=False,
            help=("Do not resolve overlapping TRC regions. By default the "
                  "clustering GFF3 is made non-overlapping across TRCs (each "
                  "contested span assigned to the TRC with the largest total "
                  "array length).")
            )
    parser_run_all.add_argument(
            "--superfamily_score", type=float, default=20,
            help=("Minimum BLASTN score, (alignment_length * percent_identity - "
                  "gap_openings) / longer_consensus_length, for a superfamily edge "
                  "between two TRC consensus sequences. Lower values give looser "
                  "superfamilies. Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--rdna_library", default=None,
            help=("rDNA reference library (RepeatMasker name#class format, "
                  "classes rDNA_45S/* and rDNA_5S/*) for rDNA identification. "
                  "Defaults to the bundled data/rdna_library.fasta.")
            )
    parser_run_all.add_argument(
            "--no_rdna", action="store_true", default=False,
            help="Disable rDNA (45S/5S) identification of TRCs."
            )
    parser_run_all.add_argument(
            "--rdna_min_coverage", type=float, default=0.7,
            help=("Minimum best-subunit reference coverage to call a TRC rDNA. "
                  "Default (%(default)s)")
            )
    parser_run_all.add_argument(
            "--rdna_min_identity", type=float, default=85.0,
            help=("Minimum percent identity for an rDNA reference hit to count. "
                  "Default (%(default)s)")
            )

    parser_rdna = subparsers.add_parser(
            'rdna', help=('Identify rDNA (45S/5S) TRCs in an existing run and '
                          'add rDNA_type/rDNA_coverage to the clustering GFF3 + '
                          'report (re-runnable, e.g. with an extended library)')
            )
    parser_rdna.add_argument(
            "-pr", "--prefix", required=True,
            help="Prefix of an existing TideCluster run (e.g. from run_all)."
            )
    parser_rdna.add_argument(
            "-f", "--fasta", required=True,
            help="Reference fasta (gzipped supported); used for the genomic "
                 "fallback when a TRC has no usable consensus."
            )
    parser_rdna.add_argument(
            "--rdna_library", default=None,
            help="rDNA reference library; defaults to bundled data/rdna_library.fasta."
            )
    parser_rdna.add_argument(
            "-c", "--cpu", type=int, default=4, help="Number of CPUs to use"
            )
    parser_rdna.add_argument(
            "--rdna_min_coverage", type=float, default=0.7,
            help="Minimum best-subunit reference coverage. Default (%(default)s)"
            )
    parser_rdna.add_argument(
            "--rdna_min_identity", type=float, default=85.0,
            help="Minimum percent identity for an rDNA hit. Default (%(default)s)"
            )

    parser.description = """Wrapper of TideHunter
    This script enable to run TideHunter on large fasta files in parallel. It splits
    fasta file into chunks and run TideHunter on each chunk. Identified tandem repeat 
    are then clustered, annotated and representative consensus sequences are extracted.
    
     
    """

    # make epilog, in epilog keep line breaks as preformatted text

    parser.epilog = ('''
    Example of usage:
    
    # first run tidehunter on fasta file to generate raw GFF3 output
    # TideCluster.py tidehunter -c 10 -f test.fasta -pr prefix 
    
    # then run clustering on the output from previous step to cluster similar tandem 
    repeats
    TideCluster.py clustering -c 10 -f test.fasta -pr prefix -m 5000
    
    # then run annotation on the clustered output to annotate clusters with reference
    # library of tandem repeats in RepeatMasker format
    TideCluster.py annotation -c 10 -pr prefix -l library.fasta
    
    # then run TAREAN on the annotated output to extract representative consensus
    # and generate html report
    TideCluster.py tarean -c 10 -f test.fasta -pr prefix
    
    Recommended parameters for TideHunter:
    short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
    long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25" (default)
    
    The increasing -p and -P values can be used to target longer monomers but it will lead
    to increased computational time. If you need to target monomers for up to 25000 bp,
    it is recommended to use --long option which will run TideHunter in three rounds
    with increasing monomer size ranges (40-3000, 3001-10000, 10001-25000). After each round,
    identified tandem repeats are masked in the input sequences for the next round. This
    approach improves detection of long monomers while keeping computational time manageable.
    
    For parallel processing include -c option before command name. 
    
    For more information about TideHunter parameters see TideHunter manual.
    
    Library of tandem repeats for annotation step are sequences in RepeatMasker format
    where header is in format:
    
    >id#clasification
    
    ''')

    cmd_args = parser.parse_args()
    validate_threshold_args(cmd_args)

    # Handle gzipped input FASTA files
    cleanup_fasta = None
    if hasattr(cmd_args, 'fasta') and cmd_args.fasta:
        # Preserve the original FASTA path before modification for reporting
        cmd_args.original_fasta = cmd_args.fasta
        cmd_args.fasta, cleanup_fasta = tc.prepare_fasta_input(cmd_args.fasta)

    # Wrap execution in try-finally to ensure cleanup
    try:
        save_args_to_file(cmd_args)
        if cmd_args.command == "tidehunter":
            # Check if --long flag is set
            if hasattr(cmd_args, 'long') and cmd_args.long:
                keep_rounds = getattr(cmd_args, 'keep_rounds', False)
                tidehunter_long(
                        cmd_args.fasta, cmd_args.prefix,
                        cmd_args.cpu, keep_rounds=keep_rounds,
                        max_memory=cmd_args.max_memory
                        )
            else:
                tidehunter(
                        cmd_args.fasta, cmd_args.tidehunter_arguments, cmd_args.prefix,
                        cmd_args.cpu, max_memory=cmd_args.max_memory
                        )
        elif cmd_args.command == "clustering":
            clustering(
                    cmd_args.fasta, cmd_args.prefix, cmd_args.gff, cmd_args.min_length,
                    not cmd_args.no_dust, cmd_args.cpu,
                    cluster_identity=cmd_args.cluster_identity,
                    cluster_coverage=cmd_args.cluster_coverage,
                    resolve_overlaps=not cmd_args.keep_overlaps
                    )
        elif cmd_args.command == "annotation":
            annotation(
                    cmd_args.prefix, cmd_args.library, cmd_args.gff,
                    cmd_args.consensus_directory,
                    cmd_args.cpu
                    )
        elif cmd_args.command == "tarean":
            tarean(
                    prefix=cmd_args.prefix,
                    gff=cmd_args.gff,
                    fasta=cmd_args.fasta,
                    cpu=cmd_args.cpu,
                    min_total_length=cmd_args.min_total_length,
                    args=cmd_args,
                    version=__version__,
                    max_memory=cmd_args.max_memory
                    )
        elif cmd_args.command == "run_all":
            # Check if --long flag is set for run_all
            if hasattr(cmd_args, 'long') and cmd_args.long:
                keep_rounds = getattr(cmd_args, 'keep_rounds', False)
                tidehunter_long(
                        cmd_args.fasta, cmd_args.prefix,
                        cmd_args.cpu, keep_rounds=keep_rounds,
                        max_memory=cmd_args.max_memory
                        )
            else:
                tidehunter(
                        cmd_args.fasta, cmd_args.tidehunter_arguments, cmd_args.prefix,
                        cmd_args.cpu, max_memory=cmd_args.max_memory
                        )
            clustering(
                    cmd_args.fasta, cmd_args.prefix,
                    min_length=cmd_args.min_length,
                    dust=not cmd_args.no_dust,
                    cpu=cmd_args.cpu,
                    cluster_identity=cmd_args.cluster_identity,
                    cluster_coverage=cmd_args.cluster_coverage,
                    resolve_overlaps=not cmd_args.keep_overlaps
                    )
            if cmd_args.library:
                annotation(
                    cmd_args.prefix, cmd_args.library,
                    cpu=cmd_args.cpu
                    )
            tarean(
                    prefix=cmd_args.prefix,
                    fasta=cmd_args.fasta,
                    gff=None,
                    cpu=cmd_args.cpu,
                    min_total_length=cmd_args.min_total_length,
                    args=cmd_args,
                    version=__version__,
                    max_memory=cmd_args.max_memory
                    )
            # Last thing in the run, and only here: every step above raises on
            # failure, so reaching this line means the run succeeded. Deliberately
            # NOT in the `finally:` below -- a failed run's intermediates are what
            # you need to diagnose it.
            if cmd_args.cleanup:
                _cleanup_outputs(cmd_args.prefix)
        elif cmd_args.command == "rdna":
            _maybe_identify_rdna(cmd_args.prefix, cmd_args.fasta, cmd_args,
                                 cmd_args.cpu)
            _build_report_v2(cmd_args.prefix)

        else:
            parser.print_help()
            sys.exit(1)
    finally:
        # Cleanup temporary uncompressed FASTA if it was created
        if cleanup_fasta:
            cleanup_fasta()
