# Specification — optional output cleanup (`--cleanup`)

**Status:** IMPLEMENTED (`tc_utils.CLEANUP_PATTERNS` / `cleanup_run_directory`,
`TideCluster.py --cleanup`; tests in `tests/test_cleanup.py` and `tests/long.sh`)
· **Date:** 2026-08-21
**Written against:** TideCluster 1.20.1
**Companion:** `docs/carp_cleanup_alignment_request.md` (the CARP-side FR)

## Problem

A TideCluster run directory is dominated by intermediates. On a real run
(*Solanum lycopersicum*: 1.2 GB, 2328 files) **82.3 % of the bytes are files
nothing reads again.** At the scale TideCluster is now used — 14 concurrent
~4 Gbp *Pisum* assemblies, and thousands of genomes through CARP — that is the
difference between archiving a batch of runs and deleting them.

Downstream consumers already work around this by deleting whole trees. CARP's
`cleanup_intermediates: maximal` removes `TideCluster_{tarean,kite,consensus}/`,
which frees the disk but silently destroys three capabilities (see the companion
FR). TideCluster should offer the cleanup itself, with a contract about what
survives.

## Goal

One optional, off-by-default flag that deletes intermediates at the end of a
successful `run_all`, and a **written guarantee** about what the pruned directory
can still do.

### Non-goals

- Cleanup levels or per-tree selection. One defined set, on or off.
- Cleanup for individual subcommands (`tidehunter`, `clustering`, `annotation`,
  `tarean`). `run_all` only — a subcommand run is by definition part of a
  workflow whose next step we cannot see.
- Deleting anything the user explicitly asked for (see `--keep_rounds` below).
- Reclaiming space from runs that already finished. `--cleanup` acts on the run
  it is part of. (A standalone tool was considered and dropped; if the 14 *Pisum*
  directories need pruning, the purge set below is a `find -delete` one-liner.)

## The contract

After `run_all --cleanup`, **all four of these must still work.** This is the
specification's core claim and every purge decision below is justified against
it.

| # | Capability | Verified by |
|---|---|---|
| G1 | Viewing the HTML report — `<prefix>_index.html`, `<prefix>_report/`, `<prefix>_report_legacy/` — with no broken links | `tests/report_linkcheck.py` (already exists) |
| G2 | Using the run as **comparative-analysis input** — all six files in the README's comparative input table resolve | new check (§Testing) |
| G3 | **Re-rendering** the report with `tc_rerender_report.py` | new check (§Testing) |
| G4 | **Per-TRA consensus** (`tc_per_tra_consensus.py`) and reannotation inputs | new check (§Testing) |

G3 is the binding constraint. Because `tc_rerender_report.py` reads from inside
`<prefix>_kite/` and `<prefix>_tarean/`, **no tree can be deleted wholesale** —
the purge is file-level within them. This is the single biggest difference from
CARP's `maximal`, and it costs surprisingly little: 82.3 % versus roughly 95 %.

## The purge set

Paths are relative to the run directory; `<prefix>` is the run's `-pr` value.
Sizes are from the *Solanum* run above.

| # | Pattern | Size | % | Justification |
|---|---|---|---|---|
| P1 | `<prefix>_kite/kitehor.periodogram` | 444 MB | 38.7 % | Consumed only by `kite_heatmaps.R`, during the run, to render `<prefix>_kite/profile_plots/*.png`. The report reads the PNGs and `monomer_size_top3_estimats.csv`, never the periodogram. |
| P2 | `<prefix>_tarean/*/*.kmers` | 246 MB | 21.4 % | Written by `tarean/kmer_counting.py`; no reader anywhere in the codebase. |
| P3 | `<prefix>_tarean/*/ggmin.RData`, `<prefix>_tarean/*/monomers.RData` | 168 MB | 14.6 % | **Write-only.** `save()`d at `tarean/methods.R:654-655`; no `load()` of either exists in the repository. |
| P4 | `<prefix>_tarean/*.fasta_tarean/TRC_*.fasta` | 76 MB | 6.6 % | Byte-identical copy of `<prefix>_tarean/fasta/TRC_*.fasta`, which is **kept** (G4). Verified identical with `cmp`. |
| P5 | `<prefix>_consensus/*_renamed.fasta`, `*_renamed.fasta.cat`, `*_renamed.fasta.masked` | 12 MB | 1.0 % | RepeatMasker leftovers from the annotation step. |
| P6 | `<prefix>_kite/_kite_input_longext.fasta`, `<prefix>_kite/_rescored_longext.*` | 0.45 MB | ~0 % | Intermediates of the selective long-period re-search; already `_`-prefixed to mark them as internal. |
| P7 | `<prefix>_clustering.gff3_1.gff3` | 72 KB | ~0 % | The mmseqs2-stage intermediate, already documented as such in the README. |
| | **total** | **~946 MB** | **82.3 %** | 1.2 GB → 204 MB |

P2/P3/P4 are globs over per-TRC directories; they must not match
`<prefix>_tarean/fasta/` (P4's glob is anchored on `*.fasta_tarean/` precisely
for that reason).

## What is kept, and why

The 204 MB residual, largest first:

| Kept | Size | Required by |
|---|---|---|
| `<prefix>_tarean/fasta/TRC_*.fasta` | 80 MB | **G4 only** — `tc_per_tra_consensus.py` reads this directory |
| `<prefix>_report/` (incl. `data/report.json` 23 MB) | 29 MB | G1 |
| `<prefix>_tarean/*/{report.html,img/,*.csv}` | 56 MB | G3 (re-render re-vendors from here) and G1 for the pre-vendored copies |
| `<prefix>_tidehunter*.gff3` | 18 MB | G4 (`_tidehunter.gff3`); the `_short*` variants are documented outputs |
| `<prefix>_consensus/consensus_sequences_all.fasta`, `TRC_*_dimers.fasta` | 6 MB | G2 |
| `<prefix>_kite/kitehor.*.tsv`, `monomer_size_*.csv`, `profile_plots/` | 7 MB | G3; the small kitehor TSVs also keep `tc_utils.build_monomer_size_csv` re-runnable (`tools/tc_regen.sh`). **`monomer_size_top3_estimats.csv` is additionally read by CARP** (`make_repeat_report.R:646`, `make_summary_plots.R:192`) to label families with their monomer length — it must not be dropped from the keep set later. |
| top-level GFF3/TSV/CSV/JSON side-cars | small | G1–G4 |

Two notes worth recording:

- **`<prefix>_tarean/fasta/` is 39 % of the residual and is kept for G4 alone.**
  Dropping G4 would take the run to ~116 MB (≈90 % freed) at this scale. Settled:
  CARP measured the same trade on a 94 Gbp run and it is 0.61 GB of 44.59 GB
  (1.4 %) there — the lever shrinks as the genome grows, because the residual is
  dominated by array FASTAs only at small scale. G4 stays in the contract.
- **Per-TRC TAREAN PNGs exist twice** — in `<prefix>_tarean/*/img/` (~31 MB here)
  and vendored into `<prefix>_report/img/` since 1.17.0. Both are kept: the
  vendored copies serve G1, the originals serve G3, since a re-render re-vendors
  from source and would otherwise emit a report with missing images. Deduplicating
  these (e.g. re-vendoring from the vendored copy) is a possible follow-up, not
  part of this spec.

## Validation of this specification

The purge set and the contract were checked against each other on the *Solanum*
run before any code was written, by enumerating what the globs match and what the
four guarantees require, and intersecting the two sets:

```
G2 comparative:   6 required files present
G4 per-TRA:      76 required files present
G3 re-render:   574 required files present
G1 report:       97 required files present

purge: 110 files, 943 MB of 1147 MB (82.3%)  -> residual 203 MB
CONFLICTS: 0
```

Zero of the 753 guaranteed files are matched by any purge pattern. In particular
P4 (`<prefix>_tarean/*.fasta_tarean/TRC_*.fasta`) does not match the 73 kept
`<prefix>_tarean/fasta/TRC_*.fasta` — the one collision the glob shapes make
plausible, and the reason P4 is anchored on `*.fasta_tarean/` rather than `*/`.

This is a static check of paths, not a substitute for the runtime tests in
§Testing: it proves the purge set does not *delete* a required file, not that
each consumer actually runs afterwards.

## Behaviour

**Flag.** `--cleanup` on `run_all`, `action="store_true"`, default `False`.
Absent = today's behaviour exactly.

**When it runs.** At the very end of the `run_all` branch, after `tarean()`
returns — i.e. after the report is built and after rDNA identification. Nothing
in the pipeline runs afterwards.

**Only on success.** If any step raised, the process has already exited non-zero
(as of 1.20.1 every whole-stage failure aborts) and cleanup never runs. It must
not be wired into a `finally:`. Rationale: a failed run's intermediates are
exactly what is needed to diagnose it, and a partially-complete directory has no
guarantees to preserve. This mirrors CARP, which cleans only after a successful
run, never on failure or dry-run.

**`--keep_rounds` is untouched.** `<prefix>_tidehunter_round{1,2,3}.gff3` are
debugging artefacts the user explicitly asked for; cleanup never deletes them,
and no purge pattern above matches them. Stated here because it is the one place
where "intermediate" and "requested output" collide.

**Idempotent and tolerant.** Every pattern is a glob; a pattern matching nothing
is not an error. Cleanup must never fail a run that otherwise succeeded — an
unlinkable file is a warning, not an abort.

**Reporting.** A summary line: how many files were removed and how much was
freed, plus the note that the run is no longer byte-complete. `--cleanup` is
recorded in `<prefix>_cmd_args.json` automatically (`save_args_to_file` persists
every non-`None` arg), so the run's own provenance says it was pruned — no new
side-car file, which would be at odds with the goal of fewer files.

*Resolved:* `cleanup_files_removed` / `cleanup_bytes_freed` are written into
`<prefix>_pipeline_stats.json`. They are deliberately **not** surfaced in the
report: the report is built before cleanup runs, so the card would only show them
after a subsequent re-render, and a number that appears on the second render but
not the first is worse than no number. The JSON is the durable record.

## Relationship to CARP (settled 2026-08-21)

CARP adopted this purge set in `e9e8494` as its own glob list (route (a) of the
companion FR), and **intends to keep maintaining it CARP-side even once
`--cleanup` ships** — their cleanup is one config key covering DANTE, DANTE_LTR,
DANTE_TIR, RepeatMasker and mmseqs scratch in the same pass, and carving the
TideCluster part out into a per-tool flag would make that contract harder to
explain, not easier.

Two consequences for this spec:

- **No double-purge risk.** CARP will not pass `--cleanup`, so the "please avoid
  applying both" warning in the FR is moot.
- **Two independent definitions of "TideCluster scratch" now exist and can
  drift.** Our tests should therefore assert the *capability contract* (G1–G4),
  not a file list: a layout change then fails loudly on both sides rather than
  silently widening one purge set. CARP's own test asserts the same property from
  their side ("maximal deleted a TideCluster capability file: …").

Their measurement at genome scale, for calibration: on a 94 Gbp assembly the same
set frees **90.2 %** of the three trees (40.23 GB of 44.59 GB), with
`monomers.RData` + `ggmin.RData` alone accounting for 62 % — versus 82.3 % and
14.6 % on the *Solanum* run used here. The trade improves with assembly size.

## Testing

`tests/unit.sh` gets one fixture-driven test that builds a miniature run
directory, applies the purge set, and asserts each guarantee:

- **G1** — reuse `tests/report_linkcheck.py` on the pruned copy. Note the
  existing `--purge` mode deletes a strict *superset* (whole `<prefix>_kite/`,
  `<prefix>_tarean/`, `dotplots/`), so G1 is already implied; the new check
  confirms it directly rather than by inference.
- **G2** — assert all six comparative inputs resolve, driven by the same list the
  README documents. This is the guarantee with no existing gate and is the reason
  a test was requested.
- **G3** — run `tc_rerender_report.py` on the pruned directory and require exit 0
  plus a clean `report_linkcheck` closure.
- **G4** — assert `tc_per_tra_consensus.py`'s four inputs resolve.

Plus a purge-set test proper: that the patterns match what they should on a
realistic tree and, specifically, that **`<prefix>_tarean/fasta/` is not matched**
by P4 — the one glob that could plausibly eat a kept path. CARP hit exactly this
trap while implementing the same set (the natural `*_tarean/*/TRC_*.fasta` matches
both copies) and guards it with a test; ours must too.

`tests/long.sh` gains a `run_all --cleanup` variant asserting the same four
guarantees end to end on a real (if small) run.

**Fixture note.** The bundled `long`/`short` fixtures produce no TAREAN output at
the default `-M` (`n_trcs_above_threshold: 0`), so the end-to-end gate passes
`-M 1000`, and `-l <library>` as well because `consensus_sequences_all.fasta` and
the RepeatMasker leftovers only exist on an annotated run.

*Resolved:* a trimmed two-sample comparative fixture is committed under
`tests/data/comparative` (576 KB), so `tests.sh determinism` runs instead of
skipping. Its README is explicit that it guards the determinism machinery but
does **not** reproduce issue #4's original non-determinism — the source dataset
is low-connectivity, so the prefilter has almost no cross-TRC edges to drop.

## Documentation to update alongside

The README **Output** section is the natural place to mark deliverable vs
intermediate, and it is currently stale independently of this work:

- `prefix_consensus_1` is documented but no longer produced (the code is
  commented out at `TideCluster.py:932`).
- `prefix_kite_report.html` is documented; nothing writes it (only
  `_move_v1_to_legacy` references it, and it skips missing files).
- The v1 HTML reports are documented at top level; since 1.17.0 they live in
  `<prefix>_report_legacy/`.
- Undocumented entirely: `<prefix>_kite/`, `<prefix>_report/`,
  `<prefix>_report_legacy/`, `dotplots/`, `<prefix>_pipeline_stats.json`,
  `<prefix>_cmd_args.json`, `<prefix>_seqid_lengths.tsv`, `<prefix>_rdna.tsv`,
  `<prefix>_per_tra_consensus/`, `<prefix>_tidehunter_short_annotation.*`.

Proposal: rewrite the section as a table with a **Kept by `--cleanup`** column,
so the purge set and the output documentation are one artefact that cannot drift
apart, and document `--cleanup` itself in Usage.

One naming issue to note while there: **`dotplots/` is not prefixed.** Two runs
sharing an output directory with different prefixes would collide in it. Out of
scope here, but it should be recorded somewhere.

## Open questions

1. Record `files_removed` / `bytes_freed` in `<prefix>_pipeline_stats.json` and
   show it in the report? (Proposed: yes.)
2. Commit the comparative fixture under `tests/data/` so G2 and the determinism
   test both stop depending on local data? (Proposed: yes, if the ~12 MB is
   acceptable in-repo.)
3. CARP alignment — route (a) or (b) of the companion FR. If (b), CARP passes
   `--cleanup` and drops its own `TideCluster_*` entries, which makes this spec
   the single definition of TideCluster scratch.
