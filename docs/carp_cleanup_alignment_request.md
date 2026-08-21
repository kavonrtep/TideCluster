# FR (CARP) — align `cleanup_intermediates: maximal` with TideCluster's purge contract

**To:** CARP maintainers · **From:** TideCluster · **Date:** 2026-08-21
**TideCluster version this is written against:** 1.20.1
**Related:** `docs/V6_carp_113_cleanup.md` (the 1.1.3 / 1.1.4 cleanup validation and FR-1/FR-2/FR-3)

## Outcome — ADOPTED (CARP `e9e8494`, 2026-08-21)

CARP took route **(a)**: `maximal` now deletes the disposable files *inside*
`TideCluster_{tarean,kite,consensus}/` instead of the three trees. Their
measurements on a live 94 Gbp run (*drVisAlbu1.1*, 5,096 sequences) came out
better than ours at tomato scale:

| | this document (*Solanum*, 1.2 GB run) | CARP (*drVisAlbu1.1*, 94 Gbp) |
|---|---|---|
| freed | 82.3 % of the run | **90.2 %** of the three trees (40.23 GB of 44.59 GB, 3,567 items) |
| dominant item | `kitehor.periodogram` (38.7 %) | `monomers.RData` + `ggmin.RData` (**62 %** of the trees, 27.8 GB) |

The direction is consistent with what we measure: `monomers.RData` already
outweighs `ggmin.RData` 4:1 on the *Solanum* run (135 MB vs 32 MB); at genome
scale the ratio is ~9:1 and the pair dominates everything else. The bigger the
assembly, the better this trade — which is the opposite of what one might fear.

**They independently hit the glob trap we flagged.** `TideCluster_tarean/*/TRC_*.fasta`
matches both the disposable per-TRC copy and the `fasta/` original that the whole
change exists to preserve; they anchor on `*_tarean/TRC_*.fasta` and assert it in
a test. Our P4 is anchored the same way and our static validation confirms it —
see `docs/output_cleanup_spec.md`.

### Their answers to the open questions

**Q2 — yes, one CARP dependency lives inside the trees.** `make_repeat_report.R:646`
and `make_summary_plots.R:192` read
`TideCluster_kite/monomer_size_top3_estimats.csv` to label tandem families with
their monomer length (`TRC_7 (172 bp)`). Under the old tree-level purge this
degraded silently — cleanup runs after the report, so a single run looked fine,
but re-rendering from an archived run dropped every `(bp)` label. The file is in
our keep set already (G3/G4 both need it); it is now recorded there as a CARP
dependency too, so it cannot be dropped later by accident.

**Q3 — no.** Dropping `tarean/fasta/` as well would add 0.61 GB of 44.59 GB
(1.4 %) on their run. Not worth losing per-array consensus for. That settles the
"biggest remaining lever" note in the spec: the lever is real at tomato scale
(39 % of the residual) but negligible at genome scale, so per-TRA consensus stays
in the contract.

**(a) vs (b) — they will keep doing it CARP-side even after `--cleanup` ships.**
Their reasoning is sound: CARP's cleanup is one post-run step the user controls
through one config key, covering DANTE / DANTE_LTR / DANTE_TIR / RepeatMasker /
mmseqs scratch in the same pass, and splitting the TideCluster part into a flag on
one tool would make the contract harder to explain. The mitigation for "a list of
paths CARP does not own" is that the globs are narrow, individually justified, and
covered by a test that fails loudly if our layout moves.

*Consequence for us:* two independent definitions of "TideCluster scratch" will
exist and can drift apart. Our test should therefore assert the capability
contract rather than a file list, so a layout change fails on both sides. It also
settles the "please avoid applying both" warning below — CARP will not pass
`--cleanup`, so there is no double-purge risk.

### Their correction on argument 1 (comparative analysis) — valid, and now stale

*(The revised reply, updated for 1.20.2 + 1.21.1, is kept at
`docs/carp_issue3_reply_draft.md`.)*

They are right that argument 1 did not apply to CARP, for two reasons:

**a) Three of the six comparative inputs are never written on a CARP run.**
`TideCluster_annotation.gff3`, `TideCluster_annotation.tsv` and
`TideCluster_consensus/consensus_sequences_all.fasta` all come from the annotation
step, which `run_all` runs only `if cmd_args.library:`. CARP passes `-l` only when
the optional `tandem_repeat_library` config key is set. Confirmed on our side —
this is the same finding that made us document "comparative analysis needs
`run_all` **with** a library" in the README.

**b) `get_seq_files()` hardcoded a `tc_` prefix**, so the two `readDNAStringSet`
calls failed on any `TideCluster_*`-named directory regardless of cleanup. Correct
for 1.20.1.

**(b) is fixed** — commit `354af1f`, unreleased at the time of their reply. Both
consensus FASTA paths now resolve from the row's `tidecluster_prefix` like the
other four. Verified against a CARP-shaped fixture (`input_dir =
TideCluster/run-000170`, `tidecluster_prefix = TideCluster`): 1.20.1 dies with
`cannot open file .../run-000170/tc_consensus_dimer_library.fasta`; with the fix
the run completes and writes both samples' `gff3/` exports. Every other `tc_`
literal left in the script is a default argument value that the real call path
always overrides.

So after that release, comparative analysis **will** run on CARP output — and as
of **1.21.1** their (a) is gone too: `consensus_sequences_all.fasta` is written by
the clustering step rather than annotation, and pre-1.21.1 runs get it rebuilt
from `TRC_*_dimers.fasta`. A CARP run without `tandem_repeat_library` is then
missing only the two annotation *reports*, both optional. Comparative analysis
runs on CARP output regardless of that config key.

That makes two files inside the purged trees newly load-bearing for CARP —
`TideCluster_consensus/consensus_sequences_all.fasta` and
`TideCluster_tarean/SSRS_summary.csv`, plus `TRC_*_dimers.fasta` for the rebuild
fallback. All three already survive `e9e8494`, but they are now capability
inputs rather than incidental survivors.

One small clarification for anyone reading their reply later: the ignored variable
was `prefix`, which carries the input table's `sample_code`, while lines 1361-1363
build their paths from `tc_code` (the `tidecluster_prefix` column). The diagnosis
is right — the paths ignored the TideCluster prefix — but the two columns are
different things, and `sample_code` is not what those paths should have used.

---

*Original request follows.*

## Summary

CARP's `cleanup_intermediates: maximal` deletes the TideCluster scratch trees
**whole**:

```
TideCluster_tarean/   TideCluster_kite/   TideCluster_consensus/
```

That frees the right amount of disk, but it removes more than scratch. Three
things a TideCluster run directory can still do after a `minimal` cleanup stop
working after a `maximal` one:

1. **It can no longer be used as comparative-analysis input.**
2. **Per-TRA consensus (`tc_per_tra_consensus.py`) can no longer be run.**
3. **The HTML report can no longer be re-rendered** (`tc_rerender_report.py`).

None of these are report-*viewing* problems — FR-3 fixed that, see below. They
are capability losses that are invisible until someone tries, months later, to
run a cross-genome comparison over a batch of archived runs.

TideCluster is adding its own `--cleanup` flag whose purge set is **file-level
inside those trees rather than the trees themselves**. It frees essentially the
same disk (**82.3 %** of a real run) while keeping all three capabilities. This
FR asks CARP to adopt that set, either directly or by delegation.

## What changed on our side since the 1.1.4 re-validation

`docs/V6_carp_113_cleanup.md` gates the flip to `maximal` (FR-2) on FR-3: the 28
dead links after a maximal purge — 9 TAREAN drill-down `report.html` links from
the served per-TRC pages, 19 in the legacy report tree.

**FR-3 shipped in TideCluster 1.18.0.** `<prefix>_report/` now vendors each
per-TRC TAREAN `report.html` plus its `img/` and `ppm_*.csv`, and the legacy
tree's image references were repointed at the vendored copies. `tests/long.sh`
gates it permanently with `tests/report_linkcheck.py --purge`, which deletes
`<prefix>_{kite,tarean}/` and `dotplots/` outright and asserts the link closure
is clean. So **report viewing after a maximal purge is now correct**, and the
FR-2 blocker as originally written is gone.

The three losses above are what remains, and they were never part of the
original FR-3 analysis because that spike only checked the *served report*.

## Detail: what `maximal` takes away

### 1. Comparative analysis (`tc_comparative_analysis.R`)

It reads exactly six files per run. Two are consensus FASTAs:

| File | Survives `maximal`? |
|---|---|
| `TideCluster_consensus_dimer_library.fasta` | yes — a sibling **file**, not inside the purged directory |
| `TideCluster_consensus/consensus_sequences_all.fasta` | **no** — inside `TideCluster_consensus/` |
| `TideCluster_clustering.gff3` | yes |
| `TideCluster_annotation.gff3` | yes |
| `TideCluster_annotation.tsv` | yes |
| `TideCluster_tarean/SSRS_summary.csv` | **no** — inside `TideCluster_tarean/` |

Missing `consensus_sequences_all.fasta` is not a graceful degradation: the R
script aborts at `readDNAStringSet` with `cannot open file …` and exit 1. The
missing SSRS table degrades quietly (no SSR grouping).

Note `V6` already observed that the dimer library survives "a sibling FILE, not
inside the purged `TideCluster_consensus/` dir" — correct, but it is only *one*
of the two required sequence pools, so surviving alone is not enough.

### 2. Per-TRA consensus (`tc_per_tra_consensus.py`)

Needs `TideCluster_tarean/fasta/`, `TideCluster_kite/monomer_size_top3_estimats.csv`,
`TideCluster_tidehunter.gff3` and `TideCluster_clustering.gff3`. `maximal` removes
the first two.

### 3. Report re-render (`tc_rerender_report.py`)

Needs the kite CSVs and `kitehor.rescored.peaks.tsv`, each per-TRC
`*.fasta_tarean/{report.html,img/,*.csv}`, and `dotplots/`. `maximal` removes all
but `dotplots/`. The already-rendered report keeps working; it simply cannot be
regenerated — which matters when a TideCluster upgrade improves the report and
the natural fix is to re-render archived runs rather than re-run them.

## Proposal

Replace the three whole-tree deletions with this file-level set. Measured on a
real run (*Solanum lycopersicum*, 1.2 GB, 2328 files):

| Purge | Size | % of run | Why it is safe |
|---|---|---|---|
| `TideCluster_kite/kitehor.periodogram` | 444 MB | 38.7 % | input to `kite_heatmaps.R` only; the report renders from `profile_plots/*.png` + `monomer_size_top3_estimats.csv` |
| `TideCluster_tarean/*/*.kmers` | 246 MB | 21.4 % | written by `kmer_counting.py`, never read back |
| `TideCluster_tarean/*/{ggmin,monomers}.RData` | 168 MB | 14.6 % | **write-only** — `save()`d at `tarean/methods.R:654-655`, never `load()`ed anywhere in the codebase |
| `TideCluster_tarean/*/TRC_*.fasta` | 76 MB | 6.6 % | byte-identical duplicate of `TideCluster_tarean/fasta/TRC_*.fasta`, which is kept |
| `TideCluster_consensus/*_renamed.fasta{,.cat,.masked}` | 12 MB | 1.0 % | RepeatMasker leftovers |
| `TideCluster_kite/_*longext*`, `TideCluster_clustering.gff3_1.gff3` | 0.5 MB | ~0 % | documented intermediates |
| **total** | **~946 MB** | **82.3 %** | run goes 1.2 GB → 204 MB |

For comparison, today's `maximal` on the same run deletes the three trees
(~1.1 GB) — about 13 percentage points more disk, at the cost of all three
capabilities above.

### Two ways to adopt it

**(a) Narrow `maximal` in `scripts/cleanup_outputs.py`.** Replace the three
`TideCluster_*` tree entries with the globs above. Self-contained in CARP; no
version coupling beyond requiring TideCluster ≥ 1.18.0 for the report half
(already required for FR-2 anyway).

**(b) Delegate to TideCluster.** TideCluster is adding `--cleanup` to `run_all`,
which applies exactly this set at the end of a successful run. CARP would pass
the flag and drop all `TideCluster_*` entries from its own `maximal` set. This
keeps one definition of "TideCluster scratch" in one place, and it stays correct
as TideCluster's file layout evolves — a real risk, since CARP's set is a list of
paths CARP does not own. Requires the TideCluster release carrying `--cleanup`.

We prefer **(b)**, with **(a)** as the interim if you want the saving before that
release lands. Either way the three capabilities are preserved and FR-2 can flip
the server default to `maximal`.

### Please avoid

Applying both — CARP running TideCluster with `--cleanup` *and* then applying its
own tree-level `maximal` — which restores exactly the losses this FR is about.

## Verification we can offer

The claims above are checked in TideCluster's own suite: `report_linkcheck.py
--purge` (report viewing after the trees are gone) and, with the `--cleanup`
work, one check per guarantee — comparative inputs resolve, per-TRA consensus
inputs resolve, and a re-render after cleanup produces a clean link closure.

If it helps, the `.carp-spike/check_cleanup_outputs.py` harness can be extended
with the same three capability checks, so the next spike run reports them
alongside the 56 server-dependency checks.

## Open questions for CARP

1. (a) or (b)?
2. Does anything in CARP's manifest or the server's `domain/carp_outputs.py` read
   from inside `TideCluster_{tarean,kite,consensus}/` beyond what V6 enumerated?
   If so those paths need adding to the keep set on either route.
3. Is there value in a CARP-side `cleanup_intermediates: maximal` that also drops
   `TideCluster_tarean/fasta/` (a further ~80 MB here) for deployments that will
   never run per-TRA consensus? TideCluster's own `--cleanup` keeps it.
