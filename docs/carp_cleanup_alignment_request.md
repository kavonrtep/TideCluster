# FR (CARP) — align `cleanup_intermediates: maximal` with TideCluster's purge contract

**To:** CARP maintainers · **From:** TideCluster · **Date:** 2026-08-21
**TideCluster version this is written against:** 1.20.1
**Related:** `docs/V6_carp_113_cleanup.md` (the 1.1.3 / 1.1.4 cleanup validation and FR-1/FR-2/FR-3)

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
