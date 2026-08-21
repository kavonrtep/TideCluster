# Reply to issue #3 — "Align cleanup_intermediates: maximal with TideCluster's purge contract"

> Draft reply to `kavonrtep/assembly_repeat_annotation_pipeline#3`. Originally
> written 2026-08-21 against TideCluster **1.20.1** and CARP `e9e8494`; **revised
> the same day against TideCluster 1.21.1**, which obsoletes the whole
> comparative-analysis correction below. **Not yet posted** — paste the body as a
> comment from the host.
>
> The authoritative copy of this draft belongs in the CARP repo; this is
> TideCluster's record of the exchange.

---

Adopted — thanks, this was a good catch and the analysis held up everywhere we
could check it.

Shipped as option **(a)**, the CARP-side glob set, in `e9e8494`. `maximal` now
deletes the disposable files inside `TideCluster_{tarean,kite,consensus}/`
instead of the three trees.

## Verification on our side

Measured on a live 94 Gbp run (*drVisAlbu1.1*, 5,096 sequences), dry-running the
new rule over the real tree rather than a fixture:

| | |
|---|---|
| three trees | **44.59 GB** (tarean 39 GB, kite 3.5 GB, consensus 134 MB) |
| freed by the file-level set | **40.23 GB across 3,567 items = 90.2 %** |
| kept | 4.36 GB, and all three capabilities |

Breakdown, with each claim checked:

| Item | Size | Check |
|---|---|---|
| `monomers.RData` | 25.02 GB | `save()`d at `tarean/methods.R:654-655`; **no `load()` of it anywhere** in TideCluster |
| `ggmin.RData` | 2.76 GB | same |
| `*.kmers` | 8.43 GB | written by `kmer_counting.py`, never read back |
| `kitehor.periodogram` | 3.46 GB | input to `kite_heatmaps.R`, already rendered |
| per-array `TRC_*.fasta` | 0.55 GB | byte-identical to `tarean/fasta/TRC_*.fasta` (confirmed with `cmp`) |

Your tomato run gave 82.3 %; at genome scale it is **90 %**, because the
write-only `.RData` pair alone is 62 % of the trees. The larger the assembly, the
better the trade — 4.4 GB is a rounding error next to a 95 GB genome.

## One warning if you implement `--cleanup` with globs

TideCluster stores each array's FASTA **twice**:

```
TideCluster_tarean/TRC_1.fasta_tarean/TRC_1.fasta   <- disposable copy
TideCluster_tarean/fasta/TRC_1.fasta                <- needed by tc_per_tra_consensus.py
```

The natural pattern `TideCluster_tarean/*/TRC_*.fasta` matches **both** — it
would delete the file the whole change exists to preserve. We had to write
`*_tarean/TRC_*.fasta`, and our test asserts it (it fails with
`maximal deleted a TideCluster capability file: .../fasta/TRC_1.fasta` if the
glob is loosened). Worth a test on your side too.

*(Since shipped: TideCluster 1.21.0 added `run_all --cleanup` with the same
anchoring, and `tests/test_cleanup.py` asserts that a loosened
`_tarean/*/TRC_*.fasta` aborts **before** deleting anything. Same trap, same
guard, independently arrived at.)*

## Your open question 2 — yes, there is one

`make_repeat_report.R:646` and `make_summary_plots.R:192` both read

```
TideCluster/<run>/TideCluster_kite/monomer_size_top3_estimats.csv
```

to label a tandem family with its monomer length (`TRC_7 (172 bp)`). Under the
old tree-level purge this degraded silently: cleanup runs after the report, so a
single run was fine, but re-rendering a report from an archived run quietly
dropped every `(bp)` label. Your proposed set already keeps this file, so (a)
fixed it as a side effect — please keep it in the `--cleanup` keep set too.

Nothing else: no CARP manifest output lives inside the three trees, and neither
`TideCluster_clustering.gff3_1.gff3` nor the kite `_*longext*` files are
referenced anywhere in CARP.

## Argument 1 (comparative analysis) — correct at 1.20.1, obsolete at 1.21.1

**This section replaces the "does not apply to CARP" correction in the first
draft.** That correction was right about 1.20.1 and both of its reasons have
since been fixed, so it should not go out as written.

At the time of writing, comparative analysis could not consume a CARP tree for
two reasons, neither of them cleanup:

**a) Three of the six inputs were never written on a CARP run.**
`TideCluster_annotation.gff3`, `TideCluster_annotation.tsv` and
`TideCluster_consensus/consensus_sequences_all.fasta` all came from the
annotation step, which `run_all` skips unless a library is supplied, and CARP
passes `-l` only when the optional `tandem_repeat_library` config key is set.

**b) `get_seq_files()` hardcoded a `tc_` prefix**, ignoring the
`tidecluster_prefix` column while the other four paths honoured it, so those two
`readDNAStringSet` calls failed on any `TideCluster_*` directory.

Both are now fixed:

- **1.20.2** — `get_seq_files()` resolves both consensus FASTAs from the row's
  `tidecluster_prefix`, so `TideCluster_*` trees are readable. (b) is gone.
- **1.21.1** — `consensus_sequences_all.fasta` is written by the **clustering**
  step rather than annotation. It was only ever the per-TRC
  `TRC_*_dimers.fasta` concatenated, and those come from clustering, so the
  concatenation simply lived in the wrong step. Runs made by earlier versions
  are covered too: when the file is absent the comparative analysis rebuilds the
  pool from `TRC_*_dimers.fasta` and logs that it did, so **archived CARP runs
  work without being re-run**.

So of the six inputs, a CARP run made *without* `tandem_repeat_library` is now
missing only `TideCluster_annotation.{gff3,tsv}` — and both are optional: the
GFF3 falls back to `TideCluster_clustering.gff3`, and without the TSV the
`_annot` columns are simply empty. **Comparative analysis runs on CARP output
regardless of that config key.**

(1.21.1 also fixed an unrelated crash this surfaced: `cluster_ssrs_sequences()`
died with `replacement has 1 row, data has 0` when *no* sample contained an SSR
TRC. That would have hit any SSR-free pair of samples, CARP or not.)

## Consequence: your file-level set is now load-bearing for this too

When `e9e8494` was written, argument 1 was the one reason of the three that did
not apply to CARP. It does now, which makes two files inside the purged trees
newly important:

| File | Inside | Status in the file-level set |
|---|---|---|
| `TideCluster_consensus/consensus_sequences_all.fasta` | `TideCluster_consensus/` | kept — only `*_renamed.fasta*` is deleted there |
| `TideCluster_tarean/SSRS_summary.csv` | `TideCluster_tarean/` | kept — only `*.kmers`, the two `.RData` and the duplicate `TRC_*.fasta` go |
| `TideCluster_consensus/TRC_*_dimers.fasta` | `TideCluster_consensus/` | kept — needed by the rebuild fallback for pre-1.21.1 runs |

All three already survive `e9e8494`, so nothing needs changing. Worth recording
that they are now capability inputs rather than incidental survivors, so a future
narrowing of the set does not quietly drop them. Under the **old tree-level**
`maximal`, the first two would both have gone — so this is one more reason the
change was the right one.

One small disk note in the other direction: from 1.21.1 a run made without a
library writes `consensus_sequences_all.fasta` where it previously did not, and
its content duplicates the `TRC_*_dimers.fasta` that remain on disk. Measured on
a tomato run the duplicate is 18 % of `TideCluster_consensus/`, so on your 94 Gbp
run expect roughly **+24 MB** against 44.59 GB of trees — noise, but not zero,
and it cannot be purged because it is a comparative input.

## Your open question 3 — no

Dropping `tarean/fasta/` as well would add **0.61 GB of 44.59 GB** on our run.
Not worth giving up per-array consensus for 1.4 %.

## On (a) vs (b)

We went with (a) because `--cleanup` did not exist in 1.20.1. *(It shipped in
1.21.0.)* To set expectations: we expect to **keep doing this CARP-side** even
now that the flag exists. CARP's cleanup is one post-run step the user controls
through a single config key (`cleanup_intermediates: minimal|maximal|none`, with
`--keep-all` overriding it), and it covers DANTE / DANTE_LTR / DANTE_TIR /
RepeatMasker / mmseqs scratch in the same pass; splitting the TideCluster part
off into a flag on one tool would make the contract harder to explain, not
easier. Your point about our set being a list of paths we do not own is fair —
the mitigation is that the globs are narrow, documented with a reason each, and
covered by a test that fails loudly if the layout moves under us.

That also settles your "please avoid" note: we will not be passing `--cleanup`,
so there is no risk of both purges being applied.
