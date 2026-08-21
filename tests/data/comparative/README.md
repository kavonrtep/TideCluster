# Comparative-analysis fixture

Two miniature TideCluster run directories, enough to run
`tc_comparative_analysis.R` end to end. Consumed automatically by
`tests.sh determinism` (which discovers this directory; override with
`$TC_COMPARATIVE_FIXTURE`).

## Contents

Each sample directory holds exactly the six files the comparative analysis reads
(see the README's "Files read from each TideCluster run"), and nothing else:

```
sample_{a,b}/
  tc_consensus_dimer_library.fasta            one consensus dimer per TRC
  tc_consensus/consensus_sequences_all.fasta  the per-array dimer pool
  tc_clustering.gff3                          array coordinates per TRC
  tc_annotation.gff3                          annotated regions
  tc_annotation.tsv                           annotation per TRC
  tc_tarean/SSRS_summary.csv                  which TRCs are SSRs
```

14 TRCs per sample: **10 shared** between the two, 4 unique to each. That yields
17 satellite families, 9 of them spanning both samples — so the fixture exercises
real cross-sample grouping rather than a degenerate one-TRC-per-family split.

## Provenance

Derived from a real *Solanum lycopersicum* TideCluster run by subsetting: at most
12 sequences per TRC, dropping sequences over 2500 bp (the multi-kb monomer tail
was almost all of the bulk and adds nothing here), and keeping only TRCs with at
least 2 surviving sequences. 576 KB total, down from 12 MB for the untrimmed
two-sample form.

## What this fixture does and does not prove

**Does:** the comparative pipeline runs end to end, and its output is stable
across repeat runs and thread counts — that is, it guards the determinism
machinery added for issue #4 (fixed `--max-seqs`, deterministic dedup, sorted
graph edges, and `--deterministic`'s byte-identical `.m8`) against regression.
Before this fixture existed the test simply SKIPped, so none of that was covered.

**Does not:** reproduce issue #4's original non-determinism. Checked directly —
the pre-fix script (`f7f019c^`) produces byte-identical canonical output on this
fixture at 1 and 15 threads. The reason is structural, not a matter of fixture
size: issue #4's mechanism is MMseqs2's `--max-seqs` prefilter dropping redundant
hits **between different TRCs**, and this source dataset is low-connectivity — of
72 families across the full 73-TRC run, 58 contain exactly one TRC per sample and
only one groups more than four. There are almost no cross-TRC edges to drop.
Raising the redundancy *within* a TRC does not help either (verified with a
601-sequences-per-sample variant: same hash), because the canonical comparison is
at TRC granularity.

Reproducing issue #4 would need a dataset where many *distinct* TRCs are mutually
similar — a genuinely redundant satellite-monomer pool. If such a run becomes
available, point `$TC_COMPARATIVE_FIXTURE` at it and re-run the pre-fix probe
before assuming the guard is live.
