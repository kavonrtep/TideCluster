# KITE per-TRA consensus extraction — design review

Status: design review only, no code changes. Purpose is to evaluate
whether and how to extract a consensus sequence **per TRA** (tandem
repeat array), on top of the monomer-size estimate that KITE already
produces, and to pick between the implementation options.

## 1. Motivation

TAREAN already reconstructs one consensus per TRC (family). KITE
already produces one monomer-size estimate per TRA. What KITE does
**not** currently produce is a per-TRA consensus sequence at the
estimated monomer length.

Surfacing per-TRA consensus would:

- Expose **within-family variation**: two arrays in the same TRC can
  share a monomer length but carry distinct sequence variants — for
  example, arrays at different chromosomal locations after paralogous
  amplification, or arrays with locally drifting monomers.
- Make the HOR signal interpretable: an HOR-dominant TRA with
  `m₁ = 2·m₂` is much more informative when the per-TRA consensus at
  `m₁` is available for inspection alongside the consensus at `m₂`.
- Provide consistent per-TRA input for downstream tools (aligners,
  phylogenetics, custom libraries) that currently have to rebuild
  consensuses from the GFF3 array by themselves.

## 2. Current KITE pipeline

Recap from `tarean/kite.R`:

1. Every *k*-mer in the TRA is extracted; for each distinct *k*-mer
   the vector of gaps between consecutive occurrences is recorded
   (`calculate_neighbor_distances`, *k* = 6 by default).
2. The aggregate distribution of gaps (over all *k*-mers) is a
   profile indexed by integer bp distance.
3. A null profile is built from shuffled sequences with matching
   nucleotide composition (`get_neighgor_distances_background`).
4. Peaks significantly above the null are extracted
   (`simplified_findpeaks`, `get_peaks_from_neighbor_distances`).
5. Each peak position corresponds to a candidate monomer period.
   The top three by weighted score become *m*₁, *m*₂, *m*₃ with
   scores *s*₁, *s*₂, *s*₃.
6. HOR classification consumes these six numbers (see
   `docs/hor_classification.md`).

Notice that the **positions** of each *k*-mer are present in the
pipeline (transiently, inside `calculate_neighbor_distances`) but
discarded after the distance distribution is built. Consensus
reconstruction needs exactly those positions.

## 3. What "consensus on the same principle" means

Given a TRA and its estimated monomer length *m*, the consensus at
period *m* can be derived from the same *k*-mer / position
information that KITE already produces. The core idea:

- For any *k*-mer that occurs multiple times in the TRA, its
  positions modulo *m* cluster tightly if *m* is the true period —
  the *k*-mer has a consistent **phase** within the monomer.
- Averaging base frequencies across corresponding positions (taken
  mod *m*) over the whole TRA yields a position-frequency matrix
  (PFM) at resolution *m*. The argmax per column is the consensus
  monomer.

The non-trivial parts are (a) picking a phase/frame, (b) handling
monomer variability (insertions / deletions / drift), and (c)
scaling to arrays with millions of bp.

## 4. Options

Five options are laid out below. Two of them — Fourier/spectral and
MSA-based — already exist as unstaged drafts in this repository; the
review treats them on the same footing as the others so the choice
can be made head-to-head.

### Option 1 — Phase-aligned stacking with a single anchor

Use the *k*-mer with the highest periodicity score at period *m* as
the anchor. Its positions define frame 0. Cut the TRA into *m*-bp
chunks starting at each anchor position; stack-align them; compute
the consensus.

Pros: simplest; minimal code; uses only what KITE already has.
Cons: a single anchor propagates its own noise; fails when the
anchor shifts mid-array due to local insertion / deletion. No indel
handling within the monomer. Quality of consensus degrades roughly
linearly with monomer drift.

### Option 2 — Phase-aligned PFM with multiple anchors (Fourier approximation)

Accumulate a position × base count matrix of size *m × 4*:
- Scan the array; for every base at position *p*, increment
  `counts[(p - φ) mod m][base]`, where *φ* is a global phase
  computed from the consensus of many anchor *k*-mers' phases.
- PFM from counts; consensus from argmax per column.

This is what `tarean/consensus_estimation.md` describes as the
"Fourier / harmonic (spectral) consensus" — folding the array by *m*
and optionally smoothing the profile in frequency space to suppress
subrepeats. Relies on *m* being correct and the array being close to
stationary over the chosen windows.

Pros: O(L + m log m) per array; handles multi-Mb arrays comfortably;
robust to small monomer-level noise via statistical averaging; no
external tool; produces PFM + consensus + optional sub-repeat
diagnostics.
Cons: loses local phase when boundaries drift; a single point
insertion shifts all downstream bases into the wrong bin; quality
degrades on arrays with frequent indels.

### Option 3 — MSA of candidate monomers (MAFFT-based)

Identify candidate monomer start positions using high-frequency,
phase-consistent *k*-mers as anchors; extract overlapping candidate
monomers around each anchor; outlier-filter (for instance via Mash
distances); run MAFFT on the surviving set; derive consensus + PSSM
from the alignment.

This is the approach drafted in
`tarean/consensus_msa_implementation_plan.md` and partly implemented
in `tarean/consensus_msa.R`.

Pros: explicitly handles indels and boundary drift; MAFFT is the
standard tool for this; alignment is itself a useful downstream
artifact (usable for phylogenetic analysis, HOR detection, custom
libraries). MAFFT is already in the conda environment.
Cons: slow on large arrays (alignment is quadratic in the worst
case; hundreds of monomers per array); needs either sampling or a
bounded-size alignment; more moving parts; dependency on
`minimap2` / `Mash` if outlier filtering is enabled. Option 3 as
drafted adds Mash to the environment.

### Option 4 — Iterative consensus refinement

Start with any seed consensus (for example, Option 2's output).
Align every candidate *m*-bp window from the array to the seed
(banded pairwise, or minimap2); update the consensus column-by-column
from the alignment; iterate to convergence.

Pros: combines Option 2's speed with Option 3's indel handling;
converges to high-quality consensus even from noisy seeds; natural
way to produce a PSSM.
Cons: requires an alignment library; more code; convergence is not
guaranteed if the seed is far off; adds minimap2 dependency if used
as the aligner.

### Option 5 — Mini-TAREAN per TRA (de Bruijn reuse)

Re-use TAREAN's *k*-mer de Bruijn graph approach at the TRA level.
Build the graph, search for the best cycle of length ≈ *m*, emit the
cycle's sequence as consensus.

Pros: consistent with TAREAN's TRC-level approach; uses infrastructure
already in the codebase; works even when the periodicity signal is
complex.
Cons: TAREAN is designed for large read collections; on a single
array, *k*-mer coverage may be uneven and the graph may branch
ambiguously; may fail on exactly the arrays where consensus is most
useful (short or noisy ones). Would require non-trivial adaptation
of `tarean.R` and `methods.R`.

## 5. Comparison

Quality is subjective; the table below is my estimate based on the
algorithmic characteristics, not on measured benchmarks.

| Option | Quality on clean arrays | Indel tolerance | Speed (Mb arrays) | External tools | New code | Integration burden |
|---|---|---|---|---|---|---|
| 1 Phase-stacking (single anchor) | medium | none | fast | — | small | low |
| 2 Phase-PFM (Fourier) | high | low | fast | — | medium | low |
| 3 MSA (MAFFT) | highest | high | slow without sampling | MAFFT (have), Mash (new) | larger | medium |
| 4 Iterative refinement | high | high | medium | minimap2 (new) | larger | medium |
| 5 Mini-TAREAN | variable | low | medium | — (uses TAREAN) | TAREAN adaptation | high |

## 6. Output format

Independent of which option is chosen, the same output contract
works:

- **Per-TRA consensus FASTA**: one sequence per TRA, in
  `<prefix>_kite/per_tra_consensus.fasta`. Sequence name:
  `TRC_<k>__<seqid>__<start>-<end>` so arrays can be matched back to
  the KITE TSV and the GFF3.
- **PSSM side-car** (optional): one TSV per TRA under
  `<prefix>_kite/pssm/TRC_<k>__<seqid>__<start>-<end>.tsv`, matrix
  columns `A C G T` rows per position. Useful for drill-down.
- **GFF3 attribute** (optional, probably a bad idea for long
  monomers): embedding a 500-bp consensus inside a GFF3 attribute is
  feasible but makes the file hard to parse; the FASTA is the better
  primary output.
- **Report v2 surface**:
  - Per-TRC dashboard arrays table gains a `Consensus` column that
    shows a short preview (first 60 bp + ellipsis) with a
    copy-to-clipboard button, like the TAREAN consensus already has.
  - Dashboard gains a `Download per-TRA consensus FASTA for this
    family` button.

## 7. Integration points

- `tarean/kite.R`
  - Keep the KITE pipeline structure. After `get_peaks_from_seq`
    returns the top peaks for an array, call the consensus helper
    with the top peak position.
  - Persist per-TRA consensus into `best_peaks_concise` as a new
    column. Write alongside the existing TSV output.
- `TideCluster.py` `tarean()`
  - No change; kite.R handles consensus in-process.
- `tc_rerender_report.py`
  - Load the new column (backwards-compatible; legacy KITE TSVs
    simply return empty consensus).
  - Render preview column + download button.

## 8. Recommendation

**Option 2 (phase-PFM / Fourier) as the initial implementation.** It
is the fastest, has no new external dependencies, and degrades
gracefully on moderately noisy arrays. The existing
`tarean/consensus_estimation.md` already lays out the algorithm in
enough detail to implement directly.

**Promote to Option 3 (MSA) selectively**, for arrays where:
- the monomer is short enough to make MAFFT tractable (say ≤ 2 000 bp
  — above that, full MSA of all monomers in a multi-Mb array becomes
  slow), **and**
- the user asks for it via a flag (`--consensus-method mafft`), or
  KITE detects that Option 2's PFM has high column entropy
  (indicating boundary drift).

Option 4 (iterative refinement) is a superset of 2 + 3 and an
interesting second phase once Option 2 is in place; it can consume
Option 2's PFM as its seed. Option 1 is included for completeness
but would deliver lower-quality consensus than 2 at no cost saving.
Option 5 is plausible but introduces TAREAN's complexity to a
context where a lighter-weight method suffices.

## 9. Open questions before implementation

1. **Which peak drives the consensus?** *m*₁ alone, or *m*₁ + *m*₂
   (two consensuses per TRA for HOR-dominant arrays to let users
   compare the monomer and the HOR period)? Recommend: *m*₁ alone as
   default, second consensus behind a flag for HOR investigations.
COMMENT - m1
2. **Phase convention.** Arrays with the same family can start at
   different positions (different rotations of the same cyclic
   consensus). Should the tool canonicalise rotation across arrays
   in the same TRC (e.g. align each per-TRA consensus to the TRC-
   level TAREAN consensus and rotate to match)? Recommend yes — it
   is cheap and makes downstream comparison sane.
COMMENT - phasing step is very usefull,nice to have, but we need to carefully plan how to do then it can be complicated when divergence between TRA is high
3. **Quality gate.** Should the tool emit a consensus for arrays
   where KITE returned no confident *m*₁? Current suggestion: no;
   emit the array's seqid and an explanation in the TSV so it stays
   visible.
COMMENT - not, but I think all TRA have at least m1, m2 could be missing
4. **Subrepeat / HOR handling.** When *m*₁ / *m*₂ integer multiple,
   the *m*₁-cycle consensus may show a banded PFM where position
   *i* and *i + m*₂ look like shifted versions of each other.
   Option 2 can surface this via a harmonics diagnostic.
   Recommendation: compute the diagnostic, report it as a column in
   the KITE TSV, do not alter the consensus itself.
COMMENT - agreed
5. **Existing drafts**: `tarean/consensus_msa.R` and
   `tarean/consensus_msa_implementation_plan.md` (Option 3) plus
   `tarean/consensus_estimation.md` and
   `tarean/consensus_implementation_plan.md` (Option 2 / 3 mix) are
   not under version control yet. Before implementation, decide
   whether to adopt one of them as the starting point or reset to
   this review and build fresh. Recommendation: start from
   `consensus_estimation.md` for Option 2 (the algorithm there is
   clean and close to what is needed); reserve `consensus_msa*.R`
   for the Option 3 promotion.
COMMENT - DECISION  we will go wit Option 2 and yes reuse the document
## 10. Out of scope

- Consensus sequences at the TRC level already come from TAREAN.
  This review does not propose changing TAREAN.
- Phylogenetic analysis of per-TRA consensuses across a family or
  across samples is a potential downstream use of Option 3's
  alignment, but it is a separate feature.
- Read-based polishing (racon / medaka-style) is not considered; the
  inputs here are assembly sequences, not reads.
