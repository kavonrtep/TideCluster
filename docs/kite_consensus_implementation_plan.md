# KITE per-TRA consensus — implementation plan

Companion to `docs/kite_consensus_review.md`. This plan covers the
prototype implementation of Option 2 (phase-aligned PFM / "Fourier"
consensus). The prototype is deliberately kept **out of the main
TideCluster pipeline** until it has been validated on real data and a
decision is made to ship it.

## 0. Decisions captured from the review

| Topic | Decision |
|---|---|
| Algorithm | Option 2 — phase-aligned PFM (Fourier-style folding) |
| Source draft to reuse | `tarean/consensus_estimation.md` |
| Peak driving the consensus | *m*₁ only. *m*₂ may be missing; *m*₂ consensus is not produced in this phase |
| Phase canonicalisation across arrays in a TRC | deferred. Useful but non-trivial when divergence between TRAs is high; revisit after validation |
| Sub-repeat / HOR diagnostic | computed, reported as a column in the output TSV, does not alter the consensus |
| Placement | standalone prototype script under `tarean/`, not called by `TideCluster.py` or `kite.R` |
| Promotion to pipeline | conditional on validation results; separate plan step |

## 1. Architecture — standalone prototype

The prototype is invoked manually by the developer; the main pipeline
stays unchanged. Two calling modes cover all testing needs:

- **Single-array mode** — one TRA FASTA file plus its *m*₁; emits one
  consensus record and one PFM. Useful for spot-checking a particular
  problematic array.
- **Batch mode** — reads an existing KITE top-3 TSV plus the TRC
  FASTA dir produced during a TideCluster run; iterates over all TRAs,
  writes one combined consensus FASTA and one combined diagnostics
  TSV.

Only batch mode is used for the Drapa validation run in §5.

### Files added (all under `tarean/consensus_prototype/`)

```
tarean/consensus_prototype/
├── consensus_pfm.R           # the algorithm + CLI
├── README.md                 # how to run; expected inputs/outputs
└── lib/
    └── windowing.R           # helpers (window selection, base tally)
```

Nothing in `tarean/kite.R`, `TideCluster.py`, or `tc_rerender_report.py`
changes in this phase.

### Language choice

R. Matches `kite.R`, reuses `Biostrings` / `Biostrings::readDNAStringSet`
which the environment already has, and keeps the prototype adjacent to
the algorithm it would eventually integrate with.

### Dependencies

All already in the `tidecluster_1.8.0` conda environment:
- `Biostrings`
- `optparse`
- `parallel`
- base R (`stats::fft`, `stats::smooth.spline`)

No new conda deps. If validation shows we want column-entropy-driven
promotion to MSA (a future step, not this plan), `MAFFT` is also
already available.

## 2. Algorithm (Option 2 in detail)

Inputs per TRA:
- `S` — DNA sequence (`Biostrings::DNAString`), length `L`.
- `m` — monomer length in bp (integer; from KITE's *m*₁).
- parameters:
  - `window_size` (default 50 000 bp)
  - `n_windows` (default 16)
  - `H` — number of harmonics to retain (default chosen to cover 90 %
    of spectral energy up to position `m/2`)
  - `pseudocount` (default 0.5)

Steps:

1. **Window selection.** If `L <= n_windows · window_size`, use the
   whole array. Otherwise pick `n_windows` starting positions uniformly
   spaced across `[1, L - window_size + 1]`. Seeded RNG for
   reproducibility of random alternatives.
2. **Base tally.** Allocate `counts[m × 4]`. For each window
   `[w_start, w_end)`, iterate position-by-position; for every base
   `b` at position `p` (1-based, absolute in the array), increment
   `counts[((p - 1) mod m) + 1, b_index]`. Ns and non-ACGT bases are
   dropped (do not contribute to any row).
3. **PFM.** Normalise per row: `pfm[i, b] = (counts[i, b] + pseudo) /
   (sum_b counts[i, b] + 4·pseudo)`. Rows sum to 1.
4. **Consensus.** `consensus[i] = argmax_b pfm[i, b]`, one of `A C G
   T`. This is the per-TRA consensus monomer, length exactly `m`.
5. **Quality metrics**, all emitted per TRA:
   - `consensus_coverage` — total bases counted / (window_size ·
     n_windows). Flags under-sampled arrays.
   - `mean_entropy` — mean column Shannon entropy of the PFM (base-2).
     Low = clean consensus; high = noisy / drifting array.
   - `entropy_per_column` — full vector, emitted as a separate
     `<prefix>_consensus_entropy.tsv` for plotting.
   - `harmonic_diagnostic` — discrete Fourier transform of the
     A-proportion signal across positions 1…m (same on C, G, T);
     ratios of harmonic amplitudes at `m/k` for `k = 2, 3, 4, 5`.
     Peaks at `m/k` indicate a sub-repeat structure (or HOR); reported
     as a compact string `k2=0.18,k3=0.04,k4=0.02,k5=0.01` so the user
     can eyeball HOR-dominant arrays. Does not alter the consensus.
6. **Output.** Consensus sequence, PFM, and diagnostics; formats in
   §4.

### What is not done in this phase

- No phase canonicalisation across arrays. Two arrays in the same TRC
  may produce consensuses that are cyclic rotations of each other;
  inter-TRA comparison is left to downstream tools.
- No gap column (the monomer length is fixed at `m`; any base at
  `(p - 1) mod m = i` contributes to row `i` regardless of local
  indel drift).
- No outlier-monomer filter. Every base contributes; statistical
  averaging is relied on.
- No MSA. Deferred to an Option 3 promotion if validation shows
  Option 2 insufficient for high-indel arrays.

## 3. CLI / interface

### Single-array mode

```
Rscript tarean/consensus_prototype/consensus_pfm.R \
  --mode single \
  --fasta <TRA.fasta>       \
  --m     <monomer_length>  \
  --out   <prefix>
```

Output files:
- `<prefix>_consensus.fasta`    — one-record FASTA.
- `<prefix>_pfm.tsv`            — `position A C G T` matrix.
- `<prefix>_diagnostics.tsv`    — one-row TSV with the quality metrics.

### Batch mode

```
Rscript tarean/consensus_prototype/consensus_pfm.R \
  --mode batch \
  --kite-tsv    <prefix>_kite/monomer_size_top3_estimats.csv \
  --tarean-dir  <prefix>_tarean/fasta/ \
  --out-dir     <out_dir>     \
  --cpu         <N>
```

Output files under `<out_dir>`:
- `per_tra_consensus.fasta` — one record per TRA, name
  `TRC_<k>__<seqid>__<start>-<end>`.
- `per_tra_diagnostics.tsv` — one row per TRA with all columns from
  §2 step 5 plus `TRC_ID, seqid, start, end, length, m`.
- `per_tra_pfm/TRC_<k>__<seqid>__<start>-<end>.tsv` — one PFM per TRA.
  Optional; controlled by `--write-pfm` (default off — writing 2 000+
  PFM files is mostly useful for spot-checks, not for headline runs).
- `summary.log` — count of TRAs processed, failed, time, memory.

Parallelism is `parallel::mclapply` on the TRA list, matching
`kite.R`'s existing pattern.

### TRA sequence sourcing

`kite.R` already stores each TRA sequence (half of the dimer) inside
`<prefix>_tarean/fasta/TRC_<k>.fasta`. The prototype reads those same
files and matches records to the TSV rows by seqid/start/end encoded
in the FASTA record names (same convention as `kite.R`). This avoids
re-extracting from the assembly FASTA and keeps the prototype's
inputs fully local to an existing pipeline output directory.

## 4. Output format (unchanged from the review §6)

Single combined FASTA is the primary artifact. PFMs and diagnostics
are side-cars. No changes proposed to the GFF3 or to the existing
KITE TSV schema; the prototype writes to its own `--out-dir`.

## 5. Testing plan

The prototype needs convincing validation before any decision on
pipeline integration. Two fixtures:

### 5.1 Short smoke — synthetic inputs

Six synthetic TRA FASTAs generated by a helper script
`tarean/consensus_prototype/tests/make_synthetic.R`:

| fixture | monomer | mutation rate | indel rate | length | expected result |
|---|---|---|---|---|---|
| clean-short | random 50 bp | 0 % | 0 % | 50 kb | exact consensus |
| noisy-short | random 50 bp | 2 % | 0 % | 50 kb | consensus close to true, low entropy |
| indel-short | random 50 bp | 2 % | 0.5 % | 50 kb | consensus recognisable; entropy elevated |
| clean-large | random 500 bp | 0 % | 0 % | 2 Mb | exact consensus, fast |
| HOR-2x | random 171 bp × 2 | 1 % | 0 % | 500 kb | consensus ~342 bp when run with m = 342; k2 harmonic ≈ 0.5 |
| n-only | N × 10 kb | — | — | 10 kb | refuses, logs reason |

Pass criteria:
- Hamming distance from the true monomer ≤ `0.02 · m` on the clean
  and `≤ 0.10 · m` on the noisy-short fixture.
- Wall time ≤ 1 s on the 50 kb fixtures; ≤ 10 s on the 2 Mb fixture
  with 4 CPUs (batch mode over all six fixtures).
- HOR-2x harmonic diagnostic reports `k2 ≥ 0.3`.

### 5.2 Real-data validation — Drapa run2

Run batch mode on the existing Drapa output at
`tmp/drapa_run2/drapa_kite/monomer_size_top3_estimats.csv` +
`tmp/drapa_run2/drapa_tarean/fasta/`. That gives a 2 073-TRA batch
across 118 TAREAN-analysed TRCs plus the below-threshold TRCs
(which have KITE peaks but no TAREAN consensus to compare against).

Validation questions:

- **Quality vs TAREAN.** For TRCs where the TRC-level TAREAN consensus
  is available and the per-TRA KITE *m*₁ matches the TAREAN monomer
  length, compute the edit-distance between each per-TRA consensus
  and the TAREAN consensus. Distribution per TRC is the key
  diagnostic; if most arrays land within a few edits, the approach
  works.
- **Within-family variation.** Pick TRC_1 (161 arrays, 19 contigs).
  Plot the all-vs-all edit distance between per-TRA consensuses. A
  mostly-near-zero matrix with a few outliers validates the
  within-family variation story; a uniform distribution is a red
  flag.
- **HOR-dominant arrays.** For the 217 HOR-strong + 377 HOR-moderate
  arrays, check that the harmonic diagnostic reports `k2` or
  `k_n_max` above 0.15. If not, either the diagnostic is wrong or
  the HOR classifier was wrong.
- **Runtime.** Batch should complete in ≤ 5 min on 4 CPUs for 2 073
  arrays. Memory ≤ 4 GB.

Results are written to `tmp/consensus_prototype_drapa/` and
summarised in a short Markdown report
`tmp/consensus_prototype_drapa/REPORT.md` with the four diagnostics
above plus a handful of screenshots / inline consensus snippets.

### 5.3 Comparison fixture — CEN6

The 180 Mb CEN6 dataset at `test_data/CEN6_ver_220406.fasta` has
well-known centromeric satellites and a known clean monomer. Running
the prototype on that TideCluster output is the lightest-weight way
to confirm the consensus output matches published centromeric
satellites for that species. Lower priority than the Drapa run.

## 6. Success criteria for promoting to the main pipeline

The prototype graduates from `tarean/consensus_prototype/` into
`tarean/kite.R` only when **all** of the following hold:

1. All short-smoke fixtures pass their criteria (§5.1).
2. On Drapa run2:
   - median edit-distance between per-TRA consensus and TAREAN
     consensus (on matching monomer lengths) is ≤ `0.05 · m`;
   - the within-family edit distance for TRC_1 has a median ≤
     `0.10 · m` and 90th-percentile ≤ `0.25 · m` (satellites drift);
   - runtime ≤ 5 min, memory ≤ 4 GB.
3. At least three spot-checked HOR-dominant arrays have harmonic
   diagnostics consistent with their KITE `n` (e.g. n = 2 ⇒ k2
   dominates).
4. No systematic failures; failures are either transparent (the TRA
   is N-rich or *m*₁ is missing) or below 2 % of all TRAs.

If any criterion fails, write an addendum to this document detailing
the failure mode, and iterate on the algorithm (e.g. add phase
canonicalisation, promote to Option 4 iterative refinement, or fall
back to Option 3 MAFFT for the failing cases).

## 7. Promotion plan (only after §6 passes)

Ordered list, one small commit each. Each commit is reviewable on its
own; none of them ships until the one above it is green.

1. Move `tarean/consensus_prototype/consensus_pfm.R` to
   `tarean/consensus_pfm.R`; keep CLI unchanged for backward
   compatibility.
2. Call `consensus_pfm.R` in batch mode from `kite.R` after
   `best_peaks_concise` is built. Write outputs alongside the
   existing KITE artefacts.
3. Add five new columns to `monomer_size_top3_estimats.csv`:
   `consensus_length`, `consensus_mean_entropy`,
   `consensus_harmonic_k2`, `consensus_harmonic_k3`,
   `consensus_coverage`.
4. Ship the combined `per_tra_consensus.fasta` under
   `<prefix>_kite/` from `kite.R`.
5. In `tc_rerender_report.py`:
   - Load the new columns (backwards-compatible; legacy TSVs return
     empty values).
   - Add a `Consensus` preview column to the per-TRC dashboard arrays
     table (first 60 bp + ellipsis, with a copy button like the
     TAREAN consensus).
   - Add a "Download per-TRA consensus FASTA" button on each
     dashboard pointing at `<prefix>_kite/per_tra_consensus.fasta`.
6. Test tier `tests/consensus.sh` exercising the short fixtures and
   asserting the consensus FASTA is emitted with the expected
   headers.
7. Version bump to `1.10.0`; changelog entry noting the new
   per-TRA consensus column.

Steps 1–4 can ship without step 5 if we want the data before the UI;
steps 6–7 always last.

## 8. Non-goals / explicit deferrals

- **Phase canonicalisation.** Deferred (user decision: "nice to have
  but careful"). If required later, the cleanest hook is to align
  each per-TRA consensus to the TRC-level TAREAN consensus via a
  short rotation search, and rotate to the minimum-edit-distance
  frame before emitting. Straightforward but requires the TAREAN
  consensus to exist for that TRC.
- **MSA / Option 3.** Deferred. Reserve `tarean/consensus_msa.R` (the
  existing unstaged draft) for a second iteration if Option 2 turns
  out to be insufficient for high-indel arrays.
- **Integration into `tc_merge_annotations.py` / `tc_reannotate.py`.**
  These tools do not currently consume per-TRA consensus; not in
  scope.
- **Comparative analysis extension.** The comparative pipeline
  (`tc_comparative_analysis.R`) does not consume per-TRA consensus.
  Extension is a separate feature.

## 9. Risks

| risk | severity | mitigation |
|---|---|---|
| Fourier / PFM approach misidentifies consensus when monomer boundaries drift | medium | harmonic diagnostic surfaces drift; entropy column flags high-noise arrays; user can drop back to Option 3 for flagged TRAs |
| 2 Mb+ arrays are slow even with windowing | low | `window_size` and `n_windows` are tunable; default 50 kb × 16 = 800 kb is sampled not processed whole |
| Sharp memory spike with many TRAs in parallel | low | batch uses `mclapply`; each worker sees one TRA; PFM is `m × 4` integer matrix |
| Per-TRA consensus is too large to embed in report JSON for Drapa-sized runs | low | consensus stays in the combined FASTA, report stores only the 60 bp preview |
| Prototype becomes a permanent fork — rotten code that never merges | medium | the §6 success criteria + §7 promotion plan are explicit and short; if we skip §6 we kill the prototype, we do not let it linger |

## 10. Effort estimate

| piece | LOC | complexity |
|---|---|---|
| `consensus_pfm.R` core (windowing + tally + PFM + consensus) | ~150 | low |
| Harmonic diagnostic | ~40 | low |
| CLI (single + batch modes) | ~80 | low |
| `lib/windowing.R` helpers | ~60 | low |
| `tests/make_synthetic.R` | ~100 | low |
| `tarean/consensus_prototype/README.md` | ~60 | trivial |
| **total (prototype)** | **~490** | **~half a day** |

Validation run + report on Drapa: ~half a day including iteration if
short-smoke fails on first try.

Total to decision on promotion: one day of active work.
