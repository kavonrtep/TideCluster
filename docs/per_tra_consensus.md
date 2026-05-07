# Per-TRA consensus generation

Companion R tools that produce one consensus monomer sequence per
Tandem Repeat Array (TRA), validate it against the source array by
self-BLAST, and label every TRA with a quality grade plus diagnostic
flags. Distinct from the per-TRC (per-cluster) consensus already
emitted by the TAREAN step.

Status: companion tooling under `tarean/consensus_prototype/`. Not yet
wired into `TideCluster.py`. Run separately on the outputs of an
existing TideCluster run.

## Vocabulary

| term | meaning |
|---|---|
| **Monomer** | One copy of the repeat unit. |
| **Tandem repeat array (TRA)** | One contiguous tandem-repeat region in the assembly, detected by TideHunter. Has genomic coordinates and an estimated monomer length `m`. |
| **Tandem repeat cluster (TRC)** | A group of TRAs whose monomers share sequence; assigned at the clustering step. Existing pipeline term. |
| **HOR (higher-order repeat)** | Periodic super-structure where the true repeat unit is `n × m` for some `n > 1`. KITE classifies arrays on a four-level scale (No HOR / weak / moderate / strong). |
| **Per-TRA consensus** | One consensus monomer for a single TRA. The output of this tool. |
| **Per-TRC consensus** | The cluster-level consensus already produced by TAREAN. Different artefact. |

### Methods

| method | description |
|---|---|
| **Array-MSA consensus** | Scans the genomic array, picks a k-mer that fires at the same offset modulo `m` (the *phase anchor*), extracts every in-phase monomer copy, aligns them with MAFFT, takes the column-wise majority. |
| **TideHunter consensus** | Uses TideHunter's own per-fragment consensus. For each TRA, the TideHunter fragment with the highest `copy_number × consensus_length` is selected and its `consensus_sequence` attribute is used verbatim. |
| **Selected consensus** | Output of the **consensus selector**: per TRA, whichever of array-MSA or TideHunter consensus self-BLASTs better against the source array. |

### Validation

The validator BLASTs each per-TRA consensus against the source array
(self-BLAST). A correct consensus tiles the array; the geometry of the
self-hits reveals consensus quality, array boundary correctness, and
internal array variability.

| term | meaning |
|---|---|
| **Coverage fraction** (`coverage_frac`) | Fraction of TRA bp covered by merged self-BLAST hits. Primary correctness metric. |
| **Core repeat region** | Span between the leftmost and rightmost self-BLAST hit. The portion of the TRA actually containing the repeat. |
| **Core coverage** (`core_coverage`) | Coverage measured against the core repeat region instead of the full TRA. Tolerant to TideHunter boundary mis-annotation. |
| **Core fraction** (`core_fraction`) | Core repeat region length / TRA length. Values < 1 indicate annotation flanks. |
| **Identity dispersion** (`pident_sd`) | SD of per-hit identity. Low = uniform array, high = heterogeneous (HOR, polymorphism). |
| **Perfect-copy fraction** (`perfect_hit_frac`) | Fraction of self-hits at `pident ≥ 99.5`. High = many monomer copies match consensus exactly; zero = the consensus averages divergent variants. |
| **Internal gap** | Bp range inside the core repeat region that no self-hit covers. Sub-repeat or insertion candidate. |

## Quality grades

A coverage tier per TRA. Independent of flag state.

| grade | rule | meaning |
|---|---|---|
| **A** — clean pass | `coverage_frac ≥ cov_strict` (default 0.90) | consensus tiles the array reliably |
| **B** — qualified pass | `cov_qualified ≤ coverage_frac < cov_strict` (default 0.80–0.90) | usable; usually accompanied by a flag explaining the deficit |
| **C** — low confidence | `cov_low ≤ coverage_frac < cov_qualified` (default 0.50–0.80) | consensus may still be informative but not for strict use |
| **D** — unusable | `coverage_frac < cov_low` (default < 0.50) or no self-hit | dominated by SSRs that should route to a motif-only path |

## Diagnostic flags

Independent boolean flags. Multiple may co-occur. A grade-A TRA can
carry flags (e.g. a HOR array tiles fully but is heterogeneous).

| flag | rule (default thresholds) | what it indicates |
|---|---|---|
| `boundary_overext` | `core_coverage ≥ 0.90` AND `core_fraction < 0.95` | TRA boundary contains flank; consensus is fine. Annotation issue, not a consensus issue. |
| `internal_gap` | `max_internal_gap_bp ≥ 1000` | sub-repeat or insertion candidate inside the array. |
| `heterogeneous` | `pident_sd ≥ 1.5` | array has divergent monomer variants. |
| `low_pident` | `pident_p05 < 92` | the worst 5 % of monomer copies are more than 8 % divergent from the consensus. |

## Output columns

`per_tra_metrics.tsv`, one row per TRA:

| column | meaning |
|---|---|
| `id` | `TRC__seqid__start-end` |
| `source` | `M_B` (array-MSA) or `TH_single` (TideHunter) — which consensus was selected |
| `array_length`, `cons_length` | TRA length and consensus length, bp |
| `hit_count` | self-BLAST HSP count |
| `coverage_frac` | primary coverage metric |
| `core_coverage`, `core_fraction`, `flank_bp` | core-region metrics |
| `mean_pident` | length-weighted mean identity |
| `pident_sd` | identity dispersion |
| `pident_p05`, `pident_p50`, `pident_p95` | identity quantiles |
| `n_perfect_hits`, `perfect_hit_frac` | perfect-copy counts |
| `max_internal_gap_bp`, `n_internal_gaps` | internal-gap metrics |
| `hit_ratio` | `hit_count / (array_length / cons_length)` (≈1 means every copy hit) |
| `cov_mb`, `cov_th` | coverage from each method (selector input) |
| `quality_grade` | A / B / C / D |
| `flag_boundary_overext`, `flag_internal_gap`, `flag_heterogeneous`, `flag_low_pident` | independent diagnostic flags |
| `flags` | comma-separated set of flag names |
| `TRC_ID`, `monomer_size`, `HOR_status`, `HOR_confidence` | KITE annotation joined back |

## Threshold parameters

All thresholds are CLI parameters, exposed by `consensus_ensemble.R`.

| parameter | default | purpose |
|---|---:|---|
| `--cov-strict` | 0.90 | grade A floor |
| `--cov-qualified` | 0.80 | grade B floor |
| `--cov-low` | 0.50 | grade C floor |
| `--core-cutoff` | 0.95 | `boundary_overext` trigger (core_fraction must drop below this) |
| `--gap-cutoff-bp` | 1000 | `internal_gap` trigger |
| `--pident-sd-cutoff` | 1.5 | `heterogeneous` trigger |
| `--pident-p05-cutoff` | 92.0 | `low_pident` trigger |

Defaults are conservative starting points; re-tune per genome if
needed (lower `cov_strict` for permissive grading, raise
`pident_sd_cutoff` to suppress the heterogeneous flag in
HOR-rich genomes, etc.).

## Usage

Inputs: an existing `TideCluster.py run_all` output directory. The
following files are read:

- `prefix_tarean/fasta/<TRC>.fasta`
- `prefix_kite/monomer_size_top3_estimats.csv`
- `prefix_tidehunter.gff3`
- `prefix_clustering.gff3`

Three steps (each can be re-run independently because outputs are
cached):

```bash
# 1. Array-MSA consensus per TRA + self-BLAST validation.
Rscript tarean/consensus_prototype/consensus_msa.R --mode batch \
  --kite-tsv   <prefix>_kite/monomer_size_top3_estimats.csv \
  --tarean-dir <prefix>_tarean/fasta \
  --out-dir    out/msa \
  --cpu 4
Rscript tarean/consensus_prototype/tests/validate_drapa_blast.R \
  --consensus-fasta out/msa/per_tra_consensus.fasta \
  --diagnostics     out/msa/per_tra_diagnostics.tsv \
  --tarean-dir      <prefix>_tarean/fasta \
  --kite-tsv        <prefix>_kite/monomer_size_top3_estimats.csv \
  --out-dir         out/msa/blast_validation \
  --cpu 4

# 2. TideHunter consensus per TRA + self-BLAST validation.
Rscript tarean/consensus_prototype/tests/validate_th_single.R \
  --tidehunter-gff3 <prefix>_tidehunter.gff3 \
  --clustering-gff3 <prefix>_clustering.gff3 \
  --tarean-dir      <prefix>_tarean/fasta \
  --kite-tsv        <prefix>_kite/monomer_size_top3_estimats.csv \
  --out-dir         out/th \
  --cpu 4

# 3. Consensus selector + quality grades + flags.
Rscript tarean/consensus_prototype/consensus_ensemble.R \
  --mb-fasta     out/msa/per_tra_consensus.fasta \
  --mb-blast-tsv out/msa/blast_validation/blast.tsv \
  --th-fasta     out/th/per_tra_th_consensus.fasta \
  --th-blast-tsv out/th/blast.tsv \
  --kite-tsv     <prefix>_kite/monomer_size_top3_estimats.csv \
  --out-dir      out/selected
```

Wall time scales with array DNA volume. The two BLAST passes (steps
1 and 2) dominate; the selector itself is a few seconds because it
re-uses the cached BLAST outputs.

## Outputs

Under `out/selected/`:

- `per_tra_consensus.fasta` — the selected consensus per TRA, with the
  source method (`M_B` / `TH_single`) annotated in the header.
- `per_tra_metrics.tsv` — every per-TRA metric column listed above.
- `summary.log` — totals, threshold values, grade counts, flag counts.

## Limits

- 7 % of TRAs land at grade D, dominated by SSRs (m ≤ 50). MSA-based
  approaches do not fit those; they are flagged in the existing
  pipeline as `repeat_type=SSR` and should be routed to a motif-only
  consensus (the motif is already detected at the clustering step).
- ~1.5 % of TRAs are well-formed but resistant: heterogeneous HOR
  arrays where a single consensus cannot tile every variant, or arrays
  with kb-scale internal gaps that look like a different repeat
  embedded inside the dominant one. The flags surface these
  explicitly so downstream tools can decide how to handle them.
- The selector is read-only relative to its inputs. Iterating on
  thresholds, flag rules, or selection logic does not require
  re-running MSA or BLAST.

## See also

- `docs/hor_classification.md` — KITE's HOR classification, source of
  the `HOR_status` annotation joined into the per-TRA metrics.
- `docs/kite_consensus_review.md`,
  `docs/kite_consensus_implementation_plan.md` — earlier design notes
  on consensus approaches that motivated this work.
