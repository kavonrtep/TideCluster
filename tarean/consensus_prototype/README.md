# KITE per-TRA consensus — prototype

Stand-alone R prototype for extracting a per-TRA consensus sequence
on top of KITE's monomer-size estimates. Implements Option 2
(phase-aligned position frequency matrix / Fourier folding) from
`docs/kite_consensus_review.md`.

**Status:** prototype. **Not** wired into `TideCluster.py` or
`tarean/kite.R`. Validation summary in
`tmp/consensus_prototype_drapa/REPORT.md` (expected to be regenerated
by the run below). The current results do not meet the promotion
criteria in `docs/kite_consensus_implementation_plan.md` §6; see the
report for the reasoning and alternative paths.

## Usage

### Single-array mode (spot-checks)

```
Rscript tarean/consensus_prototype/consensus_pfm.R \
  --mode single \
  --fasta <TRA.fasta>       \
  --m     <monomer_length>  \
  --out   <prefix>
```

Produces `<prefix>_consensus.fasta`, `<prefix>_pfm.tsv`,
`<prefix>_diagnostics.tsv`.

### Batch mode (full run)

```
Rscript tarean/consensus_prototype/consensus_pfm.R \
  --mode batch \
  --kite-tsv    <prefix>_kite/monomer_size_top3_estimats.csv \
  --tarean-dir  <prefix>_tarean/fasta/ \
  --out-dir     <out_dir>     \
  --cpu         4
```

Produces under `<out_dir>`:

- `per_tra_consensus.fasta` — one record per TRA, name
  `TRC_<k>__<seqid>__<start>-<end>`.
- `per_tra_diagnostics.tsv` — length, coverage, mean Shannon entropy
  of the PFM, harmonic k₂–k₅ magnitudes per TRA.
- `summary.log` — total / ok / failed / wall-time.
- `failed.tsv` (only if any TRA failed) — reason per TRA.
- `per_tra_pfm/TRC_*.tsv` (only with `--write-pfm`) — per-TRA PFM.

### Tunables

- `--window-size` (default 50 000 bp)
- `--n-windows` (default 16)
- `--pseudo` (default 0.5) — pseudocount for PFM smoothing
- `--cpu` (default 4) — parallel workers in batch mode

## Tests

### Synthetic short-smoke

```
Rscript tarean/consensus_prototype/tests/make_synthetic.R
```

Builds five synthetic fixtures (clean, noisy, indel-drifted,
large-clean, HOR-2×) plus one N-only. Runs the single-mode CLI on
each and asserts pass criteria described in
`docs/kite_consensus_implementation_plan.md` §5.1. Output under
`tmp/consensus_prototype_smoke/`.

### Drapa batch validation

```
Rscript tarean/consensus_prototype/tests/validate_drapa.R
```

Runs the four validation questions in plan §5.2 on the output of
an earlier Drapa `run_all` (path hard-coded to
`tmp/drapa_run2/`). Prints results to stdout; the canonical write-
up is `tmp/consensus_prototype_drapa/REPORT.md`.

## Files

- `consensus_pfm.R` — CLI (single + batch).
- `lib/consensus_core.R` — algorithm helpers (windowing, base
  tally, PFM + consensus, column entropy, harmonic diagnostic,
  FASTA record-name parser).
- `tests/make_synthetic.R` — short-smoke harness.
- `tests/validate_drapa.R` — real-data validation driver.

## Design references

- Review (options A–E head-to-head): `docs/kite_consensus_review.md`.
- Implementation plan (Option 2 details + promotion path):
  `docs/kite_consensus_implementation_plan.md`.
