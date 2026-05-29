# KITE → kitehor integration

TideCluster 1.10.0 replaces the in-tree R KITE analysis
(`tarean/kite.R`, removed) with [kavonrtep/kitehor](https://github.com/kavonrtep/kitehor),
a Rust reimplementation of the k-mer-interval principle that emits
substantially richer per-array structural information (HOR multiplicity,
subrepeat candidates, SSR scan, per-peak identity rescoring).

**Current kitehor pin:** `0.12.0` (from `-c petrnovak`, see
`conda-deps.txt`).

The integration evolved through three kitehor releases — see the
[Implementation history](#implementation-history) appendix at the
bottom for the transitions.

## Pipeline overview

`TideCluster.py` runs three kitehor subcommands per KITE-input fasta:

```
kitehor kite-periodicity --periodogram --out-peaks
    → kitehor.kite.tsv + kitehor.kite.peaks.tsv + kitehor.periodogram
# then in parallel:
kitehor rescore --max-period 10000 --top-n 20
    → kitehor.rescored.peaks.tsv
kitehor ssr-scan --kite-peaks
    → kitehor.ssr.tsv + kitehor.ssr.regions.tsv
```

No `analyze`, no `rule-classify`, no `tandem-validate`, no
`summary-merge` — those subcommands belong to the older
combined_class-cascade flow.

Two CLI flags on `run_all` / `tarean` expose the rescore caps:

| Flag                          | Default | Meaning                                         |
|-------------------------------|---------|-------------------------------------------------|
| `--kite-rescore-max-period`   | 10000   | `rescore --max-period` cap (bp)                 |
| `--kite-rescore-top-n`        | 20      | `rescore --top-n` cap (peaks per array)         |

Outputs land under `<prefix>_kite/`:

| File                              | Role                                                  |
|-----------------------------------|-------------------------------------------------------|
| `kitehor.kite.tsv`                | Per-array top-3 monomer-size + score                  |
| `kitehor.kite.peaks.tsv`          | Per-array ranked kite peaks                           |
| `kitehor.periodogram`             | FASTA-like bundle consumed by the heatmap renderer    |
| `kitehor.rescored.peaks.tsv`      | **Canonical per-peak diagnostic** (24 cols × every peak per array) |
| `kitehor.ssr.tsv`                 | SSR scan summary per array                            |
| `kitehor.ssr.regions.tsv`         | SSR per-region detail                                 |
| `monomer_size_top3_estimats.csv`  | TideCluster-derived per-array summary (see schema below). Filename kept for backward compat with `tc_per_tra_consensus.py` + the R consensus prototype. |

Heatmap PNGs are still rendered in R (`tarean/kite_heatmaps.R`, base
R only — no new conda deps) from the periodogram bundle.

## Founder reassignment

rescore's `founder_period` column is the peak with the highest
`identity_med` on a record. On HOR arrays that's typically the **HOR
tile** (e.g. 1512 bp on CEN6's classic 503/1512/×3 case), not the
base monomer. TideCluster reassigns:

- **Strongest** = rescore's `founder_period`.
- **Founder** is reassigned by a four-pass algorithm
  (`tc_utils.build_monomer_size_csv`):

  1. **Pass 1 — per-TRA strict.** Smallest divisor P with
     `id_med ≥ 0.7`, integer k = Strongest / P in `[2, 30]`, ratio
     tolerance ±0.05.
  2. **Pass 2 — TRC consensus rescue.** Per TRC, compute the
     consensus founder as the median of the largest ±10 % cluster of
     Pass-1 founders (gated by ≥ 3 arrays AND ≥ 25 % of the TRC's
     reassigned arrays). Arrays where Pass 1 left founder = Strongest
     look for a rescored peak within `max(±10 %, ±20 bp)` of the
     consensus with `id_med ≥ 0.5`. The rescue is only adopted when
     the rounded multiplicity ≥ 2. Flagged
     `irregular_multiplicity = true`; raw fractional k stored in
     `multiplicity_raw`.
  3. **Pass 3 — solo / no-consensus relaxed individual.** Per-array
     ratio tolerance relaxed to ±0.20, still `id_med ≥ 0.7`,
     k ∈ [2, 30]. Flagged irregular when it fires.
  4. **Pass 4 — SSR-founder override.** When
     `ssr_total_coverage_pct ≥ 95 %` and a dominant SSR motif is
     present, force founder = Strongest = kite top peak (rank 1),
     multiplicity = 1. On SSR-pure arrays rescore's identity rule
     picks a noise peak (every multiple-of-motif period aligns at
     id_med ≈ 1.00 on `(motif)ₙ`, but the kite-score ratio between
     the SSR period and the noise peaks is ~10³). Flagged
     `ssr_founder_override = true`.

- **Fallback** (separate concept): rescore returned NA
  `founder_period` (no peak passed 0.7, or every peak above
  `--max-period`) → founder = Strongest = kite rank-1 peak;
  `founder_fallback = true`. The SSR-override flag supersedes
  fallback on SSR-pure arrays.

**Validated on CEN6:** TRC_1 strict-only 7/33 arrays reassigned → 19
after the rescue passes; distribution `×2 (12), ×3 (2), ×4 (1),
×5 (3)`; prevalent founder 502 bp in 31/33 arrays.

**Validated on drapa** (33-scaffold, 392 Mb assembly): 127 SSR-pure
arrays correctly flagged across 9 TRCs; 15 founder periods flipped
by Pass 4 (the actually-broken ones); 3 spurious HOR calls demoted.

## Subrepeat tiering (per peak)

Per-peak classifier (`tc_utils._classify_subrepeat_tier`, mirrored
in `tc_rerender_report.classify_peak_tier`). Same logic, two copies
— kept synchronised by hand.

| Tier              | Rule                                                                                              |
|-------------------|---------------------------------------------------------------------------------------------------|
| `HIGH`            | period ≤ founder/4, `scan_occupancy_frac ≥ 0.15`, AND (`subrepeat=true` OR `phaseC ≥ 0.10` OR `autoF ≥ 0.4`) |
| `LIKELY`          | period ≤ founder/4, `scan_occupancy_frac ≥ 0.20`, `scan_n_intervals ≥ 10`                         |
| `KMER_SUPPORT`    | period ≤ founder/4, k-mer signal only (no per-base scan)                                          |
| `AMBIGUOUS`       | 0.25 < ratio ≤ 0.33 with mild support                                                             |
| `WEAK`            | period ≤ founder/4 but no signal                                                                  |
| `OBSERVATIONAL`   | founder NA but `scan_occupancy_frac ≥ 0.05`                                                       |
| `REJECT_*`        | phantom / `period > founder/3` / founder-itself / scan_occ=0 with founder known / weak ambig      |

Per-TRA report table surfaces **HIGH + LIKELY only**, top-2 by
`scan_occupancy_frac` desc. The Details child row shows every
non-rejected peak with its tier.

Subrepeats are intentionally rare: a real nested subrepeat is a
short motif tiling inside the founder monomer but occupying only
part of it. Full occupancy with integer multiplicity is HOR, not
subrepeat. Validated rare-positive profile on
`test_data/IPIP200579_2026-04-14`: 40 / 3024 records (~1.3 %) carry
HIGH/LIKELY calls; canonical TRC_104 (founder 180, sub 36,
occ ~0.5) and TRC_666 (founder 250, sub 36) all fire correctly.

## Per-array CSV schema (`monomer_size_top3_estimats.csv`)

```
TRC_ID seqid start end array_length
monomer_size score   monomer_size_2..5 score_2..5
hor_status hor_confidence              ← empty cells (R HOR_* alias shim)
founder_period strongest_period multiplicity
multiplicity_raw irregular_multiplicity
delta_id_pp founder_id_med strongest_id_med
founder_fallback ssr_founder_override
subrepeat_1_period subrepeat_1_occ subrepeat_1_tier
subrepeat_2_period subrepeat_2_occ subrepeat_2_tier
ssr_flag ssr_dominant_motif ssr_dominant_motif_coverage_pct
ssr_total_coverage_pct ssr_top_motifs consensus_period_bp
```

The full per-peak rescore diagnostic table lives alongside in
`kitehor.rescored.peaks.tsv` (24 cols × every peak per array). The
report reads it directly for the per-array Details child row.

`hor_status` / `hor_confidence` cells are intentionally emitted
empty so the per-TRA consensus R script's `HOR_*` alias shim doesn't
choke; the same shim covers `combined_class` and the 0.10.0 `tv_*`
columns on legacy CSVs.

## Report (v2) — per-TRA table

```
▸ # · seqid · start · end · length · Founder · Δid · Strongest · ×k
  · Other periods · Subrepeat · SSR
```

- **Founder** cell: red `*` for `founder_fallback`; lavender **SSR**
  pill for `ssr_founder_override` (suppresses the `*`).
- **×k** cell: amber `~` pill when `irregular_multiplicity = true`;
  hovering shows raw k (e.g. `k = 4.78`).
- **Subrepeat**: coloured tier pills (HIGH green, LIKELY teal).
- **Other periods**: suppressed on SSR-overridden rows (would be
  flooded with `(1.00)` noise peaks on `(motif)ₙ`).
- `▸` opens a DataTables child row with the founder/strongest
  summary (with `irreg (k = X.XX)` line when irregular), the SSR
  scan, and the full per-peak diagnostic table.

**Per-TRC card** (Classification section) surfaces:
- Prevalent founder (e.g. `502 bp (in 31/33 arrays)`).
- HOR multiplicities (`×2 (12), ×3 (2), ×4 (1), ×5 (3)`).
- Counts: arrays with HOR call (×k≥2) / Subrepeat / SSR / Fallback /
  Irregular.

**Per-TRC genome-distribution ideogram** colours each array by the
dominant structural signal: HOR (green), Subrepeat (teal), SSR
(purple), Fallback (red), Plain TR (grey). Legend matches the
`_class_fill()` v0.12 cascade.

## Code touched

| File                                       | Role                                                                                                                                                                  |
|--------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `TideCluster.py`                           | kitehor invocation block (kite-periodicity + rescore + ssr-scan), the two `--kite-rescore-*` CLI flags                                                                |
| `tc_utils.build_monomer_size_csv`          | Four-pass founder algorithm, subrepeat tier classifier, SSR aggregation, output schema                                                                                |
| `tc_rerender_report.py`                    | v0.12 per-TRA table, Details child row, TRC aggregates, cluster-overview bubble chart, TRC distribution plot on index, legacy load paths for 0.9.x / 0.10.0 / R-kite  |
| `tarean/assets/tidecluster.js`             | Tooltip generalisation, Details child row click handler, index distribution selector + search                                                                         |
| `tarean/assets/tidecluster.css`            | Tier pills, irregular pill, SSR pill, fallback marker, distribution-plot grid, sectioned summary card, legacy class-badge palette                                     |
| `tarean/kite_heatmaps.R`                   | Heatmap PNGs from the `.periodogram` bundle                                                                                                                           |
| `conda-deps.txt`                           | Pinned `kitehor=0.12.0`                                                                                                                                               |
| `changelog.md`                             | TideCluster 1.10.0 entry                                                                                                                                              |

## Conda + Singularity

`kitehor` is shipped via the `petrnovak` conda channel
(`-c conda-forge -c bioconda -c petrnovak`). Bioconda recipe is a
deferred opportunistic contribution.

`TideCluster.def` reads `conda-deps.txt` so the SIF picks up the
pin automatically on rebuild.

---

## Implementation history

The kitehor integration was built in three increments. Each release
introduced a breaking change to the per-array signal schema;
TideCluster keeps load paths for all three so older output dirs
still rerender via `tc_rerender_report.py`.

### kitehor 0.9.3 — initial cut (commit `cb3602e`)

Pipeline: `kitehor analyze --periodogram` (single command, parallel
internally). Canonical merged output: `.summary.tsv` (32 cols) with
a `combined_class` column over **8 categories**: `hor`,
`hor_with_ssr`, `tr`, `tr_with_ssr`, `tr_with_subrepeat`,
**`tr_with_nested_tr`**, `pure_ssr`, `unresolved`. Subrepeat scan
emitted `.subrepeat.tsv`; HOR-within-tile diagnostics emitted
`.hor_within_tile.tsv`.

### kitehor 0.10.0 — schema reshuffle (commit `14a0761`)

**Class set:** 8 → 7. `tr_with_nested_tr` removed (those arrays
fall to `tr`, `tr_with_subrepeat`, or `unresolved`); `hor_with_ssr`
retained.

**`summary.tsv`:** still 32 cols but reshuffled. Dropped:
`length_bp`, all `subrepeat_*`, `density_hint`, `founder_density`,
`phase_contrast`, `density_n_windows`. Added 10 `tv_*`
(tandem_validate) columns. Subrepeat signal moved to:
host monomer = `tv_host_period`, subrepeat period =
`tv_best_candidate_period`.

**Subcommand outputs:** 9 TSVs → 6. `.tandem_validate.tsv` replaces
`.subrepeat.tsv` + `.hor_within_tile.tsv`; `.windows.tsv` removed.
Removed subcommands: `subrepeat-scan`, `hor-validate`.

### kitehor 0.12.0 — combined_class retired, current state (commit `d358bda` + follow-ups)

The `combined_class` cascade was retired in TideCluster (v0.11.0
added `hor_with_ssr` / `unresolved_with_ssr` classes and switched
`pure_ssr` to read `ssr_raw_total_coverage_pct`, but TideCluster
never shipped with it). The current flow uses the three kitehor
subcommands documented in the [Pipeline overview](#pipeline-overview)
above; the four-pass founder algorithm
([Founder reassignment](#founder-reassignment)) replaces the cascade
entirely.

The v2 report (`tc_rerender_report.py`) was overhauled to surface
per-array structural signals (Founder / Strongest / ×k / Subrepeat /
SSR / Fallback / SSR-override) instead of the now-defunct
`combined_class` roll-ups. The legacy 0.10.x (`tv_*`,
`combined_class`) and 0.9.x (`subrepeat_*`, `tr_with_nested_tr`)
columns are still consumed by the loader when present so old output
dirs rerender unchanged.
