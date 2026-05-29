## 1.10.1 (2026-05-29)
- Fix Python 3.11 incompatibility in the report v2 builder: a
  backslash inside an f-string expression part (per-array Details
  child row Founder line) raised a `SyntaxError` on the CI runner
  (PEP 701 relaxed this in 3.12, which is why local tests passed).
  The exception was swallowed by `_build_report_v2`'s broad
  try/except, surfacing as a stderr WARNING and an empty
  `<prefix>_index.html` — the `tests/short.sh` `[ -s short_index.html ]`
  guard then failed with `FAIL: no index.html`. Targeted fix: hoist
  the conditional out of the f-string. Verified by an AST scanner
  across `tc_rerender_report.py`, `tc_utils.py`, `TideCluster.py` —
  this was the only instance.
- CI: bump `python-version` 3.11 → 3.12 across all three workflows
  (`tests.yml`, `release.yml`, `conda-release.yml`) so 3.12-only
  syntax can't slip past local-dev testing again. `conda-deps.txt`
  still has `python>=3.6` — the bump is CI-side only.

## 1.10.0 (2026-05-29)
- KITE step now runs [kitehor](https://github.com/kavonrtep/kitehor)
  `0.12.0` (Rust) instead of the R `tarean/kite.R` script. The slow
  cascade (`rule-classify` → `tandem-validate` → `summary-merge`)
  is dropped in favour of a leaner two-stage flow:
    1. `kitehor kite-periodicity --periodogram --out-peaks` — k-mer
       interval scan + periodogram bundle + long-format peaks TSV.
    2. In parallel: `kitehor rescore` (banded semi-global alignment
       per kite peak; appends `identity_med`, `identity_iqr`,
       `phantom`, `subrepeat`, `coverage_frac`, `scan_*`,
       `kmer_autocorr_founder`, `kmer_phase_contrast`, and a per-
       record `founder_period`) and `kitehor ssr-scan`.
  Two new CLI flags on `TideCluster.py run_all`/`tarean`:
  `--kite-rescore-max-period` (default 10000) and `--kite-rescore-top-n`
  (default 20).
- Founder / Strongest / Multiplicity are now derived in TideCluster in
  **three passes** that combine per-array strict reassignment with
  TRC-level cluster context:
    - **Pass 1 (per-TRA strict)**: **Strongest** = rescore's
      `founder_period` (highest-identity peak). **Founder** = the
      smallest peak P with `identity_med(P) ≥ 0.7` whose period is an
      integer divisor of Strongest (k = Strongest/P, 2 ≤ k ≤ 30,
      ±0.05 tolerance) — recovers the v0.10-style HOR base / tile
      decomposition.
    - **Pass 2 (TRC consensus rescue)**: per TRC, computes the
      consensus founder (median of the largest ±10 % cluster across
      Pass-1 successes, gated by ≥ 3 arrays AND ≥ 25 % of the TRC's
      reassigned arrays). Arrays where Pass 1 left founder = strongest
      look for a rescored peak within max(±10 %, ±20 bp) of the
      consensus with `id_med ≥ 0.5`; rescues fire only when the
      resulting multiplicity rounds to ≥ 2. Marked
      `irregular_multiplicity = true`. On CEN6 this nearly triples the
      HOR call count (7 → 19) by recovering arrays with irregular k.
    - **Pass 3 (solo / no-consensus fallback)**: when Pass 2 doesn't
      help, a relaxed per-array tolerance (±0.20 of integer, still
      `id_med ≥ 0.7`) fires and is also flagged irregular.
  When rescore returns NA founder_period (no peak passes 0.7 or every
  peak above `--kite-rescore-max-period`), TideCluster falls back to
  the top-scored kite peak (rank 1) and marks the array as a
  `founder_fallback` in the report (separate from `irregular_multiplicity`).
- Subrepeat tiering is ported from kitehor's `docs/rule_proto.md`
  cheat-sheet: per peak, tier ∈ `{HIGH, LIKELY, KMER_SUPPORT,
  AMBIGUOUS, WEAK, OBSERVATIONAL, REJECT_*}`. The per-array
  "Subrepeat" column surfaces the top 2 HIGH/LIKELY candidates by
  `scan_occupancy_frac` desc; the Details child row labels every
  non-rejected peak with its tier. Strict by design — subrepeats
  are intentionally rare.
- Report (v2) — per-TRA table redesigned:
    ▸ # · seqid · start · end · length · Founder · Δid · Strongest
    · ×k · Other periods · Subrepeat · SSR
  Coloured tier pills mark subrepeat strength; a click on the ▸
  control opens a Details child row with the full rescore diagnostic
  table (every peak: rank, period, score, id_med, id_iqr, occ,
  scan_n, cov_frac, spatial, autoF, phaseC, subrepeat, phantom,
  tier). Founder cells get a red `*` when the rescore fallback fired.
- Report (v2) — roll-ups updated: the per-TRC Classification card,
  the index summary, and the KITE overview swap the old
  `combined_class` distribution / Class mix columns for per-array
  counters (HOR ×k≥2, Subrepeat, SSR, Fallback). The per-TRC card
  additionally surfaces the **Prevalent founder** (median of the
  largest ±10 % founder cluster within the TRC, with the supporting
  array count) and the **HOR multiplicities** distribution
  (e.g. `×2 (12), ×3 (2), ×4 (1), ×5 (3)`) when v0.12 data is
  present. Genome-distribution ideograms colour by the dominant
  signal per array (HOR green, subrepeat teal, SSR purple,
  fallback red, simple TR neutral grey).
- `monomer_size_top3_estimats.csv` schema (still consumed by
  `tc_per_tra_consensus.py` and the R consensus prototype):
  - HOR columns (`hor_status`, `hor_confidence`) kept as empty cells
    for the per-TRA consensus R script's `HOR_*` alias shim.
  - Top-5 monomer-size peaks by score (`monomer_size` ..
    `monomer_size_5` + `score*`), from the rescored peaks file.
  - New columns: `founder_period`, `strongest_period`, `multiplicity`,
    `multiplicity_raw` (fractional k for irregular cases),
    `irregular_multiplicity` (true when Pass 2 / Pass 3 fired),
    `delta_id_pp`, `founder_id_med`, `strongest_id_med`,
    `founder_fallback`, `subrepeat_{1,2}_period`,
    `subrepeat_{1,2}_occ`, `subrepeat_{1,2}_tier`, and the SSR
    columns (`ssr_flag`, `ssr_dominant_motif`,
    `ssr_dominant_motif_coverage_pct`, `ssr_total_coverage_pct`,
    `ssr_top_motifs`, `consensus_period_bp`).
  - The full rescored peaks TSV (`kitehor.rescored.peaks.tsv`, 24
    cols × every peak per array) is the source of truth for the
    Details child row and stays alongside.
  - Rerender still reads the kitehor 0.10.x (`combined_class` +
    `tv_*`) and 0.9.x (`subrepeat_*` + `tr_with_nested_tr`) schemas
    plus the legacy capitalised `HOR_*` schema, so archival output
    directories from earlier versions still render with the
    appropriate roll-ups.
- Per-TRC heatmap PNGs (`profile_plots/profile_*.png`,
  `profile_plots/profile_top3_*.png`) still rendered by
  `tarean/kite_heatmaps.R` from kitehor's `--periodogram` bundle.
- **Pass 4 SSR-founder override.** On SSR-pure arrays
  (`ssr_total_coverage_pct ≥ 95 %` with a dominant motif), rescore's
  identity-driven founder pick is unreliable: the SSR motif period
  sits below rescore's min-period (so id_med is NA) while every long
  multiple-of-motif period rescore *can* evaluate has near-perfect
  identity on the underlying `(motif)ₙ` string. Pass 4 forces
  founder = strongest = kite top peak (rank 1), multiplicity = 1, and
  stamps a new `ssr_founder_override = true` column. A lavender
  **SSR** badge in the Founder cell marks the provenance and
  supersedes the red `*` fallback marker (rescore-NA on a
  sub-min-period SSR motif is the expected case, not an anomaly).
  Validated on drapa: 127 SSR-pure arrays flagged (matching the
  source data), 15 spurious founder periods flipped, 3 spurious HOR
  calls demoted, 112 redundant fallback flags cleared.
- **Index page redesign.** REPORTS tile removed; RUN SUMMARY card
  reorganised into three labelled subsections — Input / Clusters
  (TRC) / Arrays (TRA) — with SSR enrichment (# arrays with SSR
  signal, # SSR-dominant ≥ 50 %, total SSR bp, top-3 SSR motifs by
  bp covered). "TAREAN-analysed" relabelled to "with TAREAN
  consensus" since every array goes through KITE.
- **TRC distribution plot** (new section on the index page). Hand-
  rolled SVG ideogram showing every TRA across the assembly, with a
  left-side TRC selector + live search + "Show all" reset. Default
  view is muted grey; clicking a TRC lifts only its arrays to red
  with a thin outline. Contigs < 1 Mb hidden, capped at 80 longest.
- **TAREAN tab cleanup.** Dropped the H/S/X (HOR/subrep/SSR tally)
  column on v0.12 reports — the same signals live on the per-TRC
  dashboard. Type column simplified to TR or SSR. Page widened to
  1560 px so Graph / Consensus columns no longer scroll on standard
  screens.
- **Legends rewritten for v0.12.** KITE tab legend (`kite.md`) and
  per-TRC card column legend (`arrays_legend`, now a collapsible
  `<details>` block) drop combined_class references and describe the
  v0.12 founder / strongest / ×k / subrepeat / SSR signals. Per-TRC
  ideogram legend lists the five colour swatches `_class_fill`
  actually emits in v0.12 mode (HOR / Subrepeat / SSR / Fallback /
  Plain TR).
- **Bug fixes.** Other-periods identity formatter capped at `1.00`
  (was emitting `.100` for any identity in [.995, 1.00] due to the
  `02d` minimum-width format). Multi-line markdown bullets and
  blockquotes now render correctly in the About / Credits sections.
- New conda runtime dependency: `kitehor=0.12.0` (requires
  `-c petrnovak`). No new Python plotting deps.
- Removed: `tarean/kite.R` (736 lines).

## 1.9.0 (2026-04-21)
- KITE now classifies each tandem repeat array into one of three HOR
  (higher-order repeat) categories — `No HOR detected`, `HOR-visible`,
  `HOR-dominant` — using integer-multiple relationships between the
  top three monomer-size estimates. Per-TRC counts in each category
  are reported in the KITE HTML summary and TSV. Thresholds exposed
  as constants at the top of `tarean/kite.R`.
- New `tc_rerender_report.py` CLI builds an opt-in report v2 from an
  existing TideCluster output directory. Emits `<prefix>_report_v2/`
  containing modern DataTables-powered HTML pages (summary, all TRCs,
  TAREAN, KITE, superfamilies) plus one per-TRC dashboard each,
  offline-first, no new conda dependencies. The original reports are
  left untouched.
- Pipeline now writes `<prefix>_pipeline_stats.json` alongside the
  index HTML so downstream tools have the summary numbers without
  having to scrape HTML.
- Test suite gains a `rerender` tier.

## 1.8.3 (2026-04-18)
- Conda releases are now built from this repository via GitHub Actions on
  tag push; the external `kavonrtep/recipes/tidecluster` recipe is retired.
- Added `tests.sh` dispatcher with smoke/short/long tiers gating CI and
  release builds.

## 1.7.1 (2025-10-09)
- Gzip FASTA input support added
- Interactive visualization and comparative analysis HTML reports implemented
- Comparative analysis of multiple assemblies was added

## 1.6.2 (2025-07-10)
- RepeatMasker wrapper function added to handle long sequence names
- Documentation updated with additional output files
## 1.6.1 (2025-05-03)
- Bugfix in tarean run
- Extract consensus functionality added
## 1.6 (2024-11-05)
- TRC superfamily clustering added (consensus-vs-consensus BLAST analysis)
- HTML report improvements and bugfixes
- Documentation updated
## 1.5 (2024-03-07)
- K-mer interval tandem repeat estimation (KITE) added, this enables identification of higher order tandem repeats
- Analysis include estimation of monomer size for individual tandem repeat arrays
- Output of KITE analysis added to the HTML report
- Bug fix in sorting of output CSV tables

## 1.4 (2023-11-23)
- Better memory efficiency
- Limiting MMseq2 memory usage
## 1.3 (2023-10-23)
- Improved input FASTA parsing to decrease memory usage
## 1.2 (2023-10-06)
- Bug fix in SSRs detection
- Better SSRs reporting TAREAN output
## 1.1 (2023-09-07)
- Bugfix - better handling of SSR
## 1.0 (2023-08-17)
- Bugfix in HTML report
- Tidehunter run in batches to decrease memory usage
- Documentation updated
## 0.9 (2023-08-11)
- Documentation updated
- Info related to bug #1 added to documentation
- Exporting TR library added
## 0.0.8 (2023-05-02)
- Documentation update, HTML bugfix
- Workflow scheme added
- Workflow scheme formatted
- Update_gff script added
## 0.0.7 (2023-04-26)
- Output files added
## 0.0.6 (2023-04-25)
- Bug-fix - handle lower case character in DNA sequences
- --min_total_size option added, some bugfixes
## 0.0.5 (2023-04-17)
- Bug fix in logo plotting
- Clustering of SSRs added
- Format of output GFF corrected
## 0.0.4 (2023-04-05)
- Bugfixes, better defaults, documentation updated
## 0.0.3 (2023-03-23)
- Run all option added, bugfixes
## 0.0.2 (2023-03-23)
- Bug fixes
