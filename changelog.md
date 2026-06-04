## Unreleased
- Fix per-array `strongest_period` semantics in `tc_utils.build_
  monomer_size_csv()` Pass 1. The previous code trusted kitehor's
  `founder_period` column as the strongest periodicity, but kitehor's
  rescore picks by single-peak kite score on some arrays rather than
  by identity, occasionally returning the rank-1-by-score peak (low
  id_med) even when a clean integer multiple of it has higher id_med
  (the real HOR period). On *S. pimpinellifolium* TRC_8 chr12:1314680
  and chr12:1327272 this collapsed `multiplicity` from the true 10
  down to 1 — no HOR call. Strongest is now defined per the project's
  documented semantics (see
  `hermit/.../memory/reference_founder_strongest_hor.md`): the peak
  with the highest `identity_med` among peaks that pass the founder
  identity gate (`_FOUNDER_ID_MIN = 0.7`). Pass 1's divisor search
  then naturally finds the basic monomer as the smallest valid
  integer divisor of that higher strongest, recovering the HOR call.
- New `kitehor_founder_period` column on
  `<prefix>_kite/monomer_size_top3_estimats.csv` (additive; older
  consumers ignore it): captures kitehor's original rescore pick so
  that disagreement with TideCluster's recomputed `strongest_period`
  is auditable per row. When the two agree, the column simply
  duplicates `strongest_period`; when they differ, the row is one of
  the cases where the new argmax(id_med) override fired.
- Defensive guard: if no peak passes the id gate even though kitehor
  reported a founder, TideCluster falls back to kitehor's pick rather
  than failing — preserves 1.12.0 behaviour on those edge rows.
- Pass 2/3/4/5 logic is unchanged. The TRC-consensus rescue in Pass 2
  (ce233d1, kitehor 0.12 era) still handles arrays where Pass 1 left
  `multiplicity = 1` due to an irregular (fractional `k`) HOR.
- Validation on `test_data/Solanum_pimpinellifolium` (433 arrays):
    - the two cases Petr cited (TRC_8 chr12:1314680 and chr12:1327272)
      now report `strongest ≈ 1786`, `founder = 178`, `multiplicity = 10`
      with `kitehor_founder_period = 178` documenting the override
    - 71 arrays gained an HOR call (multiplicity 1 → ≥2), 1 lost
    - 109 founder shifts are within ±10 % of the prior pick (same
      monomer family, slightly different specific bp chosen by
      argmax(id_med) — cosmetic)
    - 30 structural founder shifts are downstream consequences of the
      strongest change via the strict ±0.05 divisor tolerance —
      semantically correct per the founder/strongest/HOR definitions
- New `tests/test_strongest_by_identity.py` with 12 synthetic cases
  covering the override (fires when kitehor picked low-id_med),
  the no-op (kitehor was already argmax), the kitehor-NA fallback,
  and the all-below-gate defensive fallback. Wired into
  `tests/unit.sh`.

## 1.12.0 (2026-06-04)
- Expose the per-genome clustering and superfamily similarity
  thresholds on the CLI; previously hard-coded:
    - `--cluster_identity` (default 75) and `--cluster_coverage`
      (default 0.8) on `clustering` and `run_all` — the BLASTN
      percent-identity and shorter-sequence coverage gates that
      define array-level families (TRCs).
    - `--superfamily_score` (default 20) on `tarean` and `run_all`,
      threaded into `compare_trc_by_blast.R --score_threshold` — the
      consensus-vs-consensus BLAST score gate
      `(length · pident − gap_openings) / longer_consensus_length`
      that defines superfamilies (also governs the below-threshold
      fallback rescue).
  Defaults equal the former hard-coded constants, so default runs are
  unchanged.
- Validate the new thresholds before a run starts: out-of-range values
  abort with a clear error (a coverage fraction entered as a percent,
  identity outside (0, 100], a negative score); biologically unusual
  but valid values emit a warning and proceed.
- README: new *Advanced: tuning the clustering thresholds (expert use)*
  subsection documents the options with defaults and an effect table,
  noting the defaults are the tuned operating point and should be left
  unchanged for routine use.
- Fix: the `tarean` step printed its settings block by reading
  run_all-only args (`tidehunter_arguments`, `min_length`, `no_dust`,
  `library`) directly, so invoking `tarean` standalone without an
  accumulated `<prefix>_cmd_args.json` raised `AttributeError`. These
  reads now fall back to `n/a` via `getattr`, matching the existing
  handling of `long` / `original_fasta`; behaviour is unchanged when the
  side-car file is present.

## 1.11.0 (2026-06-03)
- Comparative analysis: expose the per-edge coverage convention via
  a new `--coverage_mode` flag on `tc_comparative_analysis.R`,
  controlling how the MMseqs2 alignment coverage threshold is
  applied when building the TRC similarity graph:
    - `max` (default; current behaviour) —
      `max(qcov, tcov) ≥ --min_coverage`. A hit covering most of
      *either* sequence survives, so a short basic monomer can
      attach to a longer family as a subset.
    - `min` (new, bidirectional) —
      `min(qcov, tcov) ≥ --min_coverage`. Both sequences must be
      well covered. Rejects short-vs-long subset hits — e.g. a
      53 bp satellite folding into the ~8 kb 45S rDNA family.
  The selected metric drives both the edge filter and the edge
  weight (`selected_coverage × pident`); both `max_cov` and
  `min_cov` are retained in `mmseqs2_results.tsv` so the threshold
  can be re-tuned post-hoc without re-running the search.
  Defaults preserve 1.10.x behaviour.
- Also exposes `--min_coverage` (default 0.8) and `--min_identity`
  (default 80) on the CLI; previously these were hard-coded inside
  `process_trc_analysis()`. The MMseqs2 search-time pre-filter
  stays an OR (a valid superset for either mode); exact filtering
  happens in R.
- README: new *How similarity between TRCs is evaluated* subsection
  walks through the MMseqs2 comparison in plain language (qcov +
  tcov + percent identity, the two thresholds, the max-vs-min
  combine choice with a worked example) and the three new flags
  are added to the command-line option list. The *Purpose* bullet
  that previously hard-coded "fast-greedy" now reflects all three
  supported clustering algorithms.

## 1.10.6 (2026-06-03)
- Fix 1.10.5 regression: when the new below-TAREAN-threshold
  superfamily fallback BLASTs the small-TRC dimers against the
  big-TRC dimer DB and finds zero above-e-value hits, blastn writes
  an empty output file and `read.table(blast_out, …)` in
  `run_blast()` (compare_trc_by_blast.R) raises
  `Error in read.table(...) : no lines available in input`. The
  script aborts before writing `<prefix>_trc_superfamilies.csv`, so
  the *native* SFs (from the main big-vs-big BLAST, which did find
  hits) are silently lost too — the surrounding TideCluster.py
  pipeline still exits 0 because the wrapper sees a missing CSV,
  not a failed step. Observed during the manuscript's set1/set2
  1.10.5 rerun on *S. pimpinellifolium* (55 below-threshold TRCs,
  11 big-TRC dimer records, zero qualifying matches; diagnosed in
  `docs/tidecluster_bug_superfamily_fallback_empty_blast.md`).
- Guard `run_blast()` against `file.size(blast_out) == 0` and
  return a clean 0-row data frame with the correct 14-column
  schema (character qseqid/sseqid, numeric pident/evalue/bitscore,
  integer everything else). The fallback section's existing
  `nrow(bl_fb2) > 0` guards then short-circuit to "no
  attachments" and the CSV is written with the native SFs intact.
  Also future-proofs the main self-vs-self BLAST against the
  (extremely rare) zero-hit case on tiny datasets — the existing
  `if (nrow(bl3) == 0)` early exit fires cleanly instead of
  crashing inside `read.table()`.
- Validation:
    - R-level: mocked `system()` to write an empty file, run_blast()
      returns a 0-row 14-col df with types matching `read.table()`;
      gsub/order/duplicated/nrow all operate cleanly on it.
    - Happy path on `analysis_kite15k`: byte-identical to 1.10.5
      (4 SFs, 3 fallback rescues, exit 0, 16-row CSV).
    - Empty-fallback end-to-end on a synthetic small-TRC with random
      sequence that produces zero BLAST hits: exit 0, "Fallback:
      no qualifying matches", CSV written with the native 2 SFs
      (11 rows), HTML produced.

## 1.10.5 (2026-06-02)
- New **below-TAREAN-threshold superfamily fallback rescue** in
  `tarean/compare_trc_by_blast.R`. TRCs whose array total length
  falls below TAREAN's `min_total_length` gate (default 50 kb) are
  skipped by `tarean.R` and so produce no consensus dimer — they are
  absent from `<prefix>_consensus_dimer_library.fasta` and never
  enter the main BLAST clustering step. On real datasets this leaves
  small TRCs unlinked from their biological superfamily even when
  RepeatMasker annotates them identically to a larger member (on
  *A. thaliana* `analysis_kite15k`: the 14 kb 5S_rDNA TRC_27 was
  stranded as a singleton despite TRC_3 also being 5S_rDNA at ≥ 98 %
  coverage).
- The fallback runs *after* the main BLAST clustering and uses the
  **raw TideHunter per-array consensus dimers** (preserved by the
  `clustering` step under `<prefix>_consensus/<TRC>_dimers.fasta`)
  as the query against the existing dimer-library DB. Same BLAST
  parameters, same score formula
  `(length * pident - gapopen) / max(qlen, slen)`, same threshold
  (20) as the main clustering — no new knobs. For each
  below-threshold TRC, the single best-scoring big-TRC hit becomes
  the attachment target:
    - if the target big-TRC already belongs to a multi-TRC SF, the
      small TRC joins that SF;
    - if the target is a singleton (had no above-threshold big-big
      match), the (big, small) pair is promoted to a brand-new
      size-2 SF.
- The specific small-TRC raw dimer that produced the best BLAST hit
  is added to the SF's dotplot input set, so the per-SF pairwise
  dotplot naturally renders the similarity that triggered the
  attachment.
- Optional belt-and-braces safety: `--annotation_tsv
  <prefix>_annotation.tsv` (passed by default from `TideCluster.py`)
  rejects (small, big) merges where both TRCs have a non-NA
  annotation and the annotations disagree. Zero rejections on the
  validation dataset.
- New CSV column `fallback` (boolean) on
  `<prefix>_trc_superfamilies.csv`. TRUE for rescue-attached TRCs;
  FALSE for native SF members. Older CSVs without the column are
  still read correctly by the v2 report.
- `tc_rerender_report.py` surfaces the new flag: per-TRC dashboard
  "Superfamily" stats row carries a `‡ via fallback` annotation on
  rescue-attached TRCs; the superfamilies page adds a top-level
  callout count plus per-SF `‡` daggers with provenance tooltip
  next to each fallback member.
- Validation on `tmp/analysis_kite15k/Arabidopsis_thaliana`:
    - SF1 unchanged (8 TRCs, satellite group)
    - SF2 augmented with TRC_26 (AthSat500, attach to existing) score=22.08
    - SF3 NEW {TRC_12, TRC_21} (uncharacterized, singleton-promotion) score=34.10
    - SF4 NEW {TRC_3, TRC_27} (5S_rDNA, singleton-promotion) score=47.86
- Runtime cost is negligible: the fallback BLAST is **0.15 s** in
  isolation on the validation set (21 small-TRC dimer records vs
  389 big-TRC dimer records, 4 threads). End-to-end
  `compare_trc_by_blast.R` is within filesystem-cache noise of the
  1.10.4 baseline (~4 min on this dataset; fallback overhead
  < 0.07 % of the clustering step, < 0.01 % of the full pipeline).

## 1.10.4 (2026-06-02)
- Fix per-TRC annotation coverage in
  `tc_utils.get_repeatmasker_annotation()`. An over-broad guard
  (`if seqid not in seq_rm_info:`) collapsed multi-hit RepeatMasker
  output to only the first hit per query, so multi-monomer consensus
  sequences reported ~one-monomer's worth of reference coverage
  (e.g. *A. thaliana* CEN180 15 %, 45S rDNA 0.2–42 %) instead of the
  true ~95–99 %. The position-sorted nature of `.out` always picked
  the left-most hit. Replace with per-(seqid, annotation) interval
  collection and report the merged-interval (union) length / sequence
  length, bounded to [0, 1], so tandem hits accumulate while
  overlapping hits are not double-counted. Family assignment is
  unchanged; only the reported proportion is corrected. Also fixes
  `RepeatMaskerFeature.length` off-by-one (end − start + 1). Covered
  by new synthetic tests `tests/test_annotation_coverage.py` wired
  into `tests/unit.sh`.
- Fix `%25` leaking through to HTML in the report-v2 renderer.
  `tc_utils.get_repeatmasker_annotation()` and
  `get_ssrs_description*()` emit percentages as `(N.N%25)` because
  the strings are written into GFF3 column-9 attribute values (where
  `%` must be URL-encoded per the GFF3 spec; `tarean_report.R`
  already decodes `%25` → `%` for both Annotation and SSRs columns).
  `tc_rerender_report.py` decoded `ssr_motif` but missed four other
  render sites: the scatter-plot tooltip, the TRC table Annotation
  column, the per-TRA arrays-table SSR column, and the TRC detail
  card Annotation row. Add the same `.replace("%25","%")` at all
  four. No change to GFF3 output — encoding stays spec-correct on
  disk; only the HTML render decodes for display.
- New **Pass 5 (peak-cluster fallback rescue)** in
  `tc_utils.build_monomer_size_csv()`. When rescore returns
  `founder_period=NA` the legacy Pass-1 fallback set founder =
  rank-1 peak by single-peak kite score. On noisy centromeric
  arrays this picks a spurious short-period peak (e.g. 310 bp) over
  the real ~10 kb HOR signal, which gets fragmented across many
  adjacent rescored periods (9383, 9386, 10461, 10466, 10468,
  10689, ...): each individual peak's score is small and its
  `identity_med` fails rescore's 0.7 threshold, but collectively
  they dominate the array with high `coverage_frac` (~0.4) and
  `spatial_contrast` (~1.0). Pass 5 clusters the rescored peaks by
  period (single-link, mixed ±5 % / ±100 bp window), sums their
  kite scores, and accepts the dominant cluster as founder iff
  `n_peaks ≥ 2 AND max(coverage_frac) ≥ 0.20 AND
  max(spatial_contrast) ≥ 0.5 AND score_sum ≥ rank-1 peak score`.
  Fires only on `fallback=true && !ssr_override` rows — Pass-1
  successes, Pass-2/3 rescues, and SSR overrides are not touched.
  Covered by `tests/test_cluster_rescue.py` (9 synthetic cases).
- New CSV columns in `monomer_size_top3_estimats.csv` (all
  additive; every consumer — `tc_per_tra_consensus.py`, the R
  consensus prototype, `kite_heatmaps.R`, `extract_consensus.R`,
  `tc_rerender_report.py` — reads by column name):
    - `founder_method` enum {`strict` | `none` | `pass2` | `pass3`
      | `ssr` | `cluster` | `fallback`}, recording which pass
      produced the final founder.
    - `cluster_rescue` boolean counterpart of
      `founder_method = cluster`.
    - `alt_cluster_1..2_{period,score_sum,n_peaks,cov_frac_max}` —
      the top-2 period clusters that are NOT the founder's, so the
      report can surface secondary signals (e.g. a 310 bp peak
      that Pass 5 correctly demoted but remains real biology).
- HTML report surfaces Pass 5: a `‡` Dagger sup-tag on cluster-
  rescued founders in the per-array Details summary; an "Other
  periodicities: 310 bp (Σ 0.027 · n=1) · …" line below the founder
  cell listing the top-2 alt clusters; a new "Arrays rescued by
  cluster (Pass 5)" counter on the TRC dashboard Classification
  card. Cluster-rescue clears `founder_fallback`, so the red `*`
  fallback marker and the `‡` rescue badge are mutually exclusive.
- Cluster-overview bubble chart on `<prefix>_index.html` now uses
  the median per-array `founder_period` for the x-axis (was: median
  rank-1 kite peak `m1`), falling back to `m1` only on legacy CSVs
  that lack the column. On *A. thaliana* centromeres TRC_2 and
  TRC_4 — both rescued by Pass 5 — move from x ≈ 310 bp to
  x ≈ 10.4 kb, matching the rest of the report (which already
  surfaces `founder_period` as the canonical monomer size). TRCs
  whose founder already matched `m1` (the majority — Pass-1
  successes on clean tandem repeats) are unaffected.
- Validation on `tmp/analysis_kite15k/Arabidopsis_thaliana`
  (*A. thaliana* centromeres, kitehor v0.12): 3 of 4 NA-founder
  rows flipped to cluster rescue (the TRC_2 and TRC_4 cases —
  founders 310/312 → 9384 / 10346 / 10464 with
  `founder_method=cluster`); the remaining 1 stayed fallback
  because no cluster passed the gate (truly no signal). 0 Pass-1-
  success rows had their founder period change.

## 1.10.3 (2026-06-01)
- Fix silent data loss in the comparative analysis
  (`tc_comparative_analysis.R`). TRCs whose consensus sequences
  produced no above-threshold MMseqs2 pair were never added as
  vertices to the similarity graph and were therefore omitted from
  every comparative output — `trc_satellite_families.tsv`, per-genome
  TRC counts in the HTML report, and the `Satellite_family`
  re-annotation of per-sample GFF3. The bias was toward genome-unique,
  divergent, large-monomer or very-short-monomer satellites — the
  cases a comparative analysis most needs to retain.
- Validation on the manuscript's 11-genome *Solanum* set
  (`set_1_2`): the union of `trc_satellite_families.tsv` and
  `ssrs_groups.tsv` now matches `tc_clustering.gff3` per sample
  (e.g. *S. lycopersicum* 53 → 73 TRCs in the satellite family table;
  535 TRCs recovered across 11 samples in total).
- New behaviour for pure-SSR TRCs: regrouped across genomes by
  shared SSR pattern (via `cluster_ssrs_sequences` output) rather
  than as isolated singletons, so a single SSR family spans every
  sample that carries the motif. Implemented in a new helper
  `apply_ssr_grouping()` which also keeps `trc_graph.rds` /
  `trc_graph.graphml` in sync with the reassigned `group_id`s.
- New post-condition check `check_trc_coverage()` warns at the end
  of the run if any TRC from `tc_clustering.gff3` is missing from
  both `trc_satellite_families.tsv` and `ssrs_groups.tsv`. Future
  regressions of this silent-drop class will surface in the log.
- Known follow-ups captured in
  `docs/comparative_analysis_followups.md`: (a) threshold mismatch
  between code (80 % identity / 0.8 max-coverage) and manuscript
  (70 % / 0.2); (b) short-monomer (2–7 bp) satellites that fall
  between the SSR and satellite classifiers. Neither is changed in
  this release.

## 1.10.2 (2026-05-29)
- Fix conda packaging: `kitehor=0.12.0` was missing from
  `conda/tidecluster/meta.yaml`'s `requirements.run`, so
  `mamba install tidecluster` did not pull it in. Users had to
  install kitehor manually for the kite step to work. Root cause:
  `meta.yaml` is hand-maintained and drifted from `conda-deps.txt`
  when the kitehor integration landed in 1.10.0.
- Add a drift-check helper `tests/check_conda_recipe_deps.py` that
  asserts every package in `conda-deps.txt` is also present in the
  recipe's `requirements.run`. Wired into `tests/smoke.sh` so future
  drift is caught on every push/PR.
- `tests/smoke.sh` now also calls `kitehor --version` so a missing
  KITE backend on PATH fails the source-tree smoke gate.
- `conda/tidecluster/meta.yaml` `test.commands` adds
  `kitehor --version` so conda-build's post-install test phase also
  catches a missing kitehor at tag time.

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
