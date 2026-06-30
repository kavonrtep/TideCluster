## 1.16.1 (2026-06-30)
- **Comparative analysis: each TRC now maps to a single satellite family.**
  Fixed a bug where SSR-grouped TRCs split across communities were assigned to
  two families in `trc_satellite_families.tsv` (inflating family counts and
  double-counting length/annotation). `apply_ssr_grouping()` now reassigns every
  matching row, and `cluster_trc_sequences()` collapses `trc_groups` to one row
  per TRC.
- **New `--lowcomplexity_mask {0,1}` option for `tc_comparative_analysis.R`**
  (default `1` = unchanged). Set `0` to disable MMseqs2 tantan low-complexity
  masking in the all-vs-all search, so AT-rich / simple-repeat satellites whose
  monomers are mostly low-complexity stop fragmenting across many families.
- **Index report: clearer TRC highlight on the "TRC distribution across the
  assembly" ideogram.** The selected TRC's arrays are now drawn as a top-layer
  red overlay with a white halo (no more brownish blend), widened to a minimum
  width with a slight vertical overhang, while non-selected arrays dim.

## 1.16.0 (2026-06-26)
- **rDNA identification.** TRCs are now flagged as ribosomal DNA, distinguishing
  45S (18S/5.8S/25S) from 5S, by blastn similarity to a bundled reference
  library (`data/rdna_library.fasta`). A TRC is labelled from the best
  single-subunit reference coverage (robust to ITS/IGS dilution and ~90 %
  cross-species identity), consensus-first with a genomic-array fallback.
  Writes a clean `rDNA_type=45S|5S;rDNA_coverage=<frac>` into the clustering
  GFF3 plus a `<prefix>_rdna.tsv`, and surfaces it in the report (index list +
  per-TRC page). Default-on in `run_all`/`tarean` (`--no_rdna` to disable;
  `--rdna_library`/`--rdna_min_coverage`/`--rdna_min_identity` to tune), plus a
  standalone `rdna` subcommand to re-run on an existing run.
- **Cross-TRC overlap resolution.** The clustering GFF3 is now made
  non-overlapping across TRCs by default ("annotate each region once"): where
  variant arrays of a satellite cluster into separate TRCs that overlap at their
  boundaries, each contested span is assigned to the dominant TRC (largest total
  array length). No base of the union is lost. Pass `--keep_overlaps` to retain
  the raw overlapping regions.
- **`tc_reannotate` parallelism.** Replaced the single whole-genome
  `RepeatMasker -pa` (ineffective with a custom `-lib` on RMBlast; Dfam #274)
  with a chunked, pooled run — one single-threaded RepeatMasker per genome chunk
  in a process pool — so both the search and the single-threaded ProcessRepeats
  tail parallelise (`--chunk_size`/`--overlap`, default 50 Mb/100 kb). Sequences
  below `2*chunk_size` are only packed (byte-identical); larger ones split, with
  <0.15 % masked-bp drift at the cuts.
- Made the post-RepeatMasker `filter_intervals` sub-quadratic (per-(seqid,name)
  sweep with prefix-max + binary search); output byte-identical.

## 1.15.2 (2026-06-12)
- **Dominant-score ladder founder** (founder Pass 7b) for divergent satellites
  whose highest-identity peak is a long, *off-ladder* period. Pass 1 anchors the
  strongest peak on `argmax(identity_med)`; when that lands on a period that is
  not a clean integer multiple of the true monomer (e.g. strongest=3150,
  monomer=188, 3150/188=16.76), the divisor search either collapsed the founder
  to the long period (mult=1) or to a non-monomer `strongest/round(k)` value.
  Pass 7b instead anchors the ladder on the **rank-1-by-score** peak and adopts
  it as the founder when (a) it is incommensurate with the current founder,
  (b) it carries a genuine low harmonic ladder (m=1 and m=2 present, ≥3 distinct
  rungs, tight IQR/occupancy), and (c) kitehor's own `founder_period`
  corroborates it. `founder_method = "dominant_ladder"` (HOR-order tier `weak`).
  Recovers e.g. *Pisum* FabTR-53 arrays 3150/347/1808/168/165 → ~188, and aligns
  off-ladder outliers to their family monomer on *Solanum lycopersicum* (TRC_4
  48/50 → 53) and *Arabidopsis* (160 → 178). Validated genome-wide with
  `tools/founder_diff.py`; every flip is kitehor-corroborated with
  equal-or-higher identity, and the previously-validated 179×6 / 111×3 / 178×10
  calls are unchanged.
- **Annotation no longer crashes when RepeatMasker writes no `.out`** (empty
  input FASTA or no library hits, e.g. annotating omitted short regions on a
  contig with no hits). The missing output is now treated as zero annotations
  instead of aborting `run_all` before TAREAN.

## 1.15.1 (2026-06-11)
- **Harmonic-ladder founder** (founder Pass 7) for divergent-HOR satellites:
  when no single short peak clears the identity gate but the surviving peaks
  form a clean harmonic ladder (`p, 2p, 3p, …`), the founder is recovered as the
  ladder's fundamental `p` and tagged `founder_method="ladder"`. Recovers the
  basic monomer of high-divergence HOR arrays the strict/divisor passes miss
  (validated *A. thaliana* TRC_47 179×6 and *D. rapa* TRC_158 111×3). Purely
  additive — strict-path founder values are unchanged.
- **Report fix**: TRCs below the TAREAN size threshold (combined array length
  `< -M`, skipped by TAREAN) now still surface their KITE per-array founder
  results in the report instead of an empty section.
- **Report fix**: the per-array SSR section now trusts the pipeline's SSR call
  and shows **raw** per-array SSR coverage (`ssr_raw_*`) rather than the inflated
  kitehor consensus coverage (`ssr_consensus_*`).

## 1.15.0 (2026-06-10)
- **Problematic founder/SSR calls fixed** (root-cause analysis in
  `docs/problematic_calls_analysis.md`):
  - SSR founder now follows the **clustering** classification: a TRC tagged
    `repeat_type=SSR` (TRA consensus > 0.9 simple-sequence) gets its fundamental
    motif length as the founder for every array, regardless of kitehor coverage.
    The old kitehor per-array 95 % coverage override is removed (it mis-set a
    lone ATC-rich array in a satellite TRC). Fixes e.g. *S. lycopersicum*
    TRC_18/TRC_5 → 3 bp, TRC_2 → ~9491 bp.
  - **Harmonic-basis founder deepening**: when the rank-1-by-kite-score peak
    passes the identity gate, is shorter than the founder, and the founder is a
    clean integer multiple of it, the founder is deepened to that basic monomer
    (bypassing the divisor-search ceiling). Recovers e.g. TRC_4 → 53 bp.
  - **Prevalent-founder anchor** for long HORs: when an array's dominant peak is
    the TRC consensus founder, it is adopted under a **k-scaled** integer
    tolerance (`|k − round(k)| ≤ frac·k`); above k ≈ 50 the integer test is
    vacuous and the family consensus carries the call.
  - **Selective long-period rescore**: arrays whose dominant period exceeds the
    rescore cap are re-rescored at a higher cap (`--kite_rescore_max_period_ext`,
    default 25000) only when needed, recovering large monomers (e.g. TRC_10 →
    16326 bp) without the cost of a global high cap.
- **HOR-order confidence tier** (`hor_order_confidence` in
  `monomer_size_top3_estimats.csv`): `none | strict | supported | weak`,
  distinguishing a recovered founder from a confidently-ordered higher-order
  structure. The HTML report counts "HOR" as `strict ∪ supported`; `weak`
  arrays render as *founder recovered (HOR ≈)* — irregular, very high k, or a
  relaxed rescue. Purely additive (founder values unchanged).
- **Spectral HOR classifier retired** from the report. The unreliable 4-bin
  spectral HOR (strong/moderate/weak/none from the monomer score spectrum) is
  removed everywhere; the raw spectral scores remain as table columns / input to
  detailed analysis. **Report schema bumped to v3** (`hor_order_confidence`
  counts replace the spectral `hor_status` fields).
- **HTML report**: per-TRC distribution ideograms now show array **positions
  only** (single neutral colour, no structural class colour-coding, so a non-HOR
  TRC no longer reads as HOR).
- **Lighter TRC dashboards (satellite-rich genomes)**: the per-array Details
  child row is no longer embedded as HTML per row — each row carries a small key
  and the diagnostics live in one compact embedded per-TRC JSON block, with the
  child row built lazily in JS on unfold (self-contained, works from `file://`).
  A large centromeric-satellite dashboard drops ~4× (e.g. 5.3 MB → 1.3 MB) and
  no longer stalls the browser on expand.
- **Accessible TSV bundle** in `<prefix>_report/data/`, written from the same
  model the HTML renders from (so it mirrors the report): `tra_table.tsv`
  (per array), `trc_table.tsv` (per TRC), `tra_peaks.tsv` (per rescored peak),
  and `columns.tsv` (a data dictionary). A *Download tables (TSV)* link sits
  under each arrays table. `monomer_size_top3_estimats.csv` is unchanged.
- **Unfold column legend**: the per-peak diagnostics table now documents itself
  — column tooltips plus a shared *column meanings* glossary in the expanded row.
- **Report content single-sourced** (maintainability): one `COLUMN_DICT` drives
  every report table's headers, descriptions, the HTML legend, the TSV bundle and
  `columns.tsv`; one `TIER_DEFS` defines the HOR-order tier text once; structural
  blurbs live in one place. A CI check guards legend ≡ dictionary ≡ tier text.

## 1.14.1 (2026-06-10)
- **kitehor 0.13.2** (pinned in `conda-deps.txt`): improved short-monomer
  detection that limits false positives. kitehor 0.13.1's short-monomer
  scanning emitted a pairwise `identity_med` for very short periods (P≈5–20)
  that was high *by chance* on low-complexity arrays (skewed composition →
  frequent matches) with no real periodicity; those cleared TideCluster's
  `identity_med ≥ 0.7` founder gate and collapsed real 180/309 bp satellite
  monomers to 11–14 bp microsatellite harmonics. 0.13.2 fixes this at source —
  the phantom short periods now get a realistic low `identity_med` (≈0.54,
  below the gate), so the founder logic no longer adopts them and needs no
  suppression of its own. Measured on drapa + *S. tuberosum*: founder
  collapses vs the previous baseline fall from 55/68 to 6/11 with **no
  TideCluster-side change**, and those residuals are mostly genuine short
  monomers or short telomere-like families.
- **Short-founder review diagnostics** (`monomer_size_top3_estimats.csv`): new
  `weak_short_founder_flag` (a short founder ≤ 30 bp whose founder peak has weak
  kite support, score < 0.20) and `alt_longer_period` (the dominant credible
  longer-period alternative). These are purely additive — they never change the
  founder call — and surface a surviving weak short founder for manual review.
- **HTML report**: amber `⚠short` badge on flagged arrays (Founder cell +
  Details child row; hover for the longer alternative), and an expanded per-TRA
  legend that now documents the full subrepeat-evidence tier vocabulary and
  every rescored-peak column in the Details child row.
- **Conda recipe pin sync (packaging fix)**: `conda/tidecluster/meta.yaml`
  `requirements.run` had drifted from `conda-deps.txt` — it pinned
  `kitehor=0.12.0` (three minor versions behind) and left `r-igraph` unpinned,
  so `mamba install tidecluster` shipped the wrong KITE backend and an
  unpinned (non-deterministic) igraph. The recipe now mirrors `conda-deps.txt`
  (`kitehor=0.13.2`, `r-igraph=2.0.3`), and `tests/check_conda_recipe_deps.py`
  now fails on a version-pin mismatch (not just a missing package), so future
  drift is caught in CI.

## 1.13.1 (2026-06-05)
- **Deterministic comparative analysis** (`tc_comparative_analysis.R`).
  The MMseqs2-driven satellite-family clustering was not reproducible
  across thread counts or run-to-run at threads > 1. Fixes (refs #4):
  raise `--max-seqs` to `max(10000, n_seqs)` so the prefilter no longer
  drops redundant hits thread-dependently; deterministic best-hit dedup
  (`LC_ALL=C sort -s` by pair then descending pident/qcov/tcov, keep
  first); order `df_pass` by `(query, target)` before building the
  igraph graph so Leiden tie-breaking is invariant to edge order. The
  community table is now identical across thread counts and repeated
  runs. A new `--deterministic` flag additionally forces MMseqs2
  `--threads 1` for a byte-identical `mmseqs2_results.tsv`. Pinned
  `r-igraph=2.0.3` (Leiden tie-breaking is igraph-version sensitive).
  Added `tests/determinism.sh` (`./tests.sh determinism`).

## 1.13.0 (2026-06-05)
- Bumped kitehor to **0.13.0** (`conda-deps.txt`), adopting its two new
  capabilities while TideCluster keeps full control of the founder/HOR
  call. The kite step now runs `rule-classify` then
  `rescore || ssr-scan || tandem-validate` in parallel.
- **Founder go-deeper** (new Pass 6 in `build_monomer_size_csv`):
  consult kitehor's own per-array basic monomer (`hor_basic_period`)
  and adopt it as founder only when it is *deeper* than TideCluster's
  final founder by ≥ 20 %, a clean integer multiple of `strongest`
  (±0.10), and the rescored peak at that period is a near-full-coverage
  array-wide tandem (`coverage_frac ≥ 0.75`, `scan_occupancy_frac ≥
  0.90`). Recovers real deep HORs that the strict ±0.05 divisor gate
  rejects at large `k` while the hard coverage gate rejects kitehor's
  spurious near-2× over-splits; runs last so it never walks a correct
  TRC-consensus founder back to an intermediate divisor. Marked
  `founder_method = "kh_deeper"`. On drapa it recovers 4 confirmed deep
  HORs (e.g. TRC_173 496→165 ×27, TRC_203 3148→1522 ×4).
- **Subrepeats from `tandem-validate`** (kitehor's unified spec-v5
  nested-TR detector) replace the per-peak rescore-column tier
  heuristic for the main subrepeat cell, gated against TideCluster's
  founder: founder-coincident / HOR-subunit candidates are suppressed,
  and only genuine partial-occupancy motifs (density ≥ 0.10, presence
  ≥ 10 %) are surfaced. New `subrepeat_founder_divisor_flag` /
  `subrepeat_founder_divisor_k` columns flag clean high-occupancy
  founder divisors as a possible-founder-miscall watch-list. The legacy
  `_classify_subrepeat_tier` path remains as a fallback for run
  directories without a tandem-validate file.

## 1.12.2 (2026-06-04)
- Pass 1 now generates founder candidates from peak **clusters** rather
  than individual peaks. Within each cluster (single-link `±5%/±100 bp`,
  same `_cluster_peaks_by_period` Pass 5 uses), the score-weighted mean
  of cluster members that pass the founder identity gate
  (`_FOUNDER_ID_MIN = 0.7`) is tested as a divisor of `strongest_period`
  under the existing strict `|k − round(k)| ≤ 0.05` gate. This recovers
  arrays where the basic monomer is supported by multiple nearby peaks
  but no single peak forms a clean integer divisor on its own — e.g.
  drapa TRC_26 `chr9:1816989` with strongest 1880 and peaks
  154/158/161/165 each at `id_med ≈ 0.95` but `|k − kr|` up to 0.21
  individually; the score-weighted mean ≈ 156.8 gives `k = 11.99`
  (clean) → founder 157 / multiplicity 12 (previously 315 / ×6, the
  multiple-of-basic Pass 1 settled on as the only single-peak clean
  divisor).
- Defensive: when a cluster's mean fails the strict gate but an
  individual member of that cluster passes it (e.g. drapa TRC_45
  `chr12:19629243` cluster {334, 337, 340, 344} where the mean has
  `|k − kr| = 0.32` but member 344 alone has `|k − kr| = 0.035`), Pass 1
  falls back to per-member candidate generation inside that cluster —
  the pre-existing single-peak path. Each cluster contributes either
  its mean (when the mean is clean) or its qualifying single members
  (when the mean is not). The per-kr-winner sort from 9254c93 then
  arbitrates across all candidates.
- Pass 2 expansion (Lever 3): the
  `if e["multiplicity"] > 1 or e["fallback"]: continue` gate is widened
  to also re-evaluate Pass-1-success rows whose founder sits at an
  integer multiple of the TRC consensus and is significantly off it
  (`founder / consensus ∈ [2..5]` with `|ratio − round(ratio)| ≤ 0.10`
  AND `|fp − consensus| / consensus > 0.20`). Those rows hand off to
  the existing `_rescue_founder_from_trc` helper, which finds a near-
  consensus peak with the relaxed `id_med ≥ 0.5` gate and sets
  `irregular_multiplicity = true` because the resulting `k` is
  fractional. Example: drapa TRC_138 `chr4:26085301` rescues from
  founder 262 (= 2×131, mult 13) to founder 131 (mult 26 irregular,
  raw k 25.92) — the basic monomer the TRC's consensus says is right.
- TRC_20-class arrays (TRCs with fewer than `_TRC_CONSENSUS_MIN_ARRAYS
  = 3` Pass-1 successes) still cannot be rescued because the TRC
  consensus is not computed at all — known limitation, applies to
  tiny / fragmented TRCs only.
- Validation on `test_data/drapa` (2073 arrays):
    - 0 strongest_period changes (Lever 2/3 leave strongest untouched
      by design).
    - 613 founder shifts (mostly within-monomer-family micro-shifts
      from cluster mean vs single peak — cosmetic; same basic
      monomer).
    - 38 multiplicity changes: 27 ↑ recoveries (×6→×12, ×7→×14,
      ×5→×10, ×4→×8 etc.) and ~11 ↓ adjustments where the cluster
      mean differs slightly from the single-peak pick (same
      decomposition depth, founder bp nudged a few % within the
      monomer family).
    - User's cited case `chr9:1816989`: founder 315 → 157,
      mult ×6 → ×12, `mult_raw 11.987` (clean integer).
    - `TRC_138 chr4:26085301` rescued by Lever 3: founder 262 → 131,
      mult ×13 → ×26 (`irregular_multiplicity=true`,
      `mult_raw=25.916`, `founder_method=pass2`).
- Validation on `test_data/Solanum_pimpinellifolium` (433 arrays):
    - Flagship 1.12.1 cases (TRC_8 chr12:1314680 / 1327272) unchanged
      — founder 178, mult ×10.
    - 86 founder shifts, 28 multiplicity changes (22 recoveries, 6
      minor down-shifts), 0 strongest changes.
- Two new unit cases in `tests/test_strongest_by_identity.py`:
    - `(g)` Lever 2 — cluster-mean founder lands in basic-monomer band
      with multiplicity recovered.
    - `(h)` Lever 3 — multi-array TRC where a single mult>1 stuck row
      gets rescued to the TRC consensus founder.
- Design notes in `docs/founder_cluster_mean_plan.md`.

## 1.12.1 (2026-06-04)
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

- Fix within-kr founder tiebreaker in `tc_utils.build_monomer_size_csv()`
  Pass 1. After the strongest-by-identity change landed (944735d), several
  rows whose strongest had multiple valid divisors at the *same* `kr =
  round(strongest/P)` level were still picking the wrong representative
  — the old "smallest period wins" sort favoured numerically smaller P
  over the cleanest integer relationship. On drapa TRC_2 chr3:56557:
  strongest=2176 had three kr=2 candidates {1066 (k=2.041, id_med
  0.972), 1088 (k=2.000, id_med 0.995), 1091 (k=1.995, id_med 0.995)};
  1066 won purely because numerically smallest, despite being the
  fuzziest divisor with the lowest identity. The fix groups candidates
  by `kr` and within each group picks the best representative (highest
  id_med → cleanest `|k − kr|` → smallest P), then across kr groups
  the smallest period still wins (preserves the basic-monomer
  semantic so a kr=10 divisor outranks a kr=2 divisor when both
  qualify). On drapa: 490 of 2073 rows now report a slightly better
  founder at the same multiplicity (0 strongest changes, 0
  multiplicity changes). Covered by two new cases in
  `tests/test_strongest_by_identity.py`.

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
