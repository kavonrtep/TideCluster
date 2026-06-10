# KITE/HOR Founder Review

Date: 2026-06-10

Scope: read-only review of TRA characterization, with focus on `kitehor`
integration, founder selection, HOR multiplicity calls, SSR handling, and
miscalled founder risk. No code was changed.

Primary code paths reviewed:

- `TideCluster.py` KITE step: `kite-periodicity`, `rule-classify`, `rescore`,
  `ssr-scan`, `tandem-validate`, selective long-period rescore, and
  `build_monomer_size_csv`.
- `tc_utils.py`: founder/strongest/multiplicity derivation, TRC-consensus
  rescue, SSR-family override, fallback cluster rescue, `kh_deeper` adoption,
  tandem-validate subrepeat gating, and report-facing diagnostics.
- `tc_rerender_report.py`: how multiplicity, founder method, fallback, SSR, and
  subrepeat flags are interpreted in the HTML report.
- Example output:
  `test_data/Solanum_lycopersicum/tc_kite/monomer_size_top3_estimats.csv`
  and raw `kitehor.*.tsv` files.

Commands run:

```bash
conda run -n tidecluster=1.14.1 python3 tests/test_strongest_by_identity.py
conda run -n tidecluster=1.14.1 python3 tests/test_cluster_rescue.py
conda run -n tidecluster=1.14.1 python3 -c "import tc_utils as tc; base='test_data/Solanum_lycopersicum'; tc.build_monomer_size_csv(kite_tsv=base+'/tc_kite/kitehor.kite.tsv', ssr_tsv=base+'/tc_kite/kitehor.ssr.tsv', rescored_peaks_tsv=base+'/tc_kite/kitehor.rescored.peaks.tsv', tandem_validate_tsv=base+'/tc_kite/kitehor.tandem_validate.tsv', out_csv='/tmp/slyco_current_monomer_size_top3.tsv', trc_repeat_type=tc.parse_trc_ssr_motif_len(base+'/tc_clustering.gff3'))"
```

Both focused founder/HOR regression tests passed.

## Current Solanum Status

Regenerating `monomer_size_top3_estimats.csv` from the existing raw Solanum
`kitehor` TSVs with current code produced 748 rows, matching all 748 clustering
arrays by `TRC_ID/seqid/start/end`. There is no row-loss issue in this example.

The regenerated table is not byte-identical to the checked-in Solanum table:
66 founder periods and 61 multiplicities differ. Most differences are explained
by current founder-core changes that deepen `TRC_4` intermediate founders
such as 315, 477, or 1128 bp to the prevalent 53 bp founder, plus a few
SSR-family override changes.

Several previously documented Solanum issues are now fixed by current code:

- `TRC_2 chr2:15535707-15614778`: current output keeps a long founder around
  9520 bp instead of incorrectly forcing the ATC 3 bp SSR motif.
- `TRC_18`: current output forces the family founder to 3 bp via the clustering
  `repeat_type=SSR` motif.
- `TRC_5`: current output gives a uniform 3 bp SSR founder and avoids the old
  confusing fallback marker on the cited array.
- `TRC_4 chr6:9719494-9750762`: current output deepens the founder to 53 bp.

The old `TRC_10 chr10:59071587-59183483` above-cap case is still fallback when
only the existing raw TSVs are reprocessed, because the raw rescore table has
the 16326 bp period with `identity_med=NA`. The selective long-period extension
in the current pipeline should fix this on a fresh run, but the summary table
could make the above-cap condition more explicit.

## Main Risks

### 1. Very high multiplicities are counted as ordinary HOR calls

In the regenerated Solanum table, 134 arrays have `multiplicity > 30`; 130 are
`TRC_4`, and the maximum observed call is `x280`. Many of these probably reflect
real founder deepening to the 53 bp basic monomer, but the report treats every
`multiplicity > 1` as an HOR array.

That conflates two ideas:

- founder recovery: "the best basic monomer is 53 bp";
- confident HOR order: "this array has a supported xN higher-order structure".

For `pass2`/consensus-rescued or irregular high-k calls, the second statement is
weaker than the first. The report and downstream summaries should distinguish
these.

Recommended mitigation:

- Add a high-k or irregular confidence tier, for example
  `hor_order_confidence = strict | consensus_rescued | high_k | irregular`.
- Count `multiplicity > 1` as "HOR" only when the supporting path is strict or
  otherwise strongly supported.
- Surface high-k rescue counts separately in the report.

### 2. TRC-consensus rescue is useful but under-gated

Pass 2 looks for a peak near the TRC consensus founder with relaxed
`id_med >= 0.5`, then sets `multiplicity = round(strongest/founder)`.
It does not currently require sufficient kite score, coverage, scan occupancy,
or a maximum supported k for the final HOR-order interpretation.

This is useful for rescuing real basic monomers, but it can inflate high-order
HOR calls when `strongest` is a long/noisy period and the consensus founder is
only weakly supported in that particular array.

Recommended mitigation:

- Keep consensus founder rescue, but record a distinct founder method and
  confidence tier.
- For reporting an HOR order, add support gates using available columns:
  `score`, `coverage_frac`, `scan_occupancy_frac`, `identity_n`,
  `delta_id_pp`, and `multiplicity_raw` distance from integer.
- Treat high-k rescue as founder evidence first, HOR-order evidence second.

### 3. Period clustering can merge short-period harmonics

`_cluster_peaks_by_period()` uses a mixed window:

```text
max(member_period * 0.05, 100 bp)
```

For short founders this 100 bp minimum is very large. Peaks at 53, 106, and
159 bp can be merged into one single-link cluster even though they are distinct
harmonics. This can hide alternative periodicities in `alt_cluster_*` and can
affect fallback rescue interpretation.

Recommended mitigation:

- Use a smaller minimum window for short periods.
- Or make cluster windows period-dependent, for example a few bp below 200 bp
  and percentage-based for larger periods.
- Or use harmonic-aware grouping: nearby-period clusters and harmonic-series
  clusters should be separate concepts.

### 4. SSR-family override is biologically useful but coarse

Current code parses `repeat_type=SSR` and the `ssr` motif from the clustering
GFF3, then forces every array in that TRC to the motif length. This fixes the
known `TRC_18` and `TRC_5` short-founder problems and avoids the old per-array
SSR coverage bug that incorrectly clobbered `TRC_2`.

However, in Solanum some arrays inside an SSR TRC have low per-array SSR scan
coverage, while still receiving the family SSR founder. This may be correct,
because clustering consensus is the authoritative family decision, but it is a
discordance that should be visible.

Recommended mitigation:

- Keep clustering `repeat_type=SSR` as authoritative for founder selection.
- Add/report a discordance flag when `ssr_founder_override=true` but per-array
  `ssr_flag=no` or per-array SSR coverage is low.
- Continue showing kitehor SSR scan as annotation, not as a founder override
  for `repeat_type=TR` families.

### 5. Above-cap long periods need clearer diagnostics

The selective long-period extension is present in current code, but a stale
or partially rerendered output can still show fallback founder calls when a real
long period was above the original rescore cap. `TRC_10` is the example:
rank 2 at 16326 bp has `identity_med=NA`, while the fallback 35 bp peak has low
identity and low coverage.

Recommended mitigation:

- Add an output flag such as `above_rescore_cap_candidate=true`.
- Report the best above-cap period and score even when it was not rescored.
- On fresh pipeline runs, validate that `--long` plus selective extension
  rescues `TRC_10` to the 16326 bp period.

### 6. Report color semantics can overstate structure

The report treats `multiplicity > 1` as HOR for array coloring and counts.
For non-HOR or ambiguous TRCs, this can make the genome-distribution view look
like a structural classifier when the user's preferred interpretation is simply
TRA/TRC position.

Recommended mitigation:

- Render per-TRC distribution ideograms with neutral position-only coloring.
- Keep structural signals in the array table/details, where founder method and
  confidence can be shown.

### 7. Documentation drift

There is stale documentation around the executable KITE path. For example,
comments in `TideCluster.py` still say the `rule-classify + tandem-validate`
cascade is dropped, but the current code runs both and consumes
`tandem_validate.tsv`.

Recommended mitigation:

- Update comments/docstrings near the KITE step and `build_monomer_size_csv`.
- Document the actual current pass order:
  strict per-array selection, consensus rescue, prevalent-founder anchor,
  SSR-family override, cluster fallback rescue, and `kh_deeper` adoption.

## Recommended Priority

1. Split "founder recovered" from "confident HOR order" in output and report
   semantics, especially for high-k/irregular/consensus-rescued calls.
2. Tighten or annotate Pass 2 consensus rescue with support evidence.
3. Fix short-period clustering windows so harmonics are not merged as one
   alternative-period cluster.
4. Add diagnostics for SSR-family/per-array discordance and above-cap long
   period candidates.
5. Neutralize per-TRC ideogram coloring if the view is intended to show
   positions only.
6. Update stale KITE integration comments/docs.

## Validation Suggestions

For future changes, validate at three levels:

- Unit fixtures: keep `tests/test_strongest_by_identity.py` and
  `tests/test_cluster_rescue.py` green, and add targeted fixtures for high-k
  consensus rescue and short-period clustering.
- Solanum targeted rows: check `TRC_2`, `TRC_4`, `TRC_5`, `TRC_10`, `TRC_18`.
- Cross-genome regression: rerun on drapa, S. tuberosum, S. lycopersicum, and
  a CEN-style short fixture to ensure short founders and large monomers do not
  regress.
