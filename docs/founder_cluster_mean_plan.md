# Pass 1 founder picking: cluster-mean + Pass 2 expansion

Branch: `feat/founder-cluster-mean` (off `main` at 2647f6e)

## Motivation

Rescore's per-array peaks are spread by bp-level imprecision. On
arrays with a clear basic monomer (e.g. drapa TRC_26 ~154 bp) the
peak list often shows several nearby periods (154, 158, 161, …) each
with high `id_med`, none of which forms a clean integer divisor of
the array's `strongest_period` by itself, but whose score-weighted
mean does. The current strict `|k − round(k)| ≤ 0.05` divisor gate
rejects all of them and Pass 1 falls back to an integer multiple of
the basic (e.g. founder = 315 = 2·154 for `chr9:1816989`, mult = 6
instead of the biologically true mult = 12).

Empirically on drapa: of 20 "suspicious" arrays (founder ≠ TRC
dominant founder), **17 of 20 (85 %)** have a score-weighted mean
period across nearby peaks that gives a CLEANER divisor of strongest
than any single peak. **15 of 20 (75 %)** become clean (`|k − kr| ≤
0.05`) under the mean. The remaining 5: 2 have only one peak in the
band (so clustering can't help), 1 has a worse mean than the best
single peak, and 2 are improved but still slightly fuzzy.

## Lever 2: peak-cluster mean as founder candidate (Pass 1)

Where the change happens: `tc_utils.build_monomer_size_csv()` Pass 1,
the divisor search currently at lines ~2593-2630 (after the
strongest-by-identity block).

### Algorithm

After `strongest_period` is fixed (existing argmax(id_med) logic
from 944735d / 9254c93):

1. **Cluster the array's rescored peaks by period** using the
   existing `_cluster_peaks_by_period` helper (single-link, mixed
   `max(±5%, ±100bp)` window). Each cluster carries `members`,
   `median_period`, `score_sum`, `cov_frac_max`, etc.

2. **For each cluster, compute a representative period**
   = score-weighted mean of member periods. Helper:

   ```python
   def _cluster_weighted_mean(cluster):
       members = cluster["members"]
       num = 0.0; den = 0.0
       for m in members:
           p  = _num(m.get("period"))
           sc = _num(m.get("score"))
           if p is None or sc is None or sc <= 0: continue
           num += p * sc
           den += sc
       return num / den if den > 0 else None
   ```

3. **Filter to clusters that qualify as founder candidates**:
   - At least one cluster member passes `id_med ≥ _FOUNDER_ID_MIN`
   - cluster mean `P_mean` satisfies `2 ≤ round(strongest/P_mean) ≤ _KMAX`
     and `|k − round(k)| ≤ _RATIO_TOL` (strict ±0.05 — no relaxation)
   - cluster mean < strongest (only divisors)

4. **Pick the winning cluster**:
   - Across clusters: smallest cluster mean wins (basic-monomer rule —
     preserves "kr=10 beats kr=2 when both qualify").
   - Within-kr ties: highest cluster `id_med_max`, then cleanest
     `|k − round(k)|`, then smallest mean (mirrors the within-kr
     tiebreaker from 9254c93 but applied to clusters).

5. **Set founder**:
   - `founder_period = round(cluster.weighted_mean)`
   - `founder_peak = cluster.member with highest id_med` (so
     `founder_id_med`, `delta_id_pp`, etc. read off the cluster's
     best representative — matches what `founder_id_med` semantically
     means).
   - `multiplicity = round(strongest / founder_period)`
   - `multiplicity_raw = strongest / founder_period` (the actual ratio
     of the rounded founder; preserves the existing semantic).
   - `irregular_multiplicity` stays false at this point — the
     clustering already absorbed the bp-level noise so the resulting
     k is integer-clean by construction.

### Compatibility with existing single-peak behaviour

A cluster of size 1 has `weighted_mean = the peak's period` exactly,
so on arrays with a single clean basic-monomer peak this algorithm
produces the same founder as today. The new behaviour only differs
on arrays where multiple peaks cluster.

### Effect on previously-correct rows

For each of the 2073 drapa rows, the cluster mean of the
(single-peak) winning cluster equals the peak's period itself, so
**founder values for rows already at the basic monomer are
unchanged**. Only the ~20 "suspicious" rows are expected to shift.

### Inline note in the docstring

The Pass-1 section gets a new bullet describing the cluster-mean
candidate generation.

## Lever 3: Pass 2 expansion for mult>1 rows far from TRC consensus

The 2 single-peak-in-band cases (TRC_138 `chr4:26085301`, TRC_20
`chr10:11263647`) can't be helped by Lever 2 — there's only one peak
near the dominant founder, so the cluster is a singleton and the
mean equals that peak's period (which alone doesn't divide cleanly).

Today Pass 2 gates on `if e["multiplicity"] > 1 or e["fallback"]:
continue`. Replace the gate with an additional escape clause: also
re-evaluate rows where:
- `multiplicity > 1`, and
- `founder_period` is significantly far from `trc_consensus_founder`
  (e.g., `|founder − consensus| / consensus > 0.20`), and
- the founder is approximately an integer multiple of consensus
  (`round(founder / consensus) ∈ [2, 5]`, `|ratio − round(ratio)| ≤
  0.10`) — the signature of "Pass 1 settled on a multiple of the
  basic".

For those rows: look for a rescored peak within `max(±10%, ±20bp)`
of `trc_consensus`, with `id_med ≥ _RESCUE_ID_MIN (0.5)`, and accept
it as the new founder even if the resulting multiplicity is
fractional (`irregular_multiplicity = true`, raw k in
`multiplicity_raw`).

This re-uses the existing `_rescue_founder_from_trc` helper; only
the entry gate needs to widen.

### Why two layers

Lever 2 handles arrays whose own peak list contains the basic
monomer's full evidence (the within-array clustering recovers it).
Lever 3 handles arrays where only ONE peak survived near the basic
(the TRC consensus across other arrays is the only oracle that can
say "we know the basic should be ~154 even if your single peak is at
154 with no neighbours"). Both layers are needed for full coverage
of the 20 drapa cases; together they should fix all 20.

## Provenance / diagnostic columns

No new CSV columns required. Existing columns capture provenance:
- `irregular_multiplicity = true` when Pass 2/3 fired (Lever 3 path)
- `founder_method` enum extended? Maybe add a new value `cluster_mean`
  for Lever 2's path, distinguishing it from `strict`. Tentative — to
  decide during implementation.
- `kitehor_founder_period` unchanged (still kitehor's original pick)

## Test strategy

Synthetic peaks in `tests/test_strongest_by_identity.py` (extending
existing file rather than a new file, since the logic lives in the
same Pass-1 region):

1. **chr9:1816989 mirror** — strongest=1880, peaks at 154/158/161/165
   with the actual rescore scores. Expected: founder ≈ 157
   (score-weighted mean), multiplicity = 12.

2. **Single peak in band** (Lever 3 territory) — strongest=1860,
   only one peak in the basic band at 154 (k=12.08, fuzzy). TRC
   consensus = 154. Pass 1 falls to multiple (308, mult=6); Pass 2
   expansion rescues to founder=154, mult=12, irregular_multiplicity=true.

3. **Regression** — strongest=2176 with kr=2 candidates 1066/1088/1091
   (the existing test (e)). Cluster mean should still pick a clean
   ~1088 representative. Confirm founder=1088 mult=2 unchanged.

4. **No-cluster regression** — single clean peak at 178 with kr=10
   relative to strongest=1786. Cluster of size 1 → mean=178 →
   founder=178, mult=10. Confirms single-peak case isn't broken.

## Validation on drapa

After implementation, regenerate
`monomer_size_top3_estimats.csv` on the existing kitehor outputs.
Compare against the current 1.12.1 output:

- **Required**: all 20 previously-suspicious rows resolve to a
  founder within ±5% of the TRC dominant founder, with mult equal
  to round(strongest / dominant_founder).
- **Required**: 0 regressions on rows whose founder was already at
  the TRC dominant.
- **Acceptable**: some "well-behaved" rows may shift founder by ±5%
  (within-monomer-family noise) because the cluster mean is slightly
  different from the single highest-score peak. These don't change
  the biological interpretation (same monomer family).
- **Disallowed**: structural shifts on previously-correct rows.

Compare per-row HOR multiplicity distribution:
- before: 1016 HOR calls
- after: expected to remain ~1016 (Lever 2/3 changes founder bp but
  not multiplicity on already-correct rows; on suspicious rows
  multiplicity goes UP, e.g. mult=6 → mult=12, which doesn't
  change the HOR count but does affect the per-TRC ×k breakdown).

## Files touched

- `tc_utils.py`: `build_monomer_size_csv` Pass 1 + Pass 2 gate
- `tests/test_strongest_by_identity.py`: 2-3 new cases
- `changelog.md`: Unreleased entry
- `docs/founder_cluster_mean_plan.md`: this doc (delete or keep as
  rationale archive — I'll keep it)

## Out of scope

- TAREAN's per-TRC `monomer_length` as a hard oracle (requires
  pipeline reorder; explicitly deferred per the data — within-array
  clustering plus Pass 2/3 TRC consensus covers what TAREAN's
  monomer would add anyway).
- Tolerance relaxation by identity (Lever 1 from the original
  three-option analysis) — superseded by the cluster-mean approach.

## Rollback plan

If validation finds an unacceptable number of regressions, revert
the branch and we discuss design adjustments. The CSV column schema
isn't extended, so no consumer migration is needed in either
direction.
