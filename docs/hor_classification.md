# HOR classification (KITE)

## 1. Background

A **higher-order repeat** (HOR) is a tandem repeat in which a block of
several basic monomers has itself been duplicated. The characteristic
signature on a raw sequence level is that the array shows periodicity
at the basic monomer length *m* **and** at an integer multiple
*n · m* corresponding to the HOR period. In centromeric satellites
HOR structures with *n* up to 30 or more are common.

TideCluster's KITE module (`tarean/kite.R`) detects these periods by
computing the distribution of distances between every pair of
repetitive *k*-mers in each tandem repeat array (TRA) and extracting
significant peaks. Two or more peaks that sit at integer multiples of
a common base are evidence of HOR structure; their relative heights
and the tightness of the integer-multiple relationship determine how
strong that evidence is.

## 2. Inputs from KITE

For each TRA the KITE pipeline returns up to three candidate monomer
sizes, ranked by score:

| Symbol | Meaning |
|---|---|
| *m₁, m₂, m₃* | monomer-size estimates in bp (top three peaks) |
| *s₁, s₂, s₃* | weighted *k*-mer-interval scores for each peak |

Peaks with no supporting evidence return `NA` for both the position
and the score. Any combination of the three positions can be missing.

## 3. Algorithm

The classifier searches for the base monomer *m\** that best explains
the observed peaks as a harmonic series, then scores how well the
harmonic structure is actually supported by the data.

### 3.1 Candidate search

Candidate *m\** values are drawn from

- every observed peak *mᵢ*;
- every *mᵢ / k* for *k = 2 … HOR_MAX_K*; and
- a fine 1 bp grid between `min(m) / 2` and `max(m) / 2`.

The fine grid is essential — it lets us pick a real-valued base that
fits all three peaks to within a few percent even when the observed
peaks drift by several percent due to natural monomer variability
(e.g. `(1692, 925, 2617)` fits `(2, 1, 3) × 873`).

### 3.2 Per-candidate scoring

For each candidate *m\**:

    kᵢ      = round(mᵢ / m*)                         integer multiple
    errᵢ    = |mᵢ − kᵢ · m*| / m*                     fractional error
    closeᵢ  = max(0, 1 − errᵢ / HOR_TOL)              1 at exact, 0 at tol
    skip unless ∃ kᵢ = 1  AND  ∃ kᵢ ≥ 2
    f_base  = Σ sᵢ · closeᵢ   over i with kᵢ = 1   / Σ sᵢ
    f_harm  = Σ sᵢ · closeᵢ   over i with kᵢ ≥ 2  / Σ sᵢ
    support = √(f_base · f_harm)
    bonus   = 1 + HOR_HARMONIC_BONUS · max(0, n_distinct_harmonics − 1)
    conf    = support · bonus

where *n_distinct_harmonics* is the number of distinct *kᵢ* ≥ 2 that
still had `closeᵢ > 0`. The candidate with the highest `conf` wins.

### 3.3 Why a geometric mean

Using `√(f_base · f_harm)` rather than a simple sum or product has
two consequences that matter in practice:

- If the base peak has no real score (e.g. it sits at *k = 1* of a
  noise peak while the "harmonic" is actually the primary signal),
  `f_base → 0` and the support collapses to 0. This blocks the noise-
  peak false positives that defeated the previous hard-rule
  `HOR-dominant` classifier.
- If the harmonic peak has a score much smaller than the base (a
  "weak but real" HOR-visible pattern), `f_harm` is small but not
  zero, and the support stays positive in a way that is smooth in
  score rather than gated on a fixed ratio like the previous
  `s₂ ≥ 0.8 · s₁` rule — which, as a detectability calculation
  shows, precludes any real HOR with *n ≥ 2*.

### 3.4 Category bins

The continuous `conf` is the primary output; a categorical label is
derived from it for display:

| confidence | label |
|---|---|
| `< HOR_BIN_WEAK`                        | `No HOR`        |
| `[HOR_BIN_WEAK, HOR_BIN_MODERATE)`      | `HOR weak`      |
| `[HOR_BIN_MODERATE, HOR_BIN_STRONG)`    | `HOR moderate`  |
| `≥ HOR_BIN_STRONG`                      | `HOR strong`    |

The thresholds were picked so that the three Drapa TRC_1 reference
arrays (§7) land at the biologically intuitive labels. They are not
calibrated against a statistical null; if future analyses show the
distribution shifting, re-tune the constants rather than invent a
data-driven quantile scheme (genome-to-genome distributions differ
too much for one-size-fits-all quantiles to be meaningful).

### 3.5 Reported values

For each TRA the classifier reports:

| Field | Meaning |
|---|---|
| `HOR_status` | categorical label (`No HOR` / `HOR weak` / `HOR moderate` / `HOR strong`) |
| `HOR_confidence` | continuous score in [0, ~1] |
| `HOR_base_monomer` | fitted *m\** rounded to bp |
| `HOR_hor_period` | largest supported harmonic, i.e. *k_max · m\**, rounded to bp |
| `HOR_n_harmonics` | number of distinct *k ≥ 2* that contributed to the support |

If no candidate *m\** satisfied the `base + harmonic` requirement
(for instance when *m₂* and *m₃* are both missing or when the array
only supports a single dominant monomer period), all five HOR fields
are reported as `NA`/`0` appropriately and the status is `No HOR`.

## 4. Per-TRC rollup

The KITE report aggregates array-level calls per TRC into four
counts:

    N_no_HOR, N_HOR_weak, N_HOR_moderate, N_HOR_strong

plus a median confidence across the TRC's arrays. **No single TRC-
level HOR label is assigned** — a TRC can easily contain arrays with
different HOR structures, and a single label would be reductive.

## 5. Outputs

### 5.1 TSV schema (`<prefix>_kite/monomer_size_top3_estimats.csv`)

```
TRC_ID · seqid · start · end
monomer_size · score · array_length
monomer_size_2 · score_2
monomer_size_3 · score_3
HOR_status · HOR_confidence · HOR_base_monomer · HOR_hor_period · HOR_n_harmonics
```

### 5.2 Report v2

- KITE tab: aggregate 4-segment bar chart across all arrays + per-TRC
  table with four tinted count columns and median-confidence column.
- Merged TAREAN tab: four-count HOR cell per TRC.
- Per-TRC dashboard arrays table: `m₁ · s₁ · m₂ · s₂ · m₃ · s₃ ·
  HOR status · Confidence · Base (bp) · HOR period (bp)`.
  `HOR_n_harmonics` is exposed in the TSV only, not in the dashboard
  (keeps the row narrow enough for readability).

## 6. Tuning

Constants live at the top of `tarean/kite.R`:

| Constant | Default | Effect |
|---|---|---|
| `HOR_TOL` | 0.10 | Fractional tolerance in the integer-multiple fit. Increase to accept more monomer-length variability at the cost of false positives. |
| `HOR_HARMONIC_BONUS` | 0.5 | Weight given per extra distinct harmonic. Increase to reward richer harmonic series. |
| `HOR_BIN_WEAK / MODERATE / STRONG` | 0.10 / 0.20 / 0.40 | Categorical bin thresholds on the continuous confidence. |
| `HOR_MAX_K` | 5 | Highest denominator when generating candidate *m\** = *mᵢ / k*. Increase if you expect very-high-order HORs. |
| `HOR_GRID_STEP` | 1 | bp resolution of the fine candidate grid. Keep at 1 unless you need speed on very long TRAs. |

## 7. Worked examples (Drapa TRC_1)

### #18 — noise-peak false positive in the previous rule

Input: `m = (1692, 66, 401)`, `s = (0.533, 0.004, 0.004)`.

Old rule: `m₁ / m₂ = 25.6 ≈ 26`, `m₁ / m₃ = 4.22 ≈ 4`, both within
10 %, both ≥ 2 → **HOR-dominant, n = 26**. The secondary peaks are
noise (score ≈ 1 / 130 of the primary), so this is wrong.

New rule: best *m\** found at 65 bp, fits `k = (26, 1, 6)` with
small errors. `f_base = 0.85 · 0.004 / 0.541 = 0.006` — the base
peak has effectively no score. `f_harm ≈ 0.68`. Support =
`√(0.006 · 0.68) = 0.064`. No bonus, confidence ≈ 0.066.
**Label: `No HOR`.** ✓

### #124 — approximate HOR missed by the old rule

Input: `m = (1692, 925, 2617)`, `s = (0.223, 0.208, 0.148)`.

Old rule: `m₁ / m₂ = 1.83` — not within 10 % of 2. → **No HOR**.
But visual inspection suggests that 925 is a candidate base monomer
and 1692, 2617 are 2× and 3× harmonics (with the underlying monomer
slightly shorter than 925 due to monomer variability).

New rule: best *m\** found at 873 bp, fits `k = (2, 1, 3)` with
errors ≈ 6 %, 6 %, 0.2 %. `f_base = 0.40 · 0.208 / 0.579 = 0.144`,
`f_harm = 0.38 · 0.223 / 0.579 + 0.98 · 0.148 / 0.579 = 0.397`.
Support = `√(0.144 · 0.397) = 0.239`. Two distinct harmonics (2, 3),
bonus = 1.5. Confidence ≈ 0.360. **Label: `HOR moderate`.** ✓

### #119 — weak but real 2× HOR

Input: `m = (925, 1692, 719)`, `s = (0.522, 0.063, 0.005)`.

Old rule: `m₂ / m₁ = 1.83` — within 10 % of 2, but
`s₂ = 0.063 ≪ 0.8 · s₁ = 0.418`, so the visible rule fails.
→ **No HOR**. But visually a weak 2× HOR is present.

New rule: best *m\** found at 865 bp, fits `k = (1, 2, 1)` with
errors ≈ 7 %, 4 %, 17 %. `f_base = 0.31 · 0.522 / 0.590 = 0.275`,
`f_harm = 0.56 · 0.063 / 0.590 = 0.060`. Support =
`√(0.275 · 0.060) = 0.128`. No bonus. Confidence ≈ 0.127.
**Label: `HOR weak`.** ✓

## 8. Known limitations

- **Tolerance grows with *n***. Because `closeᵢ` uses
  `errᵢ / HOR_TOL`, and `errᵢ` is fractional (already scaled by
  *m\**), the absolute bp tolerance widens with higher *k*. For *n*
  up to ~30 this is biologically reasonable (measurement noise is
  relative); for extreme *n* values it may allow harmonic fits that
  are not real.
- **Base + harmonic both required**. Arrays that only produce a
  single strong peak (*m₂*, *m₃* both missing) can never be labelled
  HOR — the confidence formula needs both a *k = 1* and a *k ≥ 2*
  peak.
- **1 bp grid resolution**. The fine grid searches at 1 bp steps,
  not sub-bp. In practice monomer lengths are integer bp so this is
  fine, but the fitted *m\** loses any fractional-bp accuracy that
  could in principle be recovered by least-squares.
- **Harmonic independence is assumed but not verified**. A strong
  *m\** naturally produces harmonics at 2 *m\**, 3 *m\**, …; the
  classifier treats their co-occurrence as confirming evidence, which
  is correct for a noisy real HOR but also fires on pure harmonic
  echoes of a single underlying period.

## 9. Previous design

The previous classifier used a two-category hard rule
(`HOR-dominant` / `HOR-visible` / `No HOR detected`) with fixed
thresholds:

- `HOR-dominant`: `m₁ / m₂ ≈ integer ≥ 2` within 10 % **and**
  `m₁ / m₃ ≈ integer ≥ 2`.
- `HOR-visible`: `m₂ / m₁ ≈ integer ≥ 2` within 10 % **and**
  `s₂ ≥ 0.8 · s₁`.

Three failure modes motivated the rework:

1. **Noise-peak false positives in HOR-dominant.** The rule checked
   positional ratios only, so an array with a dominant primary peak
   and two tiny noise peaks at integer divisors was called HOR-
   dominant (Drapa TRC_1 #18).
2. **Unreachable HOR-visible score ratio.** For a real HOR with
   *n ≥ 2* the secondary peak score is approximately `1 / n` of the
   primary, so `s₂ ≥ 0.8 · s₁` effectively required `n ≤ 1.25` —
   excluding every real HOR from the HOR-visible category (Drapa
   TRC_1 #119).
3. **Approximate multiples discarded.** A strict 10 % tolerance on
   pairs of observed peaks could not accommodate monomer-length
   variability — if the base monomer drifted and the harmonic did
   not drift in the same direction, no pair of peaks sat on the
   integer grid even when the underlying structure was clear
   (Drapa TRC_1 #124).

The confidence-score design described above addresses all three by
(a) weighting by real score mass rather than ratios alone,
(b) smoothing the asymmetry between base and harmonic scores, and
(c) fitting *m\** to the whole peak set rather than insisting on
pairwise integer fits of observed peaks.
