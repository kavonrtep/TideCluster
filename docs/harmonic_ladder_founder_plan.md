# Implementation plan — harmonic-ladder founder (Path 1)

Date: 2026-06-11
Status: IMPLEMENTED (recommended gates + exceptionally-clean ×2). Pass 7 in
tc_utils.build_monomer_size_csv. Blast radius (tools/founder_diff.py): exactly 2
arrays across slyco/arabidopsis/drapa — slyco TRC_47 (179 ×6) and drapa TRC_158
(111 ×3), both legitimate ladders (3× / 3× units conserved at 0.91 / 0.9955);
arabidopsis 0 changes, all confirmed satellites untouched. Fixture + unit test +
suite green.
Source: TRC_47 analysis (this session). Target release: 1.15.1 / next.
Scope: `tc_utils.build_monomer_size_csv` (founder logic) + a fixture/test +
before/after comparison. No report/schema change beyond a new `founder_method`.

## Goal

Recover the **basic monomer** of a divergent HOR satellite when it sits below the
identity gate but is corroborated by a clean **harmonic ladder**. TRC_47
`chr5:68075763` is the motivating case: the basic 179 bp monomer has
`id_med 0.687` (just under `_FOUNDER_ID_MIN = 0.70`), so the current logic settles
on the conserved 6×179 ≈ 1073 HOR unit. KITE, TideHunter, and TAREAN all rank
1073/1079; only the **ladder** evidence distinguishes 179 as the true fundamental.

### The signal (measured on TRC_47 — distinct rungs at integer multiples of 179)

| rung m | period | id_med | id_iqr | scan_occ |
|---|---|---|---|---|
| **1×** | 179 | 0.687 | **0.045** | 0.91 |
| 2× | 358 | 0.686 | 0.027 | 0.91 |
| 3× | 537 | 0.674 | 0.017 | 0.99 |
| 6× (= strongest) | 1073 | 0.913 | 0.119 | 1.00 |
| 12× | 2135 | 0.895 | 0.136 | 1.00 |

Five distinct rungs. The fundamental and low rungs share a consistent ~0.68
identity with **tight IQR** (≤ 0.05) and **high occupancy** (≥ 0.9) — the
signature of a real-but-diverged monomer, not a noise period. The strongest is a
clean integer multiple (1073/179 = 5.99 ≈ ×6).

## Where it fits

A new per-array rescue inside `build_monomer_size_csv`, applied **only** to arrays
the existing passes left at `founder == strongest` (`multiplicity == 1`), i.e.
where founder recovery failed to go below the strongest because the basic monomer
was under the id gate. It does **not** touch:
- `repeat_type=SSR` arrays (founder = motif length),
- `fallback` arrays (Pass 5 handles those),
- arrays that already found a deeper basic (harmonic-basis Pass 1 / Lever 3/4 /
  kh_deeper) — those have `multiplicity > 1` already.

So it runs after the existing per-array passes, as a final "the strongest is
actually a HOR unit on a diverged base" check.

## Algorithm

```
for each array with founder == strongest (mult == 1), not fallback, not SSR:
    S = strongest_period
    # candidate fundamentals: peaks shorter than S of which S is ~an integer
    # multiple, with real-but-diverged support (id below the strict gate is OK).
    for P0 in peaks sorted ascending by period, P0 < S:
        k = round(S / P0)
        if k < _LADDER_KMIN:                continue          # excludes ×2 (k<3)
        if abs(S/P0 - k) > _RATIO_TOL:      continue          # S must be a clean rung
        if id_med(P0)  <  _LADDER_ID_FLOOR: continue          # not noise (>=0.60)
        if id_iqr(P0)  >  _LADDER_IQR_MAX:  continue          # consistent (<=0.10)
        if scan_occ(P0) < _LADDER_OCC_MIN:  continue          # tiles the array (>=0.50)
        # count DISTINCT integer multiples m in [1..k] that have a supporting peak
        rungs = { m for m in 1..k if any peak ~ m*P0 (within _RATIO_TOL)
                                     with id_med >= _LADDER_ID_FLOOR }
        if len(rungs) >= _LADDER_MIN_RUNGS and 1 in rungs and k in rungs:
            adopt founder = P0, multiplicity = k, method = "ladder";  break
```

Pick the **smallest** qualifying P0 (the true fundamental). Requiring `S` to be a
clean multiple of `P0` disambiguates against half-period artifacts (e.g. 90:
1073/90 = 11.9, not integer → rejected). For TRC_47 → P0 = 179, k = 6, rungs
{1,2,3,6,12} ≥ 3 → founder = 179, ×6.

## Constants (initial, tunable — grounded in TRC_47)

| const | value | rationale |
|---|---|---|
| `_LADDER_MIN_RUNGS` | 3 | ≥3 distinct integer multiples present (incl. P0 and S). **Excludes ×2** (only 2 rungs) — a 2-rung "ladder" is too easily coincidental. (Accepted: ×2 cases are not rescued.) |
| `_LADDER_KMIN` | 3 | the strongest must be ≥3× the fundamental (consistent with ≥3 rungs). |
| `_LADDER_ID_FLOOR` | 0.60 | P0 / rung id_med floor — a diverged real monomer, not noise. (TRC_47: 0.67–0.69.) Below the strict `_FOUNDER_ID_MIN = 0.70`; that gap is the whole point. |
| `_LADDER_IQR_MAX` | 0.10 | tight, consistent identity distinguishes a real monomer from a random period. (TRC_47: 0.017–0.045.) |
| `_LADDER_OCC_MIN` | 0.50 | the fundamental tiles a real fraction of the array. (TRC_47: 0.91.) |
| `_RATIO_TOL` (existing, 0.05) | — | integer tolerance for "S = k·P0" and rung matching. |

These are deliberately conservative; the before/after diff (below) is the safety
net for tuning.

## Output

- `founder_period = P0`, `strongest_period = S` (unchanged), `multiplicity = k`,
  `multiplicity_raw = S/P0`, `irregular = |S/P0 − k| > _RATIO_TOL` (TRC_47: false).
- New `founder_method = "ladder"`.
- **HOR-order tier**: map `ladder` → `supported` in `_hor_order_confidence`
  (the ×k order is well-supported by the clean integer ladder; it relies on a
  relaxed-id fundamental, so not `strict`). Open question 2.
- `founder_id_med = id_med(P0)` (so the report honestly shows the diverged
  ~0.69 identity of the basic monomer).

## Validation (incl. the before/after exploration you asked for)

1. **Before/after diff** with the committed helper:
   `tools/founder_diff.py test_data/Solanum_lycopersicum test_data/Arabidopsis_thaliana/run_e2e …`
   — lists every array whose founder columns changed. Expected: TRC_47 →
   founder 179 ×6; **inspect every other changed array** to classify it as a
   correct ladder rescue vs a false positive, and tune the gates accordingly.
2. **No-regression:** confirmed satellites must not change — TRC_2 (9491), TRC_4
   (53, already mult>1 so untouched), drapa/potato/arabidopsis satellites. The
   diff makes this explicit.
3. **Fixture + unit test:** extend `tests/data/kite_fixture/` with TRC_47's
   kitehor rows (+ clustering line) and add an assertion to
   `tests/test_ssr_raw_fixture.py` (or a new `test_harmonic_ladder.py`): TRC_47 →
   `founder_period == 179, multiplicity == 6, founder_method == "ladder"`. Also a
   negative fixture: an array with a strong strongest but **no** ≥3-rung ladder
   must keep `founder == strongest`.
4. `tests/unit.sh` green.

## Risks
- **False positives** (a coincidental ladder on a noise array). Mitigated by the
  four gates (≥3 rungs at clean multiples + id floor + tight IQR + occupancy) and
  caught by the before/after diff. The IQR + occupancy gates are the strongest
  discriminators — noise periods are wide-IQR / low-occupancy.
- **Does not help ×2 cases** (only 2 rungs) — accepted limitation per your call.
- Additive to the pipeline; report/TSV schema unchanged except the new
  `founder_method` value. `tools/founder_diff.py` proves the blast radius.

## Open questions for review
1. **Gate values** — start with the table above (tuned to TRC_47, conservative),
   then adjust from the before/after diff? Or different starting points?
2. **Tier mapping** for `founder_method = "ladder"`: `supported` (recommended —
   the order is ladder-proven but the base is relaxed-id) or `strict`?
3. ~~`_LADDER_MIN_RUNGS = 3` firm?~~ **Resolved: allow an *exceptionally-clean*
   ×2 (defined below).**

## Exceptionally-clean ×2 (Petr, 2026-06-11)

A ×2 ladder (P0 + 2·P0 = strongest, only 2 rungs) is allowed **only** when all of:

| gate | const | value | why |
|---|---|---|---|
| P0 real, not noise | `_LADDER_ID_FLOOR` | ≥ 0.60 | a half-harmonic of a real 2·P0 monomer compares unrelated monomer-halves → id falls below this, self-excluding |
| exceptionally consistent | `_LADDER_X2_IQR_MAX` | ≤ 0.05 | stricter than the 0.10 ≥3-rung gate |
| exceptionally complete tiling | `_LADDER_X2_OCC_MIN` | ≥ 0.85 | stricter than 0.50; P0 tiles ~the whole array |
| **double is a genuine HOR unit** | `_LADDER_X2_HOR_DELTA` | `id_med(2·P0) − id_med(P0) ≥ 0.10` | **decisive** — the 2× unit must be *meaningfully more conserved* than P0 (the HOR-conservation signature that P0 is the diverged base, not 2·P0 the real monomer with P0 a sidelobe) |

Sanity: a hypothetical TRC_47 179+358-only ladder **fails** the HOR-delta gate
(id 358 ≈ id 179), correctly — there 358 is just the 2× harmonic and the real
conserved unit is 1073 (×6). The ×2 path fires only for a true 2-monomer HOR
(monomer id ~0.65, conserved dimer id ~0.85).

Algorithm: `_LADDER_KMIN = 2`; for k≥3 use the ≥3-rung gate, for k==2 use the
four exceptional-×2 gates above.
