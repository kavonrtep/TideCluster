# Design note — HOR-order confidence tier (founder review Point 1)

Date: 2026-06-10
Status: IMPLEMENTED (branch fix/problematic-calls) — tier values
`none|strict|supported|weak`; HOR = strict ∪ supported; spectral classifier
removed; report schema v3. Decisions locked per Petr 2026-06-10: remove spectral
HOR everywhere, "weak" label, supported kept distinct but counted as HOR, high-k
constant = round(0.5/_RATIO_TOL_FRAC) = 50.
Branch context: `fix/problematic-calls` (off `main` @ 1.14.1)
Source: `docs/kitehor_founder_review_2026-06-10.md` Point 1 + Petr's high-k
TRC_4 instinct ("when HOR is really long compared to founder … we are already
reporting `~`").

## Problem in one sentence

The report treats **every** `multiplicity > 1` as an HOR array, which conflates
two claims of very different strength:

1. **founder recovered** — "the best basic monomer is 53 bp" (strong, this is
   what the founder pipeline is good at);
2. **confident HOR order** — "this array has a supported ×N higher-order
   structure" (weak when N is large, irregular, or rescue-derived).

## Why it matters (measured on the live `run_e2e` table)

`test_data/Solanum_lycopersicum/run_e2e/tc_kite/monomer_size_top3_estimats.csv`,
748 arrays:

| group | count | note |
|---|---|---|
| `multiplicity = 1` | 388 | no HOR claim; founder == strongest |
| `multiplicity > 1` | 360 | **all currently painted "HOR"** |
| ↳ strict, clean, low-k | 180 | confident HOR order |
| ↳ `irregular_multiplicity=true` | 170 | order approximate (already `~`-marked) |
| ↳ clean but k ≥ 50 | 5 | integer test vacuous at this k |
| ↳ pass3 / cluster / fallback | 4 | weakly supported rescue |

`mult>1` k: min 2, **median 15.5, max 368**. So ~**179 of 360** "HOR" arrays
(≈ half) are really *founder-recovered, order-approximate*. The 130 TRC_4 calls
the review flags (`x280` etc.) live here — they are correct *founder* calls
(53 bp basic) but their literal ×k is not a confidently supported HOR order.

**Key point: the discriminator already exists.** `irregular_multiplicity` plus
`founder_method` plus a high-k threshold separate the two groups with no new
computation on the founder side. This is a labelling/derivation change, not a
founder-selection change.

## Non-goals (explicit, per Petr's "no complicated rules")

- **No founder values change.** Founder periods, multiplicities, and
  `multiplicity_raw` stay exactly as they are. We only *annotate* them.
- **No new gates that suppress calls.** Nothing is dropped or down-weighted in
  the pipeline; this is report/diagnostic semantics (same philosophy as the
  0.13.2 "keep flag/badge, drop suppression" decision).
- **No re-tuning of Pass 1–5.** The tier reads their *outcome* columns.

## Proposed change

### 1. One derived column in `build_monomer_size_csv` — `hor_order_confidence`

A pure function of columns the row already carries
(`multiplicity`, `multiplicity_raw`, `irregular_multiplicity`,
`founder_method`, `founder_fallback`, `founder_id_med`). Four values:

| value | meaning | rule (mult>1 unless noted) |
|---|---|---|
| `none` | no higher-order structure | `multiplicity == 1` (founder == strongest) |
| `strict` | confident HOR order | `founder_method == strict` **and** not irregular **and** `k < _HOR_ORDER_HIGH_K` |
| `supported` | order plausible, family-backed | `founder_method in {pass2, kh_deeper}` **and** not irregular **and** `k < _HOR_ORDER_HIGH_K` |
| `approximate` | founder solid, **order not confidently integer** | everything else with `mult>1`: `irregular_multiplicity=true`, **or** `k ≥ _HOR_ORDER_HIGH_K`, **or** `founder_method in {pass3, cluster, fallback}` |

`_HOR_ORDER_HIGH_K` is **not a new free parameter** — it is the same vacuity
threshold the k-tolerance work already derived: above `k ≈ 0.5 / _RATIO_TOL_FRAC`
(= 0.5 / 0.01 = **50**) every k is within 0.5 of an integer, so the integer test
carries no information and the ×k is a *scaling estimate*, not a verified order.
Define `_HOR_ORDER_HIGH_K = round(0.5 / _RATIO_TOL_FRAC)` so it tracks the
constant rather than hard-coding 50.

SSR rows force `multiplicity = 1`, so they land in `none` automatically — no
special case.

Rationale for the buckets:
- `strict` is the existing strict-divisor success that is *also* clean and
  low-k — the only group where "×N HOR" is a verified statement.
- `supported` = the well-gated consensus paths (Lever-4 prevalent-founder anchor
  `pass2`, and `kh_deeper`). These have family/anchor support and passed the
  k-scaled tolerance, so the order is *plausible* but rests on cross-array
  evidence rather than this array's own clean divisor.
- `approximate` is deliberately the catch-all. It is the honest label for "we
  trust the founder; the ×k is an estimate." It already overlaps almost exactly
  with `irregular_multiplicity=true` (170 of the 179), so the report's existing
  `~` prefix is the natural visual partner.

### 2. Report relabel (`tc_rerender_report.py`) — the only user-visible change

Today the v0.12+ path treats `multiplicity > 1` as HOR for array
colouring/counts. Change the *framing*, not the data:

- **"HOR" count/badge** = arrays with `hor_order_confidence in {strict,
  supported}`. These are the confidently-ordered arrays.
- **`approximate`** arrays render as **"founder recovered (HOR order ≈)"** —
  keep the `~k` they already show, but do not count them as confident HOR.
- Per-TRC ideogram is already neutral (`a3f438b`), so no colour change needed
  there; this is about the **counts, badges, and the per-array detail table**.
- Add the tier to the per-array detail table and to the column legend (ties into
  the documentation work already requested).

### 3. Retire the spectral `compute_hor` HOR classification from reporting

**Decided (Petr, 2026-06-10): the spectral HOR notion is unreliable and should
be removed from reporting.** Spectral analysis (the top-3 monomer score
spectrum) is an **input to more detailed analysis, not a final HOR call.**

`tc_rerender_report.py` currently has a second, independent HOR notion derived
from the spectral top-3 monomer scores (m1/m2/m3):

- `compute_hor()` + `HOR_BIN_WEAK/MODERATE/STRONG` → "HOR strong/moderate/
  weak/none" (`tc_rerender_report.py:97-120`);
- `hor_badge()` renders those four states (`:341-348`);
- the genome/kite aggregates count `hor_strong/hor_moderate/hor_weak/no_hor`
  (`:1320-1325`);
- arrays are sorted/labelled by `hor_status` / `hor_confidence`
  (`:2338, :2386, :2433`).

This whole surface is the unreliable spectral classifier. The plan:

- **Remove the spectral HOR *classification* from the report** — drop
  `compute_hor`'s strong/moderate/weak/none status, `hor_badge`, the 4-bin
  aggregate counts, and `hor_status`/`hor_confidence`-driven sorting/labels.
- **Keep the raw spectral scores** (m1..m5 monomer sizes + scores) as table
  columns and as input to downstream/detailed analysis (dimer library, TAREAN,
  manual inspection). Only the *binned HOR verdict* leaves the report.
- **`hor_order_confidence` becomes the single HOR notion in the report**, on
  both modern and legacy paths. (Modern v0.12+ rows already carry the founder
  columns it needs; the legacy `compute_hor` recompute path at `:1080` is
  removed along with the rest.)

This is a **report schema bump** (`__version__` is currently `"2"`, described as
"4-bin HOR"; removing the 4-bin spectral HOR → `"3"`). `report.json` consumers
that read `hor_status`/`hor_confidence` must move to `hor_order_confidence`.

**CSV caveat — keep the vestigial `hor_status`/`hor_confidence` columns.** The
`monomer_size_top3_estimats.csv` already emits these two columns **empty** today,
and `tarean/consensus_prototype/consensus_ensemble.R`'s HOR_* alias shim selects
`HOR_status`/`HOR_confidence` after aliasing — it would error if the source
columns vanished. So on the **producer (`tc_utils.build_monomer_size_csv`)** side
we **keep** the empty `hor_status`/`hor_confidence` columns (zero change for the R
shim) and **add** `hor_order_confidence`. The spectral retirement is entirely a
**report-renderer** change; the dead `_kitehor_status`/`_KITE_HOR_BIN_*` helper
(defined but never called) can be removed as cleanup.

Note: this is a strictly *larger* simplification than the original note assumed
— there is now **one** HOR axis (founder arithmetic), not two. The reconciliation
problem disappears.

## Files touched (scope estimate)

- `tc_utils.py` — add `hor_order_confidence` derivation in
  `build_monomer_size_csv` (one helper + one column in the output dict at
  ~`tc_utils.py:3420`); add `_HOR_ORDER_HIGH_K` constant next to
  `_RATIO_TOL_FRAC`. ~25 lines, no change to founder logic.
- `tc_rerender_report.py` — read the new column; restrict the HOR
  count/badge to `{strict, supported}`; label `approximate`; add legend entry.
  **Also retire the spectral HOR classifier** (Section 3): remove
  `compute_hor`/`hor_badge`/the 4-bin aggregates/`hor_status`-driven sort, keep
  the raw m1..m5 scores; bump report `__version__` 2 → 3.
- `tests/test_strongest_by_identity.py` — assert tier values on the existing
  fixtures (strict low-k → `strict`; the case-(h) irregular → `approximate`;
  a high-k TRC_4-like row → `approximate`; mult=1 → `none`).
- Docs: column legend / detailed-table documentation (already on the to-do
  list) gains the tier description.

## Validation plan

1. **Unit:** extend the founder fixtures with explicit `hor_order_confidence`
   assertions across all four values; keep `tests/test_cluster_rescue.py` green.
2. **Solanum live table:** regenerate from `run_e2e` and confirm the bucket
   counts above (≈ 180 strict / ≈ 179 approximate / rest none); spot-check that
   TRC_4 high-k arrays are `approximate` and a clean TRC_2 9491 / strict
   satellite is `strict`.
3. **Cross-genome no-regression:** drapa / S. tuberosum / A. thaliana — confirm
   **founder periods and multiplicities are byte-identical** (only the new
   column appears); the change must be additive.
4. **Report:** rerender `run_e2e`; confirm the HOR count drops to the
   confident subset, `approximate` arrays read as founder-recovered, the spectral
   HOR badges/aggregates are gone, the raw m1..m5 score columns remain, and the
   legend documents the tier.
5. **Schema:** confirm `report.json` v3 has no `hor_status`/`hor_confidence`
   spectral fields and carries `hor_order_confidence` instead; check no
   downstream consumer still reads the removed fields.

## Open questions for Petr

1. **Label wording.** `approximate` / "founder recovered (HOR order ≈)" — or do
   you prefer `founder_only`, `weak`, etc.? ("weak" is now free since the
   spectral "HOR weak" is being retired.)
2. **`supported` vs folding into `strict`.** Should the well-gated `pass2`
   (Lever-4 anchor) / `kh_deeper` arrays count as confident HOR (merge into
   `strict`), or stay a distinct middle tier? They passed the k-scaled tolerance
   but lean on cross-array support. Default proposal: keep them distinct but
   *count them as HOR* in the report (i.e. HOR = strict ∪ supported).
3. ~~`compute_hor` reconciliation~~ — **resolved: retire the spectral HOR
   classifier from reporting (Section 3); keep raw spectral scores as input.**
4. **High-k constant.** Tie `_HOR_ORDER_HIGH_K` to `0.5/_RATIO_TOL_FRAC` (= 50,
   self-adjusting) as proposed, or pin an explicit value?
5. **Scope of removal.** Retire the spectral HOR classifier on the **modern
   path only** and leave it for legacy CSVs, or remove it **everywhere**
   (cleaner; legacy reports lose the spectral HOR badge)? Default proposal:
   remove everywhere — one HOR notion, less code.
