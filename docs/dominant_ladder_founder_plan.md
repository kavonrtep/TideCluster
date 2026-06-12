# Dominant-score harmonic-ladder founder (Pass 7b)

Status: **implemented** (`tc_utils.build_monomer_size_csv`, branch
`feat/ladder-anchored-founder`).

## Problem

For a **divergent satellite**, the basic monomer has only moderate per-copy
identity (FabTR-53: ~188 bp at `id_med ≈ 0.90`) while long, multi-monomer
periods are *more* self-similar (`id_med ≈ 0.98`). TideCluster's Pass 1 anchors
`strongest = argmax(identity_med)`, then recovers the monomer by a divisor search
requiring `strongest = k·P0` with `|k − round(k)| ≤ _RATIO_TOL`. That works only
when the highest-identity peak is a **clean integer multiple** of the monomer.

When `argmax(id_med)` lands on a long, **off-ladder** period (FabTR-53
`ptg000030l_rc:2-44995`: strongest = 3150, monomer = 188, `3150/188 = 16.76`),
the decomposition fails and the founder either:

- **collapses** to the long period (`founder = strongest`, `mult = 1`, method
  `none`) — reported as a clean 3150 bp monomer with no HOR, or
- lands on a **non-monomer** `strongest/round(k)` value (165, 168, 347, 1808)
  with low identity.

Pass 7 (`_harmonic_ladder_founder`) does not help: it also anchors the ladder on
`strongest` (`strongest = k·P0`), so an off-ladder strongest disarms it too.

Evidence: 6 of 43 FabTR-53 arrays. In 5 of the 6, the **rank-1-by-score** peak is
the monomer (187–189) **and** kitehor's own `founder_period` = the monomer — the
right answer is present in the data but discarded.

## Fix — Pass 7b

A pass complementary to Pass 7. Anchor the ladder on the **rank-1-by-score** peak
P0 (the true monomer of these satellites), independent of `strongest`. Adopt P0
as founder when **all** hold:

1. **Incommensurate guard** — the current founder is *not* commensurate with P0
   (neither `founder = k·P0` nor `P0 = k·founder` within `_RATIO_TOL`, absolute).
   Leaves every array already consistent with its dominant monomer untouched:
   the 37 correct FabTR-53 arrays, the canonical 178×10 case (founder = P0 = 178,
   k = 1), Pass 7's own 179×6 / 111×3 (commensurate), and ordinary HOR calls.
2. **Self-standing ladder** — P0 carries `≥ _DOM_LADDER_MIN_RUNGS` (3) distinct
   integer-multiple rungs `m·P0` present among the peaks (incl. m = 1), with
   `id_med(P0) ≥ _LADDER_ID_FLOOR`, `identity_iqr ≤ _LADDER_IQR_MAX`,
   `scan_occupancy_frac ≥ _LADDER_OCC_MIN` (the k ≥ 3 gates of Pass 7).
3. **kitehor corroboration** — kitehor's own `founder_period` agrees that P0 is
   the monomer (`|kitehor_founder − P0| ≤ _DOM_LADDER_KH_TOL · max`). An
   independent (k-mer-autocorrelation) method must concur, so a lone high-score
   peak that kitehor does *not* call a founder (e.g. an array with a genuine
   ~12 kb HOR unit — slyco TRC_13) is not overridden.

On adoption: `founder = P0`, `strongest` unchanged, `multiplicity = round(S/P0)`
(`irregular` when non-integer — the honest state for an off-ladder strongest),
`founder_method = "dominant_ladder"`. The HOR-order tier is **weak** (founder
recovered, but the ×k order is not confidently integer) — `dominant_ladder` is
deliberately not in the `supported` set of `_hor_order_confidence`.

## Validation

`tools/founder_diff.py` (working tree vs HEAD) across slyco + arabidopsis (run_e2e,
run_e2e_long) + drapa (run_0131, run_0132):

- **5 FabTR-53** target arrays flip to 187–189 (3150→188, 347→188, 1808→187,
  168→189, 165→188).
- **17 further arrays** flip on the other genomes — **every one corroborated by
  kitehor** (`new = kitehor_founder_period`) with equal-or-higher identity:
  arabidopsis TRC_1 160→178, slyco TRC_4 48/50→53, drapa TRC_103 4864→69, etc.
- The kitehor-corroboration gate excludes the one outlier (slyco TRC_13
  12252→53, where kitehor's founder = 12254 and identity drops 0.99→0.64).
- `tests/unit.sh` stays green — 179×6, 111×3, 178×10 and drapa TRC_158 unchanged.

Regression guard: `tests/test_dominant_ladder.py` (synthetic peaks; no external
data). See also `docs/harmonic_ladder_founder_plan.md` (Pass 7) and the founder
semantics note.
