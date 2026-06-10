# Implementation plan — Solanum_lycopersicum problematic calls

Fixes for the six items in `problematic_calls.md` (root causes in
`problematic_calls_analysis.md`). Target: kitehor 0.13.2 / post-1.14.1.

## Findings that shaped the design (measured)

- **SSR is classified at clustering, not in kite.** `TideCluster.py:552-660`:
  TRA consensus with **SSR proportion > 0.9** is pulled out of clustering,
  grouped by motif, and written to the GFF3 as `repeat_type=SSR` + an `ssr`
  motif (e.g. TRC_5 `AAT (90.6%)`, TRC_18 `ATC (98.6%)`). The % is on the
  **consensus**, so it hides per-array noise (TRC_18 TRAs are 16–100 % ATC on
  the full sequence but the family is one ATC cluster). `repeat_type=TR` TRCs
  (TRC_2 = 9491 satellite, TRC_4 = 53 bp) are normal. **`build_monomer_size_csv`
  never reads `repeat_type`** — that blindness is the root of the SSR cases.
- **kitehor's ssr-scan is secondary** — "this 500 bp founder is also SSR-rich"
  annotation, NOT the classifier. Driving the founder off its per-array coverage
  (the 95 % override) is the bug.
- **min-period:** rescoring TRC_18 at `--min-period 2` makes P=3 scorable but
  ambiguous (chr8 id 1.0, chr3 id 0.667) — **not pursued**; the SSR classifier
  already names the motif.
- **max-period:** TRC_10 (112 kb) at cap 15000 → junk `35 (cov 0.03)` in 2.8 s;
  at 25000 → real `16326 (id 0.9995, cov 0.985)` in 5.9 s; cost is all in the
  long-period O(period²) work → confine it to the few above-cap arrays.

## Sequencing (by risk)

- **PR-A** reporting (low risk): neutral per-TRC ideogram + legend.
- **PR-B** case 4 (isolated): selective max-period extension.
- **PR-C** founder core: C1 SSR-via-repeat_type (cases 1/3/5) + C2 harmonic-basis
  (case 2). Re-validate satellites.

---

## PR-A — reporting

**A1 — HTML ideogram → positions only (HTML bug).** `tc_rerender_report.py`
`_class_fill`/`_render_ideogram`/`render_trc_distribution`: render the **per-TRC**
distribution ideogram with a single neutral fill (drop the HOR/subrepeat/SSR/
fallback colour cascade for that view; TRC_3 has no HOR yet shows structural
colours). *Open: keep colour on the genome-wide index overview, or strip
everywhere?*

**A2 — `*` vs SSR legend.** Surface the `*` (rescore-NA fallback) / SSR key on
the index + genome views, not just the collapsed per-TRA legend. (Mostly
dissolves once PR-C gives SSR TRCs a uniform SSR badge.)

---

## PR-B — selective max-period extension (case 4)

`TideCluster.py` kite step. Two-pass, selective:
1. First `rescore` at the normal cap (`kite_rescore_max_period`).
2. Flag arrays whose **rank-1-by-kite-score** peak has `period > cap` (→ NA) with
   no credible founder below the cap.
3. If any: build a small multifasta of just those + their peaks, re-`rescore` at
   a high cap (`--kite_rescore_max_period_ext`, default 25000 under `--long`, or
   `0` = unlimited), **merge** those `case_id` rows back into
   `rescored.peaks.tsv`.

New arg `--kite_rescore_max_period_ext`. Empty flag set → zero extra cost.
**Validation:** TRC_10 → 16326; confirm small wall-time delta on a full `--long`
run.

---

## PR-C — founder core

### C1 — SSR founder from `repeat_type` (cases 1 + 3 + 5)

**The whole SSR cluster collapses to one change: make the founder pipeline read
the clustering SSR classification instead of re-deriving it from kitehor coverage.**

- **Plumb it in.** Pass a per-TRC `{repeat_type, ssr_motif}` map (parsed from
  `tc_clustering.gff3`'s `repeat_type` + `ssr` attributes) into
  `build_monomer_size_csv` (new input; the kite step runs after clustering so the
  GFF3 exists).
- **`repeat_type=SSR` TRC → founder = SSR motif length** (fundamental: `ATC`→3,
  `AAT`→3) for **every** array, `method=ssr`, SSR badge — regardless of kite
  peaks / coverage / rescore. No HOR decomposition. → **fixes case 5 (TRC_18 all
  3) and case 3 (TRC_5 all 3, uniform badge, no `*`).** Motif length taken from
  the clustering `ssr` attribute (authoritative family call).
- **`repeat_type=TR` TRC → normal kite logic, and DELETE the kitehor-coverage
  Pass-4 override** (`_SSR_OVERRIDE_COVERAGE = 95`). TRC_2's lone 96 %-ATC array
  keeps the family founder 9490 (wins on id_med 0.9875 once the override is gone).
  kitehor's SSR info stays in the row as **annotation** ("9490 founder, ATC-rich
  96 %"). → **fixes case 1 (TRC_2 → 9491).**

**Validation:** TRC_18 all 3; TRC_5 all 3 (badged, no `*`); TRC_2 all ~9491
(method none, no ssr override); the other 3 SSR TRCs → motif length; TR
satellites elsewhere unchanged.

### C2 — Harmonic-basis founder selection (case 2 — independent of SSR)

TRC_4 basic is **53 bp** (a satellite monomer, not an SSR — `repeat_type=TR`). 53
is the top kite peak and the peaks are a harmonic series (53, 106, 159…), but
`strongest = argmax(id_med)` wanders onto a long noisy 7617 (id 0.768 vs 0.755)
and the divisor search (`_KMAX = 30`) can't reach 53 (k≈144) → settles on 318.

- **Change:** `tc_utils.build_monomer_size_csv` Pass 1 — let `P0` = rank-1-by-kite-
  **score** peak; if `P0` passes the id gate and the strong peaks are ~integer
  multiples of `P0`, adopt `founder = P0` (multiplicity = round(strongest/P0)),
  bypassing the `_KMAX` ceiling on this path. Use kite *score* to certify the
  basic.
- **Risk:** highest of the set (touches founder core). **Validation:** TRC_4 → 53
  for the ~93 stragglers AND the 276 already-53 unchanged; unit fixtures green;
  drapa/S. tuberosum/A. thaliana satellites (178/309/CEN178) unchanged.

---

## Cross-cutting validation (after PR-C)

1. `tests/unit.sh` green; add fixtures for repeat_type=SSR founder (C1) and
   harmonic-basis (C2).
2. Re-run the kite step on drapa + S. tuberosum + A. thaliana + S. lycopersicum:
   - **Unchanged:** drapa TRC_6/16/62, S. tuberosum TRC_1, A. thaliana CEN178.
   - **Fixed:** S. lycopersicum TRC_2 → 9491, TRC_4 → 53, TRC_10 → 16326, TRC_18
     → 3, TRC_5 → 3.
3. Confirm the selective max-period pass adds little wall time.

## Open decisions for Petr
- A1: neutral per-TRC ideogram only, or strip class colour on the genome overview
  too?
- C1: motif length from the clustering `ssr` attribute (recommended) vs kitehor
  `dominant_motif_length`; confirm "no HOR for SSR TRCs."
- B: confirm selective max-period (recommended) vs a `--long`-scaled global cap.
