# Solanum_lycopersicum problematic calls — root-cause analysis

Analysis of `problematic_calls.md` against `test_data/Solanum_lycopersicum`
(kitehor 0.13.2, 1.14.1 code; `--long`, `kite_rescore_max_period 15000`).
**No code changed** — observations, root causes, and proposed fixes only.

Two themes run through these: **(A) short-period handling** (cases 1, 3, 5 — the
SSR override and the `rescore --min-period 5` floor) and **(B) basic-monomer /
long-period founder selection** (cases 2, 4 — the `_KMAX` divisor cap and the
`max-period` cap). The user's own one-liner — *"unreported id_med for very short
founders"* — is the crux of theme A.

---

## Case 1 — SSR override masks a real long-period family founder  *(founder bug, high priority)*

**Observation.** `TRC_2 chr2:15535707-15614778` → `founder=3` (method `ssr`),
`ssr_dominant_motif=ATC`, `ssr_total_coverage_pct=96.23`. But the rescored peaks
show the real founder is there and strong:

| rank | period | score | id_med | cov | occ |
|--|--|--|--|--|--|
| 1 | **3** | 0.058 | **NA** | NA | NA |
| 4 | **9490** | 0.022 | **0.9875** | 0.93 | 0.96 |

The TRC_2 family founder is **9491** (12 of 24 arrays); this array is the lone
outlier at 3.

**Root cause.** Pass 4 SSR override (`tc_utils.py` ~3011, `_SSR_OVERRIDE_COVERAGE
= 95`) fires unconditionally when `ssr_total_coverage_pct ≥ 95 %` and forces
`founder = strongest = kite top peak` (= the 3 bp ATC). The 9490 bp satellite
monomer is itself ATC-rich, so the SSR scanner reports 96 % ATC coverage even
though the true tandem unit is 9490 bp. **Without Pass 4 this array is correct:**
Pass 1 picks `strongest = argmax(id_med) = 9490` (id 0.9875) → `founder = 9490`;
Pass 4 then clobbers it. The override conflates *SSR composition* with *the repeat
being an SSR*.

**Discriminator.** In a genuine `(ATC)ₙ`, every high-id_med period is a multiple
of 3. Here `9490 / 3 = 3163.3` (not a multiple) yet has id_med 0.9875 → there is
real non-SSR higher-order structure.

**Proposed fix.** Gate the override: **skip Pass 4 when the array already has a
confident non-SSR founder** — e.g. the pre-override `founder_id_med ≥ 0.85` with
good coverage AND the founder period is *not* an integer multiple of the SSR
motif length (or it matches the TRC consensus `consensus_by_trc[trc]`). Pure SSRs
(every strong period a motif multiple) still override; real ATC-rich satellites
like TRC_2 keep their 9490 founder. The SSR signal can still be surfaced as an
annotation rather than as the founder.

---

## Case 2 — basic monomer (53) lost when `strongest` is a long noisy period beyond the divisor cap  *(founder bug, high priority)*

**Observation.** `TRC_4` prevalent founder is **53** (276 of 369 arrays). For
`chr6:9719494` the top kite peak IS 53 and the peaks form a clean harmonic
series, yet `founder=318, strongest=7617`:

| rank | period | score | id_med |
|--|--|--|--|
| 1 | **53** | 0.139 | 0.7547 |
| 2 | 106 (2×53) | 0.052 | 0.7453 |
| 3 | 159 (3×53) | 0.026 | 0.7421 |
| … | 212,265,318 (n×53) | … | ~0.74 |

`strongest = argmax(id_med) = 7617` (id 0.768, *marginally* above 53's 0.755).
`7617` is not a clean multiple of 53 (143.7×); the divisor search caps at
`_KMAX = 30`, so it reaches 318 (`7617/24`) but never 53 (needs k≈144). ~93
TRC_4 arrays land on a multiple of 53 instead of 53.

**Root cause.** Two compounding issues in `build_monomer_size_csv` Pass 1:
1. `strongest = argmax(identity_med)` picks a long, noisy period on a tiny
   id_med margin (0.768 vs 0.755), even though the **kite score** overwhelmingly
   identifies 53 as the basic (0.139 vs ~0).
2. The divisor search (`_KMAX = 30`, `tc_utils.py:2106`) cannot bridge a long
   `strongest` down to a short basic, so it settles on an intermediate multiple.

**Proposed fix (recommended: harmonic-basis detection).** When the rank-1-by-kite-
score peak `P` is short, passes the id gate (`id_med ≥ 0.7`), and the other strong
peaks are (approximately) integer multiples of `P`, adopt `P` as the founder
directly — the kite-score-dominant base of a harmonic series. This uses the kite
score for what it is best at (periodicity/basic detection) and is robust to the
`strongest`-by-id_med pick wandering onto a long noisy period. (Alternatives:
extend the divisor search to test the top kite peak as a divisor of `strongest`
with large `k` under the id gate; or raise `_KMAX` — blunter, with side effects.)

---

## Case 3 — `*` (fallback) vs SSR badge is opaque for the same short repeat  *(reporting; founder value is correct)*

**Observation.** `TRC_5 chr4:32851204-32864884` → `founder=3` (method
`fallback`, so a red `*` in the report), `ssr_dominant_motif=AAT`,
`ssr_total_coverage_pct=12.08`. Its TRC_5 siblings (high SSR coverage) get
`founder=3` via method `ssr` (SSR badge). Same 3 bp AAT repeat, two different
markers.

**Root cause.** A real short (3 bp AAT) repeat whose SSR coverage is *below* the
95 % override threshold falls through to the fallback path (kitehor returns NA at
P=3 < `min-period 5`) → `*` marker, while its high-coverage siblings get the SSR
badge. The `*` (rescore-NA fallback) and the SSR badge are distinct mechanisms;
the per-TRA legend (`arrays_legend()`) explains `*` but it is collapsed/absent at
the genome/index view where the user was looking.

**Proposed fix.** (a) Surface the `*` / SSR legend where these markers actually
appear (index + genome-distribution views, not only the collapsed per-TRA
legend). (b) Better: when a fallback founder period equals a detected SSR motif
length with non-trivial coverage, show the **SSR badge** instead of `*` so the
same biological repeat reads consistently. (Tied to case 5 — fixing the short-
period scoring removes most of these fallbacks anyway.)

---

## Case 4 — real long monomer above `max-period` → spurious short fallback founder  *(founder bug + parameter, medium priority)*

**Observation.** `TRC_10 chr10:59071587-59183483` (112 kb array) → `founder=35`
(method `fallback`, id_med 0.57, **cov 0.03** — junk). The real monomer is right
there but uncomputable:

| rank | period | score | id_med |
|--|--|--|--|
| 1 | 35 | 0.0157 | 0.57 (cov 0.03) |
| 2 | **16326** | 0.0151 | **NA** (> max-period 15000) |

**Root cause.** `kite_rescore_max_period` (15000) caps rescore; periods above it
get NA, so a real ~16.3 kb monomer can't be scored and the founder falls back to
a near-zero-coverage 35 bp peak. With `--long`, TideHunter searches monomers up
to 25 000 nt but the rescore cap is **not** raised to match (`TideCluster.py:148`,
default 10000; the run used 15000).

**Proposed fix.**
1. **Scale the cap with `--long`**: default `kite_rescore_max_period` to ≥ 25000
   when `--long` is set (match TideHunter's `--long` monomer ceiling).
2. **Selective extension**: when the top-scored peak is above the cap (NA) and no
   credible founder exists below it, do one targeted rescore for that array with a
   raised cap — cheaper than globally raising it (rescore is O(period²)). This is
   exactly the user's "extend search selectively."
3. At minimum, never prefer a `cov ≈ 0` fallback peak over an above-cap peak that
   dominates the kite score — surface that a strong peak was above the cap.

---

## Case 5 — very short real founders (P = 3, 4) missed: `rescore --min-period 5` → NA  *(founder bug, high priority)*

**Observation.** `TRC_18` is an ATC 3 bp family, but founders are called as
multiples of 3:

| array | top kite peak | id_med(P=3) | called founder |
|--|--|--|--|
| chr3:21052539 | **3** (score 0.257) | **NA** | 9 |
| chr8:29790032 | **3** (score 0.329) | **NA** | 6 |

The harmonic peaks 6/9/12/15 all get id_med ~0.83–0.92; `strongest = argmax →
9/12`, founder lands on 6/9.

**Root cause.** TideCluster's kite step does **not** pass `--min-period` to
`kitehor rescore` (`TideCluster.py` ~167), so it uses kitehor's default **5** →
periods 3 and 4 get NA id_med → TC cannot select them and lands on the shortest
*scored* multiple. This is precisely the user's "unreported id_med for very short
founders."

**Proposed fix.**
1. **Pass `--min-period 2` (or 3)** to `kitehor rescore` so P=3,4 get real
   id_med. **This is now safe** because kitehor 0.13.2 gives realistic id_med to
   short periods — a real 3 bp ATC tandem scores high, phantoms score low (the
   0.13.2 fix we just shipped). Verify 0.13.2 handles P=2/3 cleanly (the slop
   kernel clamps to P−1 for tiny P).
2. **Combine with case-2 harmonic-basis detection**: even once P=3 has id_med,
   the founder logic must prefer the shortest basic (3) over its multiples (6,9).
   The two fixes are complementary — (1) makes the signal available, (2) selects
   the basic.

---

## HTML report bug (TRC_3) — ideogram colours a non-HOR TRC by structural class  *(reporting; user wants positions only)*

**Observation.** `TRC_3` has 26 arrays, **all multiplicity = 1** (no HOR), 17
with a subrepeat. The per-TRC distribution ideogram nonetheless shows structural
colours, reading as "HOR + plain TR" for a TRC that has no HOR.

**Root cause.** `_class_fill` (`tc_rerender_report.py`) applies per-array semantic
colouring on the per-TRC ideogram via the v0.12 fallback cascade
(`multiplicity>1 → HOR green`, `subrepeat → teal`, `ssr → purple`, `fallback →
red`). These per-array signals don't reflect the TRC-level classification, so a
non-HOR TRC shows HOR/subrepeat colours (here the 17 subrepeat arrays → teal,
plus any fallback → red).

**Proposed fix (user's stated preference).** Render the per-TRC distribution
ideogram with a **single neutral fill** — show only TRA/TRC *positions*, drop the
class colour-coding. Minimal change in `_render_ideogram` / `_class_fill`
(neutral colour for the per-TRC ideogram); optionally keep the coloured class-mix
only on the genome-wide overview, or remove it everywhere per preference.

---

## Suggested priority

1. **Case 5** (`--min-period`) — one-line wiring change, unblocks all P=3/4
   founders; safe under 0.13.2.
2. **Case 1** (SSR-override gate) — clear discriminator, fixes a confident
   miscall on a major family.
3. **Case 2** (harmonic-basis founder) — biggest correctness win (≈93 TRC_4
   arrays), but the most delicate change to the founder logic; pairs with case 5.
4. **Case 4** (max-period under `--long`) — simple default change + optional
   selective extension.
5. **HTML ideogram** + **Case 3 legend** — reporting-only, low risk.

Cases 1/2/5 should be validated against the existing `tests/unit.sh`
founder fixtures and re-checked on drapa + S. tuberosum + S. lycopersicum so the
short-period and basic-monomer changes don't regress the satellites we just
validated for 1.14.1.
