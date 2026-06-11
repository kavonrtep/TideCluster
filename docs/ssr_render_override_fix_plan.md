# Fix plan — retire the render-time SSR founder override + show *raw* SSR coverage

Date: 2026-06-11
Status: IMPLEMENTED. Recommended options taken (per-TRA raw only; consensus
columns kept + renamed `ssr_consensus_*`; shipped together). Verified: TRC_2
chr2:15535707 → Founder 9520, SSR "ATC 6%"; TRC_18 per-array raw varies
(68/82/16…%) with consensus callout "ATC 98.6%"; founders byte-identical on
slyco + arabidopsis (SSR columns only changed); regression guards + suite green.
Root cause: `docs/` analysis of TRC_2 `chr2:15535707-15614778` (this session).
Scope: `tc_rerender_report.py`, `tc_utils.py` (+ regenerate the bundle/CSV and
re-render). v2 report only.

## Problem recap

A clean ~9.5 kb satellite (`founder_period=9520`, `id_med 0.99`, clustering
`repeat_type=TR`) is shown in the report as **`founder=3`, SSR (ATC 96.2 %)**.
Two independent faults:

1. **Render-time override.** `tc_rerender_report._apply_ssr_founder_override()`
   (called in `load_kite_top3`, line ~1109) re-applies the *removed* pipeline
   override: when `ssr_total_coverage_pct ≥ 95`, it overwrites
   `founder = strongest = m1` (top kite peak = 3 bp ATC). It ignores the
   clustering `repeat_type=TR` and the CSV's `ssr_founder_override=false`, so it
   clobbers the correct `9520` at render time. (Pipeline fix `92d3db8` removed
   the equivalent override from `build_monomer_size_csv`; the report kept a
   private mirror `_SSR_OVERRIDE_COVERAGE_PCT = 95.0`.)
2. **Misleading coverage.** The `96.23 %` is a `consensus_single` artifact
   (kitehor measures ATC against an ATC-periodic *consensus*). The **raw**
   per-array coverage is `ssr_raw_dominant_motif_coverage_pct = 6.39 %`.
   TideCluster's monomer CSV carries only the consensus figure, so the ≥95 % gate
   trips on an artifact and the SSR column reads "96 %" for a 6 %-ATC array.

## Conceptual model — three SSR figures at two levels (Petr, 2026-06-11)

SSR information lives at two distinct levels and must be reported as such:

| # | figure | where it comes from | level | report it? |
|---|---|---|---|---|
| **1. TRC consensus call** | TideHunter **consensus** SSR detection **precedes clustering**; clear SSRs are clustered into a TRC by prevalent motif → clustering GFF3 `repeat_type=SSR; ssr=MOTIF (pct%)` (the `%` is on the consensus). Carried as `trc.ssr_motif`. | **per-TRC** | **Yes — for SSR-typed TRCs**, as the family-defining call (e.g. "SSR family ATC, 98.6 % consensus"). |
| **2. per-array consensus** | kitehor `ssr-scan` `consensus_single` `dominant_motif_coverage_pct` — measured against an artificial motif-periodic consensus. Inflated (96 % on a 6 %-ATC satellite; 100 % on real SSR arrays). | per-array | **No — drop from display** (the misleading artifact, and the thing the removed override keyed on). |
| **3. per-array raw** | kitehor `ssr-scan` `ssr_raw_*` — the actual per-sequence composition. Varies array to array (TRC_18: 68–82 %; the TR satellite: 6 %), reflecting real variability. | **per-TRA** | **Yes — for every TRA**, as the honest per-array SSR content. |

Design consequence:
- **Every TRA** reports **#3 (raw)** — real content, often much lower than the
  consensus, reflecting array variability.
- **SSR-typed TRCs** *additionally* report **#1 (consensus call)** — the
  authoritative "this TRC is an SSR family of motif X" statement that defined the
  cluster. (Already surfaced today via the "Simple Sequence Repeat" callout +
  `trc_table.ssr_motif`; the fix only needs to **label it as consensus** so it
  reads distinctly from the per-array raw figure.)
- **#2 is dropped** from display (kept in the raw TSV only, if at all).

---

## Part 1 — retire the render-time SSR founder override (primary)

The pipeline already makes the authoritative SSR decision (clustering
`repeat_type` → `build_monomer_size_csv`), writing the final `founder_period`,
`ssr_founder_override`, etc. into the CSV. The report must **trust those values**,
not re-derive them.

- **Remove** `_apply_ssr_founder_override()`, its call in `load_kite_top3`
  (the `if is_v012: _apply_ssr_founder_override(...)` block), and the
  `_SSR_OVERRIDE_COVERAGE_PCT` constant (`tc_rerender_report.py:~1105-1140`).
- **Also remove** the now-dead `tc_utils._SSR_OVERRIDE_COVERAGE = 95.0`
  (`tc_utils.py:2387` — defined, no longer referenced).
- The report then renders `founder_period`/`ssr_founder_override` exactly as the
  CSV holds them.

**Legacy note.** The override's stated purpose was "render pre-fix (drapa-style)
CSVs correctly without re-running". After the `repeat_type` fix, current CSVs
already carry correct SSR founders, so the override is now net-harmful (false
positives like this one). Pre-1.13 CSVs that genuinely need an SSR founder should
be **re-run through `build_monomer_size_csv`** rather than patched at render time.
(Acceptable: those are stale archival dirs; the curated path is a re-render from
fresh kite output.)

**Effect:** TRC_2 `chr2:15535707` renders `founder=9520`, no SSR badge.

---

## Part 2 — carry & display *raw* SSR coverage (secondary)

`kitehor.ssr.tsv` already has both views per array:

| view | columns (kitehor.ssr.tsv) | this array |
|---|---|---|
| consensus (current) | `dominant_motif`, `dominant_motif_coverage_pct`, `total_ssr_coverage_pct`, `top_motifs` | ATC, **96.23 %** |
| **raw** (wanted) | `ssr_raw_dominant_motif`, `ssr_raw_dominant_motif_coverage_pct`, `ssr_raw_total_coverage_pct`, `ssr_raw_n_regions`, `ssr_raw_top_motifs` | ATC, **6.39 %** |

`build_monomer_size_csv` currently lifts only the consensus columns
(`tc_utils.py:3517-3522`).

### 2a — carry raw fields into the monomer CSV (`tc_utils.build_monomer_size_csv`)
Add columns from the `ssr_raw_*` group:
`ssr_raw_dominant_motif`, `ssr_raw_dominant_motif_coverage_pct`,
`ssr_raw_total_coverage_pct`, `ssr_raw_n_regions`. Keep the consensus columns
(transparency + `consensus_period_bp` is still used) — purely additive, so the
CSV gains 4 columns and existing consumers are unaffected. Missing on
pre-0.13.2 kite output ⇒ blank via `.get(default "")`.

### 2b — per-TRA display = raw (#3), for every array
Switch every per-array SSR *display* site from the consensus figure (#2) to the
raw one (#3):
- `_ssr_cell()` (`:2678`) — the per-array table **SSR column**: show
  `ssr_raw_dominant_motif` + `ssr_raw_dominant_motif_coverage_pct` (e.g.
  "ATC 6 %", or for a TRC_18 array "ATC 82 %"); "—" when raw is absent/zero.
  Fall back to the consensus value only on legacy CSVs that lack `ssr_raw_*`.
- `structure_cell()` (`:538` uses `ssr_total_coverage_pct`) — use the raw total.
- The Details child-row SSR line: `_array_detail_payload()` builds
  `ssr = [motif, dom_pct, tot_pct, top]` from the consensus fields → switch to
  raw. (JS unchanged; it just renders the payload numbers.)
- `load_kite_top3()` (`:1045-1049`) — read the new raw fields into the model.
- **Drop #2** (per-array consensus) from display entirely; keep it in the CSV/TSV
  for transparency only (or drop — open question 2).

### 2c — per-TRC display = consensus call (#1), labelled, for SSR TRCs
The SSR-TRC callout already prints `trc.ssr_motif` (`tc_rerender_report.py:3415`).
Make it read explicitly as the **consensus / family** call so it can't be
confused with the per-array raw figure, e.g.:

> **Simple Sequence Repeat.** Consensus motif (TideHunter consensus, used to
> cluster this family): `ATC` — **98.6 %** on the consensus. Per-array SSR
> content (raw, in the table below) is typically lower and varies between arrays.

`trc.ssr_motif` already carries `ATC (98.6%)`; no new data needed — wording only.

### 2d — update the data dictionary (`COLUMN_DICT["tra_table"]`)
Add ColSpecs for the new raw columns; descriptions name the level explicitly
("raw per-array SSR content" vs the TRC consensus call). `columns.tsv` + the
legend update automatically (R1 single-source).

### Note — genuine SSR TRCs read correctly
For a real microsatellite array (`repeat_type=SSR`, TRC_18 ATC) the per-array raw
coverage is **68–82 %** (high, but variable — the honest signal), while the TRC
callout still states the **98.6 % consensus** family call. The motif-rich TR
satellite drops from a misleading 96 % to an honest ~6 % and keeps `founder=9520`.

---

## Validation

1. **Re-render** the bare `test_data/Solanum_lycopersicum` report and confirm
   TRC_2 `chr2:15535707-15614778` now shows **Founder 9520, no SSR badge**, and
   its per-array SSR column reads the **raw** ~6 % (or "—").
2. **Genuine SSR, two levels:** TRC_18 / TRC_5 still show the SSR founder
   (= motif length); the **TRC callout** states the **consensus** family call
   (ATC 98.6 %); the **per-array** SSR column shows **raw** (e.g. 68–82 %, varying
   between arrays) — not a flat 100 %.
3. **No founder regressions:** regenerate the monomer CSV on slyco + arabidopsis
   and assert `founder_period`/`multiplicity`/`hor_order_confidence` columns
   **byte-identical** (Part 2 is additive; Part 1 is render-only); only the 4 new
   `ssr_raw_*` columns appear.
4. **Report bundle:** TSV headers gain exactly the 4 raw columns; `columns.tsv`
   documents them; unit suite (`tests/test_report_exports.py` incl. the
   legend ≡ COLUMN_DICT guard) green.
5. Add a small unit fixture: an array with consensus-cov ≥95 but raw-cov low and
   `repeat_type=TR` must render `founder = founder_period` (not m1) — locks the
   regression so the render override can't creep back.

## Risks / rollback
- **Low.** Part 1 is a deletion (render trusts the CSV); Part 2 is additive
  columns + display source-swap. Pre-0.13.2 kite output lacks `ssr_raw_*` →
  graceful blank + consensus fallback in display.
- Old pre-1.13 CSVs that relied on the render override for an SSR founder now need
  a `build_monomer_size_csv` re-run (documented above).
- Rollback = revert the commits; no data-format break (additive columns).

## Open questions for review
1. **Per-TRA column shows raw only** (recommended — the per-array headline is the
   honest content; the consensus call lives at the TRC level), or also echo the
   consensus on hover?
2. The **per-array consensus (#2)** columns in the monomer CSV / TSV: keep for
   transparency (clearly named `..._consensus_...`) or drop entirely? (Recommend
   keep but rename so they can't be mistaken for raw.)
3. Land Part 1 immediately (quick win, re-render) and Part 2 as a follow-up, or
   ship both together so the re-render happens once? (Recommend together.)
