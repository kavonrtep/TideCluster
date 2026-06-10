# Implementation plan — report content single-source refactor (R1 + R2 + R3)

Date: 2026-06-10
Status: IMPLEMENTED on main (3 commits: R1 4ffd637, R2 ecfbf51, R3 3c450fc).
All open questions approved (canonicalised tier wording; BLURBS in Python;
legend≡COLUMN_DICT CI assertion added). TSV bundle headers verified
byte-identical on slyco + arabidopsis; unit suite green.
Source: `docs/report_generation_refactor_assessment.md` (R1 + R2 + R3 approved;
R4 module-split and R5 dashboard-decompose deferred).
Scope: `tc_rerender_report.py`, `tarean/assets/tidecluster.{js,css}`,
`tests/test_report_exports.py`. **v2 report only; legacy/v1 untouched.**

## Goal

Every column header, column description, tier meaning, and reusable explanatory
blurb has **one definition**, and the HTML legend / tooltips / `columns.tsv` /
TSV bundle / embedded JSON are all **derived** from it — so the drift that bit us
during the HOR-tier work (a concept edited in up to six places) cannot recur.

## Format decision (the "least friction" call)

- **Column + tier descriptions → structured Python** (extend the working
  `_PEAK_COLUMNS` pattern). They pair with `fmt`/`source`, feed the
  machine-readable `columns.tsv`, and are emitted into both server-side HTML and
  the JS `#tc-coldict`. Markdown can't carry the tuple shape → structured is
  lower friction and consistent.
- **Page-level prose stays markdown** (`docs/report_content/*.md` — already good).
- **Inline explanatory blurbs that interpolate live counts** (per-TRC summary,
  callouts) → **centralized Python format-strings** (a `BLURBS` table). Markdown
  can't interpolate the dynamic numbers cleanly, so a `.format(**vals)` table in
  one place is the lowest-friction home. Pure-static blurbs may instead become a
  short `report_content/*.md` snippet if they carry no values.

Net: nothing new added to the dependency surface; everything stays stdlib.

---

## Target data structures (the single sources)

Add one cohesive block near the existing `_PEAK_COLUMNS` (~`tc_rerender_report.py:2673`):

```python
from collections import namedtuple
ColSpec = namedtuple("ColSpec", "key header fmt source desc")
#   key    machine key (JSON / internal)
#   header HTML/TSV column label (HTML-facing name)
#   fmt    int | f3 | f4 | yn | raw | bool | str   (display formatting)
#   source model/array field name, or kitehor raw column, or None if computed
#   desc   one-line description (drives columns.tsv + legend + <th title>)

# (1) Export/data tables — atomic columns. Replaces _TRC_TABLE_COLS,
#     _TRA_TABLE_COLS, _PEAK_COLUMNS (the peaks list is unchanged in content).
COLUMN_DICT = {
    "trc_table": [ColSpec(...), ...],
    "tra_table": [ColSpec(...), ...],
    "tra_peaks": [ColSpec(...), ...],   # == today's _PEAK_COLUMNS, re-shaped
}

# (2) HOR-order tier — defined ONCE; referenced by the ×k legend, the
#     hor_order_confidence column desc, and (label only) the JS payload.
TIER_DEFS = {
    "strict":    ("HOR (strict)",    "clean integer divisor, low k — confident order"),
    "supported": ("HOR (supported)", "order backed by cross-array founder consensus"),
    "weak":      ("founder recovered (HOR ≈)",
                  "founder reliable; order approximate — irregular, very high k, or relaxed rescue"),
    "none":      ("—", "no higher-order structure (multiplicity 1)"),
}
# _HOR_TIER_LABEL becomes {k: v[0] for k,v in TIER_DEFS.items()} (back-compat).

# (3) The VISIBLE array-table columns (composite cells: Founder/Strongest,
#     Δid, ×k, Other periods, Subrepeat, SSR) — structured replacement for the
#     hardcoded <dl> in arrays_legend(). Descriptions may interpolate TIER_DEFS.
ARRAY_TABLE_LEGEND = [
    ("Founder / Strongest", "<...>"),
    ("Δid",                 "<...>"),
    ("×k",                  "<...uses TIER_DEFS...>"),
    ("Other periods",       "<...>"),
    ("Subrepeat",           "<...>"),
    ("SSR",                 "<...>"),
]

# (4) Reusable explanatory blurbs with {named} placeholders for live values.
BLURBS = {
    "trc_structural_summary": "<p class='tc-note'>... {tot_hor} ... {tot_subrep} ...</p>",
    "ssr_callout":            "...",
    "fallback_callout":       "...",
    # one entry per inline <p> currently embedded in render_* code
}
```

### What each replaces

| New | Replaces (current) | Consumers after |
|---|---|---|
| `COLUMN_DICT["tra_peaks"]` | `_PEAK_COLUMNS` (2673), `_PEAK_COL_KEYS/HEADERS`, `_DETAILS_COL_DESC` (3118) | unfold headers/tooltips, `#tc-coldict`, peak JSON, `tra_peaks.tsv` |
| `COLUMN_DICT["tra_table"]` | `_TRA_TABLE_COLS` (2973) | `tra_table.tsv`, `columns.tsv` |
| `COLUMN_DICT["trc_table"]` | `_TRC_TABLE_COLS` (2952) | `trc_table.tsv`, `columns.tsv` |
| `TIER_DEFS` | `_HOR_TIER_LABEL` (82) + inline ×k prose in `arrays_legend` (3147–3156) + `hor_order_confidence` desc in `_TRA_TABLE_COLS` | `arrays_legend`, `columns.tsv`, JS label |
| `ARRAY_TABLE_LEGEND` | hardcoded `<dl>` in `arrays_legend()` (3137–3168) | `arrays_legend()` (generated) |
| `BLURBS` | ~15 inline `<p style=…>` blocks (2126–2138, 1607, 1762, 2515, callouts) | the render_* sites, via `BLURBS[k].format(...)` |

---

## Work items

### R1a — unify the export/dictionary columns
- Define `ColSpec` + `COLUMN_DICT`. Re-express `_PEAK_COLUMNS` as
  `COLUMN_DICT["tra_peaks"]` (same content, namedtuple shape); keep
  `_PEAK_COLUMNS = COLUMN_DICT["tra_peaks"]` and
  `_PEAK_COL_KEYS/HEADERS` as thin aliases so existing call sites
  (`_peak_cell_values`, `_fmt_peak_cell`, `_peaks_payload_script`,
  `write_table_exports`) keep working with minimal edits (switch `c[0]`→`c.key`,
  `c[1]`→`c.header`, `c[3]`→`c.fmt`, `c[4]`→`c.desc`).
- `write_table_exports` reads `COLUMN_DICT[table]` for headers, values
  (`getattr`/`a.get(c.source)`), and `columns.tsv` — one loop per table instead
  of three bespoke blocks. `_tsv_cell` + `_fmt_peak_cell` stay.

### R1b — generate the HTML legend from data
- `arrays_legend()` builds its three sections from `ARRAY_TABLE_LEGEND` (visible
  columns), `_TIER_DESC` (subrepeat tiers — already structured), and
  `COLUMN_DICT["tra_peaks"]` (peak columns) — **no hardcoded column prose left**.
- `#tc-coldict` (in `_peaks_payload_script`) emitted from
  `COLUMN_DICT["tra_peaks"]` (already is, via `_PEAK_COLUMNS`).

### R1c — TIER_DEFS single source
- `_HOR_TIER_LABEL` derives from `TIER_DEFS`. The ×k legend entry in
  `ARRAY_TABLE_LEGEND` and the `tra_table` `hor_order_confidence` `desc` both
  reference `TIER_DEFS` text. JS keeps receiving the label via payload (`htl`) —
  no JS change needed for tiers.

### R2 — centralize blurbs
- Create `BLURBS`; move each inline explanatory `<p>` into it as a format-string;
  replace the render-site literal with `BLURBS[key].format(**vals)`. Pure-static
  ones with no `{...}` may instead go to a `report_content/` snippet (decide
  per-blurb; default to `BLURBS` for friction).

### R3 — `.tc-note` CSS
- Add `.tc-note { font-size:12px; color:var(--fg-muted); margin:4px 0 10px; }`
  to `tidecluster.css`; replace the ~6 inline
  `style="font-size:12px;color:var(--fg-muted)"` occurrences (and the BLURBS
  templates) with `class="tc-note"`.

---

## Validation

1. **Behaviour-preserving (the key invariant):** rerender the slyco `run_e2e`
   **and** the arabidopsis `run_e2e` before/after and assert:
   - `data/*.tsv` **headers byte-identical** (column set/order unchanged);
   - `columns.tsv` row set unchanged except intended description rewordings;
   - each TRC dashboard still contains all legend sections + the embedded
     `tc-peaks`/`tc-coldict`, and `data-details` count stays 0;
   - dashboard sizes unchanged (±, this isn't a perf change).
2. **New test (the anti-drift guard):** in `tests/test_report_exports.py`, assert
   **every rendered table header has a `COLUMN_DICT` entry and vice-versa**, and
   that `arrays_legend()` references no column label absent from the dict — so a
   future column added to the HTML but not the dict (or vice-versa) fails CI.
3. **Existing suite green:** `tests/unit.sh` (incl. the size guardrail and the
   tier test) unchanged-green.
4. **Spot visual:** open one slyco + one arabidopsis dashboard from `file://`,
   confirm the legend + unfold tooltips read identically to today.

## Risk & rollback
- **Risk: low–medium.** Mostly mechanical re-pointing of accessors; the content
  is preserved. The one content change is *intentional* de-duplication of tier
  wording (single canonical text). Covered by validation #1–#2.
- **Rollback:** the change is self-contained in `tc_rerender_report.py` (+ small
  css); revert the commit. No data-format change → no rerun needed for old
  reports.
- **Sequencing:** land R1a → R1b → R1c → R2 → R3 as **separate commits** (each
  green on the suite), so a regression bisects cleanly.

## Out of scope (deferred, per decision)
- R4 (split into `report_model/render/export.py`).
- R5 (decompose `render_trc_dashboard`).
- Any legacy/v1 report code.

## Open questions for review
1. OK to **canonicalise the HOR-tier wording** to one text (some current strings
   differ slightly between the ×k legend and the export desc)? Recommended.
2. For R2, default each blurb to **`BLURBS` (Python)**; only move a blurb to a
   `report_content/*.md` snippet if it's fully static. Agree, or prefer all-Python
   for consistency?
3. Add the **"legend ≡ COLUMN_DICT" CI assertion** (validation #2)? Recommended —
   it's what makes the single-source guarantee stick.
