# HTML report generation — maintainability assessment & refactoring recommendations

Date: 2026-06-10
Status: assessment for review (no code changed)
Scope: the **v2 report generator** (`tc_rerender_report.py`, 3680 lines) +
`tarean/assets/tidecluster.{js,css}` + `docs/report_content/*.md`. **Excludes**
legacy/v1 reports.
Trigger: descriptive/legend text had to be edited in several places by hand
during the HOR-tier + TSV-bundle work; this maps why and proposes a fix.

---

## 1. How descriptive text is handled today — three parallel mechanisms

| Mechanism | Holds | Where | Verdict |
|---|---|---|---|
| **External markdown** | page-level prose (editable, copied to output) | `docs/report_content/{overview,kite,tarean,superfamilies,credits}.md`, loaded via `load_section()` (`tc_rerender_report.py:262`), cached, shipped by `copy_report_content_source()` | **Good** — non-developers can edit; clean separation |
| **Python data constants** | structured column/tier descriptions | `_PEAK_COLUMNS` (2673), `_TIER_DESC` (3091), `_TIER_SHORT` (2542), `_TRC_TABLE_COLS` (2952), `_TRA_TABLE_COLS` (2973), `COMBINED_CLASS_LABELS` (329), `_HOR_TIER_LABEL` (82) | **Mixed** — some are single-source, some duplicate each other |
| **Inline HTML f-strings** | legend prose + explanatory blurbs baked into render code | `arrays_legend()` array-table `<dl>` (3137–3168); ~15 inline `<p style="…">` blurbs (e.g. 2126–2138, 1607, 1762, 2515) | **Problem** — high-churn copy living inside render logic |

The page-level prose is well externalized. The **dense, high-churn content —
column legends and concept descriptions — is split across constants *and*
hardcoded HTML**, with no single source per concept.

---

## 2. Core problem: the same concept is described in 3–6 places (drift risk)

### Example A — array-table columns (Founder, ×k, Subrepeat, SSR)

The *same columns* are described twice, independently, with no shared source:

- **Legend prose** is hardcoded HTML inside `arrays_legend()`
  (`tc_rerender_report.py:3137–3168`).
- **Data-dictionary descriptions** live in `_TRA_TABLE_COLS` (2973), used **only**
  to emit `columns.tsv`.

Edit the meaning of `multiplicity`/`×k` and you must change both, by hand, or the
HTML legend and the downloadable `columns.tsv` disagree.

### Example B — the HOR-order tier (worst case: six locations)

`strict | supported | weak` and its meaning are encoded in:

1. `arrays_legend()` inline HTML — the `×k` `<dd>` (3147–3156)
2. `_HOR_TIER_LABEL` (Python, 82)
3. `tidecluster.js` — a JS mirror of `_HOR_TIER_LABEL` in the child-row builder
4. `_TRA_TABLE_COLS` — the `hor_order_confidence` row description (2973+)
5. `columns.tsv` (generated from #4)
6. README "Concepts" block

This session required hand-syncing several of these. That is the smell the
refactor should remove.

### Contrast — the pattern that already works

`_PEAK_COLUMNS` (2673) is a **single source** feeding (a) the unfold table
headers, (b) the JS `<th>` tooltips + in-row legend (via the embedded
`#tc-coldict`), (c) the embedded peak JSON, and (d) `tra_peaks.tsv`. One edit
propagates everywhere. **The fix is to apply this pattern to the other tables.**

---

## 3. Secondary issues

- **One 3680-line module, ~6 responsibilities**: markdown converter
  (`md_to_html`, 118), path/model loading (`resolve_paths`, `build_model`,
  `load_*`), domain classification (`classify_peak_tier`, `_hor_tier`), HTML
  rendering (`page_shell` + `render_*`), asset copying, TSV export. Cohesive but
  hard to navigate.
- **`render_trc_dashboard` is 293 lines** — assembles many sections via nested
  f-strings; the single function most error-prone to modify.
- **Repeated inline styling**: `style="font-size:12px;color:var(--fg-muted)"`
  appears ~6× (e.g. 1607, 1762, 2116, 2126, 2164, 2515) instead of one CSS class.
- **Constraint to respect**: rendering is pure-stdlib f-strings *by design* — no
  templating dependency (see `md_to_html`'s own note, ~114: "extend the parser
  rather than reach for a dep"). A Jinja-style rewrite is **out of scope**;
  recommendations below stay within the stdlib ethos.

---

## 4. Recommendations (ordered: value ÷ risk)

### R1 — one column dictionary as the single source for every table *(highest value, well-scoped)*

Generalise the `_PEAK_COLUMNS` pattern. Define one structure, e.g.

```
COLUMN_DICT = {
  "trc_table":  [ColSpec(key, header, fmt, description, source), ...],
  "tra_table":  [...],
  "tra_peaks":  [...],   # = today's _PEAK_COLUMNS
}
```

and **derive from it**:
- the HTML legend (`arrays_legend()` generated from data, not hardcoded prose),
- the array/TRC/peak table `<th>` headers,
- the `<th title=…>` tooltips + the JS in-row legend (already via `#tc-coldict`),
- the `columns.tsv` export and the TSV column order,
- (optionally) the README concept table.

Folds `_TRC_TABLE_COLS`, `_TRA_TABLE_COLS`, `_DETAILS_COL_DESC`, and the inline
`arrays_legend()` array-table `<dl>` into one place → removes Example-A and most
of Example-B drift. **Effort: medium. Risk: low** (additive; validated by the
existing mirror/guardrail tests + a new "legend ≡ dictionary" assertion).

### R2 — one home for concept blurbs *(medium value, low risk)*

Move the recurring inline prose (HOR-order tier wording, the per-TRC structural
summary paragraph at 2126–2138, the callouts) into either the column dictionary
(structured) or short `docs/report_content/` snippets, so wording changes happen
in one place across Python + JS + README. **Effort: low–medium. Risk: low.**

### R3 — `.tc-note` CSS class *(low value, trivial)*

Replace the ~6 repeated inline `style="font-size:12px;color:var(--fg-muted)"`
with a single `.tc-note` class in `tidecluster.css`. **Effort: trivial.**

### R4 — split the module by responsibility *(optional, larger)*

Only if the file keeps growing. Candidate split (all pure-stdlib):
`report_model.py` (paths + `build_model` + `load_*`), `report_render.py`
(`page_shell` + `render_*` + helpers), `report_export.py` (`write_table_exports`),
shared `report_common.py` (`md_to_html`, `esc`, `COLUMN_DICT`, tier maps).
Improves navigability; **Effort: high. Risk: medium** (import churn, no behaviour
change). Defer unless R1–R3 aren't enough.

### R5 — decompose `render_trc_dashboard` *(optional)*

Extract its section builders (arrays section, distribution, KITE figure,
superfamily) into named helpers so the 293-line function becomes an assembler.
**Effort: medium. Risk: low–medium.** Nice-to-have alongside R4.

---

## 5. Suggested sequencing

1. **R1** (column dictionary) — biggest maintainability win, removes the drift
   that bit us this session, natural extension of `_PEAK_COLUMNS`.
2. **R3** (`.tc-note`) — trivial, bundle with R1.
3. **R2** (blurb home) — once R1 establishes the structured-text pattern.
4. **R4 / R5** — only if the module keeps growing; revisit after R1–R2.

Each step is independently shippable and covered by the existing report
unit/guardrail tests (`tests/test_report_exports.py`), to which R1 should add a
"every rendered table header/description is in COLUMN_DICT" assertion.

---

## 6. Open questions for review

1. Start with **R1 only**, or bundle **R1 + R2 + R3** as one "report content
   single-source" change?
2. For R1, keep the array/TRC table column **descriptions structured in Python**
   (consistent with `_PEAK_COLUMNS`), or push them into markdown like the page
   prose? (Recommend structured — they pair with `fmt`/`source` and feed
   `columns.tsv`.)
3. Is the **module split (R4)** wanted now, or parked until the file forces it?
