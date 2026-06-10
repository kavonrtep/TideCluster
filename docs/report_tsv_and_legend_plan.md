# Implementation plan — per-TRA legend at the unfold + a single accessible TSV

Date: 2026-06-10
Status: IMPLEMENTED (all decisions below locked + built). Parts 1+2+3 landed
together: TSV bundle in `<prefix>_report/data/` (trc/tra/tra_peaks/columns),
lazy child rows from one embedded per-TRC `tc-peaks` JSON (file://-safe),
shared in-row column legend + `<th>` tooltips, size guardrail in
`tests/test_report_exports.py`. Measured: TRC_4.html 5.3 MB → 1.3 MB.
Scope: HTML report (`tc_rerender_report.py`) + a consolidated TSV export.
Two requests + one cross-cutting constraint:
1. The per-TRA table shown **on unfold** (the Details child row, derived from
   kitehor) has cryptic columns; surface a legend **there**.
2. The per-TRC / per-TRA tables should be **downloadable as a single accessible
   TSV** instead of being scattered across several `*_kite/` files.
3. **Browser weight.** In satellite-rich genomes these tables get large and the
   report is already sluggish when unfolding rows. Nothing here may make a TRC
   dashboard heavier; ideally it gets much lighter. **This reshapes Part 1** (the
   in-row legend must be shared, not duplicated per row) and adds **Part 3**.

Decisions locked (Petr, 2026-06-10): the 4-file **bundle**, **HTML column
names** (+ dictionary), in **`<prefix>_report/data/`**, **tooltips + in-row
legend** for Part 1, keep `monomer_size_top3_estimats.csv`.

---

## Current state (measured)

The HTML report has **three table granularities**:

| HTML table | granularity | rows | current data source |
|---|---|---|---|
| TRC overview (index / genome page) | per **TRC** | 1/TRC | `report.json` model aggregates (kite block, TAREAN, superfamily) |
| "Tandem repeat arrays" table (TRC dashboard) | per **TRA** | 1/array | `*_kite/monomer_size_top3_estimats.csv` (55 cols) |
| Details child row (▸ unfold) | per **rescored peak** | N/array | `*_kite/kitehor.rescored.peaks.tsv` (30 cols) |

Scattered TSV/CSV sources today (all under `<prefix>_kite/`, raw-kitehor names):
`kitehor.kite.tsv`, `kitehor.kite.peaks.tsv`, `kitehor.rescored.peaks.tsv`,
`kitehor.ssr.tsv`, `kitehor.ssr.regions.tsv`, `kitehor.tandem_validate.tsv`,
`kitehor.verdicts.tsv`, plus the TideCluster-derived
`monomer_size_top3_estimats.csv`. None of them is a 1:1 mirror of an HTML table,
the names are kitehor-internal, and the per-peak file uses raw column names
(`scan_occupancy_frac`, `kmer_phase_contrast`, …) that the HTML renames to terse
abbreviations (`occ`, `phaseC`, …).

**A legend already exists.** `arrays_legend()` (rendered once per TRC dashboard,
just above the arrays table) documents the array-table columns, the
subrepeat-evidence tiers, AND every Details-row peak column
(`_DETAILS_COL_DESC`). The problem is **placement**: it is a *collapsed*
`<details>` at the top of the dashboard, and the peak-column glossary is the last
block inside it — so when a user expands a single TRA row lower down, the inner
peak table appears with no nearby explanation.

**`report.json` is the single source of truth** the HTML renders from (the full
model, 23 MB on slyco). Anything we export should be derived from the same model
so it cannot drift from the HTML.

### Performance findings (measured, slyco `run_e2e`)

- `trc/TRC_4.html` = **5.3 MB**, of which **~5.0 MB (93 %) is the inline
  `data-details`** child-row HTML: 369 arrays × ~13.6 KB of pre-rendered,
  HTML-escaped peak-table markup each, embedded in the `<tbody>` **upfront**.
- `deferRender: true` is already set (`tidecluster.js`), but it only defers
  *row-node creation for paging* — the full `<tbody>` (all `data-details`) is
  still parsed and the strings held in the DataTables data model at init. So the
  cost is paid even though only ~25 rows show and most are never unfolded.
- Extrapolation: a centromeric satellite TRC with a few thousand arrays → tens of
  MB on one page = the sluggishness Petr observes.
- The visible array table itself (369 rows × ~12 small cells) is cheap; **the
  child-row HTML is essentially the entire weight.**

Conclusion: the dominant lever is **not** to embed per-peak child HTML inline.
This is Part 3, and it dovetails with Part 2 (the per-peak data we export is the
same data the child row needs).

---

## Part 1 — legend at the unfold (request 1)

Goal: a user expanding a TRA row can read what `occ / scan_n / cov / spat / autoF
/ phaseC / sub / phtm / tier` mean **without scrolling away**.

### Options

- **1A — `<th title="…">` tooltips on the inner peak table.** Add the existing
  `_DETAILS_COL_DESC` text as a `title=` on each `<th>` of the Details peak table
  (hover to read). Tiny change, zero layout cost. Con: hidden until hover, not
  discoverable on touch / when printing.
- **1B — a compact one-line glossary inside each Details child row.** Append a
  small `<details class="tc-col-legend"><summary>column meanings</summary>…`
  right under the peak table, built from the same `_DETAILS_COL_DESC`. Visible,
  discoverable, self-contained per row. Con: repeats the same block in every
  child row's HTML (size — but child-row HTML is already stored per row in
  `data-details`, so it is N× regardless).
- **1C — both: tooltips on headers + one shared, always-near legend.** `title=`
  tooltips (1A) for quick hover, plus keep the full glossary but move/duplicate a
  **peak-column** mini-legend so it renders adjacent to the table. 
- **1D — full-name headers, drop the abbreviations.** Rename the inner columns to
  readable names (`occupancy`, `phase-contrast`, …) so the legend is barely
  needed. Con: wider table, breaks the dense diagnostic layout reviewers know.

**Recommendation: 1A + 1B, but rendered ONCE and shared (not per row).** Given
the performance findings, the in-row legend must **not** be a `<dl>` baked into
every row's payload (that multiplies the dominant cost). Instead:

- The child-row column **tooltips** (`title=` on each `<th>`) are built when the
  child row is rendered (in JS under Part 3), from a single column-dictionary
  object embedded once per page — cost is per *opened* row only.
- The "column meanings" **glossary** is a **single hidden block** rendered once
  in the dashboard; opening any Details row reveals/links to it (e.g. a small
  "ⓘ column meanings" affordance in the child row toggles the one shared block).
  One copy per page, not one per array.

Both still read from the single column-dictionary source (Part 2's `COLUMN_DICT`,
which also feeds `columns.tsv` and the HTML legend). Keep the existing top-of-page
`arrays_legend()`.

Touch points: `tidecluster.js` child-row builder (tooltips + reveal the shared
block), `tc_rerender_report.py` (emit the shared block + the column dictionary
once), a little CSS. **No per-row legend HTML.**

---

## Part 2 — a single accessible TSV mirroring the HTML tables (request 2)

The hard constraint is the **granularity mismatch**: TRC (1/TRC), TRA (1/array),
and peak (N/array) cannot share one flat table without either denormalizing
(repeating TRC fields on every TRA row) or going long/EAV (ugly). So "a single
TSV" needs a decision about *which* granularity is canonical.

### Design options

- **2A — One per-TRA TSV (recommended primary).**
  `<prefix>_report/data/tra_table.tsv`: **one row per array**, mirroring the
  visible "Tandem repeat arrays" table for *all* TRCs at once, with `TRC_ID` as
  the first column and the parent TRC's summary fields (type, TRC monomer,
  TAREAN, superfamily) **joined on** so each row is self-contained. A user opens
  one file, filters by `TRC_ID` to get a TRC's array table, or pivots on
  `TRC_ID` to recover the TRC overview. This is essentially today's
  `monomer_size_top3_estimats.csv` **relocated, renamed, TRC-enriched, and
  column-aligned to the HTML** (HTML-facing names + a documented header).
  The per-peak unfold detail is inherently a different shape → stays a **second**
  file (see 2B), but both live in one obvious place with a shared dictionary.

- **2B — A small documented bundle (recommended companion).**
  `<prefix>_report/data/`:
  - `trc_table.tsv`   — per-TRC overview (mirrors the index table)
  - `tra_table.tsv`   — per-TRA (mirrors the arrays table) *(= 2A)*
  - `tra_peaks.tsv`   — per-rescored-peak (mirrors the unfold table), HTML column
    names, with `TRC_ID/seqid/start/end` keys so it joins to `tra_table.tsv`
  - `columns.tsv`     — machine-readable data dictionary (column → table →
    description), generated from the **same** `_DETAILS_COL_DESC` / array-column
    descriptions that feed the HTML legend (single source of truth)
  All derived from the `report.json` model in one pass → guaranteed to mirror the
  HTML.

- **2C — One denormalized "long" TSV with a `row_type` column** (`TRC` / `TRA` /
  `PEAK`) and a union of columns. Rejected: sparse, hard to read in a spreadsheet,
  and no real analytical advantage over 2B.

**Recommendation: 2B (the bundle), with `tra_table.tsv` as the headline file
(2A).** It honestly matches the three HTML granularities, every file is a clean
rectangle openable in Excel/R/pandas, and the four files share keys so they
join. If Petr wants *literally one file*, fall back to 2A alone (`tra_table.tsv`
denormalized) and keep `tra_peaks.tsv` as the only extra, documented as "the
unfold detail".

### Where + how (architecture)

- **Emit from the report renderer.** `tc_rerender_report.py` already builds the
  full `model` and is independently re-runnable (`--input-dir`). Add a
  `write_table_exports(model, out_dir)` that writes the bundle into
  `<prefix>_report/data/` next to `report.json`, from the same model dicts the
  HTML tables iterate. Re-rendering regenerates both HTML and TSV together → they
  can never drift.
- **Single column-description source.** Promote the column descriptions to one
  structure (`COLUMN_DICT = {table: [(col, desc), …]}`) that feeds (a) the HTML
  legend, (b) the `title=` tooltips (Part 1), and (c) `columns.tsv`. The current
  `_DETAILS_COL_DESC` + the inline array-table legend `<dl>` become entries in
  it.
- **Naming / discoverability.** Link the files from the report (e.g. a
  "Download tables (TSV)" line on the index + each TRC dashboard pointing at
  `data/tra_table.tsv` etc.). Keep the raw `*_kite/kitehor.*` files as-is (they
  remain the provenance / debug source); the new bundle is the *curated, HTML-
  mirroring* export.
- **`run_all` vs rerender.** The exports are written wherever the report is
  written, so both `TideCluster.py` (which calls the renderer) and a standalone
  `tc_rerender_report.py` produce them with no extra wiring.

---

## Part 3 — lighten the report (the actual performance fix)

Stop embedding per-peak child HTML inline; build it lazily from compact data.

### Constraint: the report is viewed from `file://`

Reports are shared as a folder and opened locally. Chrome blocks `fetch()` of
local files (`file://` CORS), which is *why* everything is inlined today. So the
fix must keep the page **self-contained** — no runtime `fetch()`.

### Design (recommended)

Replace the N× escaped `data-details` attributes with **one embedded JSON block
per dashboard** + **lazy child-row rendering in JS**:

- `tc_rerender_report.py`: instead of writing pre-rendered child HTML into each
  `<tr data-details="…">`, emit a single
  `<script type="application/json" id="tc-peaks">{…}</script>` on the dashboard
  holding the per-array peak rows for *this TRC only* (keyed by
  `seqid:start:end`), columns = exactly what the unfold table shows. Each `<tr>`
  keeps only a tiny `data-array-key` instead of 13 KB of HTML.
- `tidecluster.js`: on first unfold of a row, look up its peaks in the parsed
  JSON, build the table DOM (with `<th title=…>` tooltips), and show it via
  `row.child()`. Build the shared column-legend block once, lazily.

Why this is much lighter:
- **Data, not markup:** JSON peak rows carry no repeated tags / CSS classes /
  tier-pill HTML — a fraction of the 13.6 KB/row escaped markup.
- **Parsed once** as JSON vs the browser materialising 369 escaped HTML
  sub-documents into the DataTables model.
- **DOM holds only opened rows' detail**, built on demand.
- Still 100 % self-contained (embedded `<script>`, no `fetch`).

Expected: `TRC_4.html` ~5.3 MB → a few hundred KB; the cost now scales with
*opened* rows, not total arrays.

### Alternatives (for the record)

- **3-alt-a — external per-TRC JSON loaded via `fetch()`** (`data/peaks/TRC_x.json`):
  smallest HTML, but breaks on `file://` → rejected as the default. Could be a
  future opt-in for server-hosted reports.
- **3-alt-b — keep inline `data-details` but cap/scroll** (only embed details for
  the first K arrays, link the rest to the TSV): a stop-gap, doesn't fix the big-
  TRC case → rejected.
- **3-alt-c — DataTables `Scroller`/server-side**: heavy dependency / complexity
  for a static report → rejected.

### Note (out of scope, flag only)

`report/data/report.json` is **23 MB**. Confirm no HTML page loads it at runtime
(it appears to be a data sidecar, not fetched by the pages). If something does
load it, that is a separate, larger heaviness to tackle — not part of this work.

## Validation plan

1. **Mirror check (automated):** for the slyco `run_e2e` model, assert
   `tra_table.tsv` has exactly one row per array and its founder/strongest/×k/
   tier values equal the per-array model values the HTML renders; `trc_table.tsv`
   row count == number of TRC dashboards; `tra_peaks.tsv` row count == Σ rescored
   peaks. Add as a unit/integration check.
2. **Join check:** `tra_peaks.tsv ⨝ tra_table.tsv` on `TRC_ID/seqid/start/end`
   loses no rows; `tra_table.tsv` pivoted on `TRC_ID` reproduces `trc_table.tsv`
   counts.
3. **Dictionary check:** every column in each exported TSV appears in
   `columns.tsv`, and every `_DETAILS_COL_DESC` entry is used by the HTML legend
   (no orphans) — guarantees legend ≡ export ≡ HTML.
4. **Legend-at-unfold (visual):** rerender `run_e2e`; confirm a Details child row
   shows header tooltips + the shared column legend (one block per page), and the
   values match `tra_peaks.tsv` for that array.
5. **Page-weight regression (the point of Part 3):** assert each TRC dashboard
   HTML shrinks substantially — e.g. `TRC_4.html` from ~5.3 MB to < ~0.5 MB — and
   that the embedded `tc-peaks` JSON for a TRC has exactly one entry per array,
   matching `tra_peaks.tsv`. Add a size assertion so the inline-HTML bloat can't
   silently return.
6. **`file://` smoke:** open a rerendered dashboard directly from disk (no
   server) and confirm unfold still builds the child row (no `fetch`).
7. Confirm no HTML page `fetch`/embeds `report.json` at runtime (else flag
   separately).

---

## Decisions (locked) & remaining questions

**Locked (Petr, 2026-06-10):**
- TSV = the **4-file bundle** (`trc_table` / `tra_table` / `tra_peaks` /
  `columns`) in **`<prefix>_report/data/`**.
- Export uses **HTML column names**; `columns.tsv` carries name + description
  (and may also list the raw kitehor source name).
- Part 1 = **tooltips + in-row legend**, implemented **shared/once** (Part 3
  constraint), not duplicated per row.
- **Keep** `monomer_size_top3_estimats.csv` (read by `extract_consensus.R`);
  `tra_table.tsv` is the HTML-mirroring superset.
- Report must get **lighter**, not heavier → Part 3 (lazy child rows via one
  embedded per-TRC JSON block; self-contained, `file://`-safe).

**Still to confirm before/while building:**
1. **Part 3 scope in this change.** Part 3 is a `tidecluster.js` + renderer
   refactor of the child-row mechanism. Do it **together** with Parts 1–2 (Part 1
   depends on it for the shared-legend approach), or land Parts 1–2 first with a
   minimal interim and Part 3 as an immediate follow-up? Recommend **together** —
   Part 1's "no per-row bloat" requirement is only cleanly met once Part 3 moves
   the child row to JS.
2. **`tra_peaks.tsv` columns = exactly the unfold table**, or the fuller raw
   rescored-peaks set (all 30 cols)? Recommend: the unfold-table columns by
   default (true mirror), with the raw file still available under `*_kite/` for
   power users.
3. **Embedded `tc-peaks` JSON granularity:** per-TRC dashboard (recommended,
   small, only that TRC's arrays) vs one global blob. Per-TRC keeps each page
   light and is the natural fit since dashboards are per-TRC pages.
4. **Size-regression guardrail:** OK to add an automated assertion that a TRC
   dashboard stays under, say, a few hundred KB + ~few KB/array, so inline bloat
   can't silently return?
