#!/usr/bin/env python3
"""Build report v2 (modernised HTML) from an existing TideCluster output dir.

Reads the on-disk TSV/CSV/GFF3/JSON/PNG artefacts produced by a previous
run and emits <prefix>_report_v2/ alongside, containing:

    data/report.json     unified data model (this file is what the HTML
                         pages bind to via an inlined <script> const)
    assets/              vendored DataTables + jQuery + TideCluster CSS/JS
    index.html           landing page
    trcs.html            all-TRCs table
    tarean.html          TAREAN summary
    kite.html            KITE summary + HOR aggregates
    superfamilies.html   superfamily listing
    trc/TRC_<N>.html     per-TRC dashboards

The tool reads only; the original <prefix>_*.html files are untouched
and continue to work.

This phase (phase 3 of the implementation plan) emits the data model
and copies the assets. HTML rendering lands in phases 4 and 5.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import html
import json
import os
import re
import shutil
import statistics
import sys
from collections import Counter
from pathlib import Path

__version__ = "2"  # schema version for report.json (v2: 4-bin HOR)

# ----------------------------------------------------------------------
# HOR classification — mirrors tarean/kite.R. See
# docs/hor_classification.md for the algorithm and rationale.
# ----------------------------------------------------------------------
HOR_TOL             = 0.10
HOR_HARMONIC_BONUS  = 0.5
HOR_BIN_WEAK        = 0.10
HOR_BIN_MODERATE    = 0.20
HOR_BIN_STRONG      = 0.40
HOR_MAX_K           = 5
HOR_GRID_STEP       = 1


def _score_m_star(ms, ss, m_star):
    if m_star is None or m_star <= 0: return None
    total = sum(s for s in ss if s is not None and s > 0)
    if total <= 0: return None
    ks, errs, close = [None] * len(ms), [None] * len(ms), [0.0] * len(ms)
    for i, m in enumerate(ms):
        if m is None: continue
        k = round(m / m_star)
        if k < 1: return None
        ks[i] = k
        errs[i] = abs(m - k * m_star) / m_star
        close[i] = max(0.0, 1.0 - errs[i] / HOR_TOL)
    if not any(k == 1 for k in ks if k is not None): return None
    if not any(k is not None and k >= 2 for k in ks): return None
    base_w = sum(ss[i] * close[i] for i in range(len(ms))
                 if ks[i] == 1 and ss[i] is not None)
    harm_w = sum(ss[i] * close[i] for i in range(len(ms))
                 if ks[i] is not None and ks[i] >= 2 and ss[i] is not None)
    f_base, f_harm = base_w / total, harm_w / total
    distinct = len({ks[i] for i in range(len(ms))
                    if ks[i] is not None and ks[i] >= 2 and close[i] > 0})
    bonus = 1.0 + HOR_HARMONIC_BONUS * max(0, distinct - 1)
    conf = ((max(0.0, f_base) * max(0.0, f_harm)) ** 0.5) * bonus
    contributing = [ks[i] for i in range(len(ms))
                    if ks[i] is not None and ks[i] >= 2 and close[i] > 0]
    k_max = max(contributing) if contributing else None
    return {"m_star": float(m_star), "ks": ks, "confidence": conf,
            "distinct": distinct, "k_max": k_max}


def _hor_candidates(ms):
    valid = [m for m in ms if m is not None and m > 0]
    if not valid: return []
    out = {round(m) for m in valid}
    for k in range(2, HOR_MAX_K + 1):
        out.update(round(m / k) for m in valid)
    lo = max(5, int(min(valid) * 0.5))
    hi = int(max(valid) * 0.55)
    if hi > lo:
        out.update(range(lo, hi + 1, HOR_GRID_STEP))
    return [x for x in out if x > 0]


def compute_hor(m1, m2, m3, s1, s2, s3):
    """Return dict of HOR fields for one array (mirrors kite.R).

    Fields: status, confidence, base_monomer, hor_period, n_harmonics."""
    ms, ss = [m1, m2, m3], [s1, s2, s3]
    best = None
    for m_star in _hor_candidates(ms):
        fit = _score_m_star(ms, ss, m_star)
        if fit and (best is None or fit["confidence"] > best["confidence"]):
            best = fit
    conf = best["confidence"] if best else 0.0
    if conf < HOR_BIN_WEAK:     status = "No HOR"
    elif conf < HOR_BIN_MODERATE: status = "HOR weak"
    elif conf < HOR_BIN_STRONG:   status = "HOR moderate"
    else:                         status = "HOR strong"
    if best:
        base_monomer = int(round(best["m_star"]))
        hor_period   = (int(round(best["k_max"] * best["m_star"]))
                        if best["k_max"] else None)
        n_harmonics  = int(best["distinct"])
    else:
        base_monomer = hor_period = None
        n_harmonics  = 0
    return {"hor_status": status, "hor_confidence": round(conf, 4),
            "hor_base_monomer": base_monomer, "hor_hor_period": hor_period,
            "hor_n_harmonics": n_harmonics}

# ----------------------------------------------------------------------
# Minimal Markdown → HTML converter (stdlib only).
# ----------------------------------------------------------------------
# Supported subset:
#   - ATX headings:      `# ` .. `###### ` at start of line
#   - Paragraphs:        blank-line separated
#   - Unordered lists:   lines starting with `- ` or `* `
#   - Ordered lists:     lines starting with `<n>. `
#   - Fenced code:       ```` ```lang ```` ... ```` ``` ````
#   - Bold inline:       `**text**`
#   - Italic inline:     `*text*`  (not adjacent to * or whitespace)
#   - Inline code:       `` `text` ``  (contents HTML-escaped)
#   - Links:             `[text](url)`
#   - Raw-HTML passthrough: a line beginning with `<` is emitted as-is
#     and does not participate in paragraph collection.
#
# The source markdown files we ship under docs/report_content/ deliberately
# stay inside this subset. If a future edit needs a feature the parser
# doesn't cover, extend the parser rather than reach for a dep.

def md_to_html(md_text):
    """Convert a markdown string to an HTML fragment."""
    lines = md_text.splitlines()
    out = []
    para = []
    list_mode = None          # None | 'ul' | 'ol'
    in_code = False
    code_lang = ""
    code_buf = []

    def flush_para():
        if para:
            out.append(f"<p>{_md_inline(' '.join(para))}</p>")
            para.clear()

    def close_list():
        nonlocal list_mode
        if list_mode:
            out.append(f"</{list_mode}>")
            list_mode = None

    for raw in lines:
        line = raw.rstrip()
        if in_code:
            if line.startswith("```"):
                cls = f' class="language-{code_lang}"' if code_lang else ""
                out.append(
                    f"<pre><code{cls}>"
                    f"{html.escape(chr(10).join(code_buf))}</code></pre>")
                in_code, code_lang, code_buf = False, "", []
            else:
                code_buf.append(raw)
            continue

        if line.startswith("```"):
            flush_para(); close_list()
            in_code = True
            code_lang = line[3:].strip()
            continue

        # HTML comments <!-- ... --> are dropped so the per-file editor
        # hints don't leak into rendered output.
        if line.lstrip().startswith("<!--") and line.rstrip().endswith("-->"):
            continue

        if not line.strip():
            flush_para(); close_list()
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            flush_para(); close_list()
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{_md_inline(m.group(2))}</h{lvl}>")
            continue

        m = re.match(r"^\s*[-*]\s+(.*)$", line)
        if m:
            flush_para()
            if list_mode != "ul":
                close_list()
                out.append("<ul>")
                list_mode = "ul"
            out.append(f"<li>{_md_inline(m.group(1))}</li>")
            continue

        m = re.match(r"^\s*\d+\.\s+(.*)$", line)
        if m:
            flush_para()
            if list_mode != "ol":
                close_list()
                out.append("<ol>")
                list_mode = "ol"
            out.append(f"<li>{_md_inline(m.group(1))}</li>")
            continue

        # List-item continuation: an indented line (≥ 2 leading spaces
        # or a tab) directly under a list item is appended to that item
        # so multi-line bullets render as a single <li>.
        if (list_mode and out and out[-1].endswith("</li>")
                and (raw.startswith("  ") or raw.startswith("\t"))):
            cont = _md_inline(line.strip())
            out[-1] = out[-1][:-len("</li>")] + " " + cont + "</li>"
            continue

        # Blockquote: lines starting with `>`. Consecutive quote lines
        # join into one <blockquote>. (The original parser dropped the
        # marker into paragraph text; fix lifted with the list bug.)
        m = re.match(r"^>\s?(.*)$", line)
        if m:
            flush_para(); close_list()
            qtxt = _md_inline(m.group(1))
            if out and out[-1].startswith("<blockquote>") and out[-1].endswith("</blockquote>"):
                out[-1] = out[-1][:-len("</blockquote>")] + " " + qtxt + "</blockquote>"
            else:
                out.append(f"<blockquote>{qtxt}</blockquote>")
            continue

        # Raw-HTML block passthrough.
        if line.lstrip().startswith("<"):
            flush_para(); close_list()
            out.append(line)
            continue

        para.append(line)

    flush_para()
    close_list()
    if in_code:  # unterminated code fence — emit what we have
        cls = f' class="language-{code_lang}"' if code_lang else ""
        out.append(f"<pre><code{cls}>{html.escape(chr(10).join(code_buf))}</code></pre>")
    return "\n".join(out)


def _md_inline(text):
    """Apply inline markdown transforms to `text`. Raw HTML and entities
    pass through untouched; the content files in docs/report_content/
    are our own so we can trust them."""
    # 1. Protect `code` spans so their contents are not mangled.
    codes = []
    def _stash(m):
        codes.append(m.group(1))
        return f"\x00CODE{len(codes) - 1}\x00"
    text = re.sub(r"`([^`]+)`", _stash, text)
    # 2. Links: [text](url)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
                  lambda m: f'<a href="{m.group(2)}">{m.group(1)}</a>', text)
    # 3. Bold must be matched before italic because `**x**` would
    #    otherwise be read as two italic markers.
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    # 4. Italic: single * around non-whitespace content, not adjacent
    #    to another *.
    text = re.sub(r"(?<!\*)\*(?!\*|\s)([^*]+?)(?<!\s)\*(?!\*)",
                  r"<em>\1</em>", text)
    # 5. Restore code spans (escape their contents when emitting).
    def _restore(m):
        return f"<code>{html.escape(codes[int(m.group(1))])}</code>"
    text = re.sub(r"\x00CODE(\d+)\x00", _restore, text)
    return text


_SECTION_CACHE = {}


def load_section(name):
    """Load and convert a report_content markdown file. Cached by name.

    File location is resolved relative to this script so the function
    works both from the git tree and from the installed conda package,
    where `tc_rerender_report.py` sits at `$PREFIX/share/tidecluster/`
    and `docs/report_content/` is shipped alongside it."""
    if name in _SECTION_CACHE:
        return _SECTION_CACHE[name]
    root = Path(__file__).resolve().parent
    path = root / "docs" / "report_content" / f"{name}.md"
    if not path.exists():
        _SECTION_CACHE[name] = ""
        return ""
    _SECTION_CACHE[name] = md_to_html(path.read_text())
    return _SECTION_CACHE[name]


def copy_report_content_source(out_dir):
    """Copy the markdown source of each section into <out_dir>/docs/
    so users can download or edit the editable source alongside the
    rendered HTML."""
    root = Path(__file__).resolve().parent
    src = root / "docs" / "report_content"
    if not src.is_dir():
        return
    dst = Path(out_dir) / "docs"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


# ----------------------------------------------------------------------
# Small HTML helpers
# ----------------------------------------------------------------------

def esc(x):
    """Escape a cell value for safe HTML emission. None -> empty."""
    if x is None:
        return ""
    return html.escape(str(x), quote=True)


def fmt_bp(n):
    """Human-readable length in bp/kb/Mb."""
    if n is None: return ""
    try: n = int(n)
    except (TypeError, ValueError): return esc(n)
    if n >= 1_000_000: return f"{n/1_000_000:.2f} Mb"
    if n >= 1_000:     return f"{n/1_000:.1f} kb"
    return f"{n} bp"


def hor_badge(status):
    if not status:
        return ""
    cls = {"HOR strong":   "hor-strong",
           "HOR moderate": "hor-mod",
           "HOR weak":     "hor-weak",
           "No HOR":       "hor-none"}.get(status, "hor-none")
    return f'<span class="hor-badge {cls}">{esc(status)}</span>'


def hor_count_cell(n, kind):
    cls = {"strong": "hor-cnt-strong", "mod":  "hor-cnt-mod",
           "weak":   "hor-cnt-weak",   "none": "hor-cnt-none"}[kind]
    return f'<span class="{cls}">{int(n)}</span>'


# ----------------------------------------------------------------------
# kitehor combined_class -> user-facing label / CSS / structure summary
# ----------------------------------------------------------------------
# Canonical display order for kitehor >=0.10.0 (7 classes). Used by
# every roll-up table, the class-distribution bar, and the legends.
CANONICAL_CLASS_ORDER = [
    "hor",
    "hor_with_ssr",
    "tr",
    "tr_with_ssr",
    "tr_with_subrepeat",
    "pure_ssr",
    "unresolved",
]
COMBINED_CLASS_LABELS = {
    "hor":               "HOR",
    "hor_with_ssr":      "HOR (SSR-rich)",
    "tr":                "Simple TR",
    "tr_with_ssr":       "TR with SSR",
    "tr_with_subrepeat": "TR with subrepeats",
    "pure_ssr":          "SSR",
    "unresolved":        "Unresolved TR",
    # Legacy (kitehor 0.9.x) — kept so reports rerendered from old output
    # dirs still resolve a label/colour; not part of CANONICAL_CLASS_ORDER.
    "tr_with_nested_tr": "TR with nested TR",
}
COMBINED_CLASS_CSS = {
    "hor":               "hor",
    "hor_with_ssr":      "horssr",
    "tr":                "tr",
    "tr_with_ssr":       "trssr",
    "tr_with_subrepeat": "subrep",
    "pure_ssr":          "ssr",
    "unresolved":        "unres",
    "tr_with_nested_tr": "nested",   # legacy 0.9.x
}


def class_label(raw):
    """Friendly label for a kitehor combined_class value."""
    if not raw:
        return "—"
    return COMBINED_CLASS_LABELS.get(raw, raw.replace("_", " "))


def class_badge(raw):
    """Colored pill for a combined_class value (empty string when unknown)."""
    if not raw:
        return ""
    css = COMBINED_CLASS_CSS.get(raw, "unres")
    return f'<span class="class-badge class-{css}">{esc(class_label(raw))}</span>'


def class_count_badge(raw, n):
    """Small count pill colored by class (for roll-up "class mix" cells)."""
    css = COMBINED_CLASS_CSS.get(raw, "unres")
    return (f'<span class="class-cnt class-{css}" '
            f'title="{esc(class_label(raw))}">{int(n)}</span>')


def class_mix_cell(counts):
    """Per-class count pills in canonical display order; '—' when empty.
    Any class outside the canonical set (e.g. legacy tr_with_nested_tr on
    a rerendered 0.9.x dir) is appended after the canonical ones."""
    if not counts:
        return "—"
    bits = [class_count_badge(raw, counts[raw])
            for raw in CANONICAL_CLASS_ORDER if counts.get(raw)]
    bits += [class_count_badge(raw, n) for raw, n in counts.items()
             if raw not in CANONICAL_CLASS_ORDER and n]
    return " ".join(bits) or "—"


def class_totals(model):
    """Counter of combined_class -> total arrays across all TRCs."""
    total = Counter()
    for t in model["trcs"]:
        k = t.get("kite") or {}
        for raw, n in (k.get("combined_class_counts") or {}).items():
            total[raw] += n
    return total


def class_legend_inline():
    """One-line legend of the class badges, canonical display order."""
    return " ".join(class_badge(raw) for raw in CANONICAL_CLASS_ORDER)


def _sc(s):
    """Compact score for the monomer-estimates cell: 0.31 -> '.31'."""
    try:
        txt = f"{float(s):.2f}"
    except (TypeError, ValueError):
        return ""
    return txt[1:] if txt.startswith("0.") else txt


def monomer_estimates_cell(a):
    """Top-5 KITE peaks as 'period(score)', ranked, for the array table."""
    bits = []
    for mk, sk in (("m1", "s1"), ("m2", "s2"), ("m3", "s3"),
                   ("m4", "s4"), ("m5", "s5")):
        m = a.get(mk)
        if m is None:
            continue
        s = a.get(sk)
        sc = _sc(s)
        bits.append(f'{esc(m)}<span class="tc-score">({sc})</span>' if sc
                    else f'{esc(m)}')
    return " ".join(bits)


def _fmt_motifs(top_motifs):
    """`AG:45.3%;AGAT:2.8%` -> `AG 45.3%, AGAT 2.8%` (escaped)."""
    if not top_motifs:
        return ""
    parts = []
    for item in re.split(r"[;,]", top_motifs):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            mot, pct = item.split(":", 1)
            parts.append(f"{esc(mot.strip())} {esc(pct.strip())}")
        else:
            parts.append(esc(item))
    return ", ".join(parts)


def _bp(x):
    """Render a bp value as an integer (kitehor emits some as floats)."""
    if x is None:
        return None
    try:
        return int(round(float(x)))
    except (TypeError, ValueError):
        return None


def structure_cell(a):
    """Adaptive 'Structure' cell: render only the fields meaningful for
    this array's combined_class.

    kitehor >=0.10.0: HOR base/tile/multiplicity from hor_*; subrepeat
    host/period from the tandem_validate tv_* columns; SSR coverage +
    motif list from ssr_*. The tr_with_nested_tr branch is legacy
    (0.9.x) and uses the old subrepeat_* columns."""
    cc = a.get("combined_class")
    base = a.get("hor_base_monomer")
    tile = a.get("hor_hor_period")
    mult = a.get("hor_multiplicity")
    ssr_tot = a.get("ssr_total_coverage_pct")
    motifs = _fmt_motifs(a.get("ssr_top_motifs"))
    cons = a.get("consensus_period_bp")
    tv_host = _bp(a.get("tv_host_period"))
    tv_sub = _bp(a.get("tv_best_candidate_period"))

    def _ssr_seg(decimals=0):
        if ssr_tot is None:
            return None
        s = f"SSR {ssr_tot:.{decimals}f}%"
        if motifs:
            s += f": {motifs}"
        return s

    def _hor_segs():
        bits = []
        if base: bits.append(f"base {esc(base)}")
        if tile: bits.append(f"HOR {esc(tile)}")
        if mult: bits.append(f"&times;{esc(mult)}")
        return bits

    if cc == "hor":
        return " &middot; ".join(_hor_segs())
    if cc == "hor_with_ssr":
        bits = _hor_segs()
        seg = _ssr_seg()
        if seg: bits.append(seg)
        return " &middot; ".join(bits)
    if cc == "tr":
        return f"base monomer {esc(base)} bp" if base else ""
    if cc == "tr_with_ssr":
        bits = []
        if base: bits.append(f"base {esc(base)}")
        seg = _ssr_seg()
        if seg: bits.append(seg)
        return " &middot; ".join(bits)
    if cc == "tr_with_subrepeat":
        # tandem_validate: tv_host_period is the host monomer,
        # tv_best_candidate_period the localized subrepeat.
        mono = tv_host or base or tile
        bits = []
        if mono: bits.append(f"monomer {esc(mono)}")
        if tv_sub: bits.append(f"subrepeat {esc(tv_sub)}")
        return " &middot; ".join(bits)
    if cc == "pure_ssr":
        return _ssr_seg(decimals=1) or motifs
    if cc == "unresolved":
        return f"period ~{esc(cons)} bp (unresolved)" if cons else "unresolved"
    # Legacy 0.9.x nested-TR (subrepeat_* columns; None on 0.10.0 data).
    if cc == "tr_with_nested_tr":
        sub = a.get("subrepeat_period_bp")
        subcov = a.get("subrepeat_coverage_pct")
        bits = []
        if base: bits.append(f"base {esc(base)}")
        if sub:
            cov = f" ({subcov:.0f}%)" if subcov is not None else ""
            bits.append(f"sub-TR {esc(sub)}{cov}")
        return " &middot; ".join(bits)
    # No combined_class (legacy/below-threshold): base monomer if any.
    return f"base monomer {esc(base)} bp" if base else ""


def page_shell(title, active, run_meta, main_html, *,
               extra_head="", up="", root_href="index.html",
               legacy_href=None, assets_href="assets/"):
    """Return a full HTML document string wrapping `main_html`.

    `up`          — prefix added to sibling-page nav links
                    (tarean / kite / superfamilies). "" for top-level
                    pages in <prefix>_report/, "../" for per-TRC
                    dashboards, "<prefix>_report/" for the root index.
    `root_href`   — full href for the Summary nav entry; this is the
                    root <prefix>_index.html and its depth differs at
                    each page level.
    `legacy_href` — full href for the legacy-report nav link; computed
                    by the caller based on page depth.
    `assets_href` — path to the assets/ directory."""
    nav_items = [
        (root_href,                "Summary",       "index"),
        (f"{up}tarean.html",       "TAREAN",        "tarean"),
        (f"{up}kite.html",         "KITE",          "kite"),
        (f"{up}superfamilies.html","Superfamilies", "sf"),
    ]
    def _nav_link(href, label, key):
        cls = ' class="active"' if key == active else ""
        return f'<a href="{href}"{cls}>{label}</a>'
    nav_html = "".join(_nav_link(*item) for item in nav_items)
    legacy_link = ""
    if legacy_href:
        legacy_link = (f'<a href="{legacy_href}" class="tc-nav-ext" '
                       f'title="Open the original v1 TideCluster report">'
                       f'Legacy report &#x2197;</a>')
    run_strip = (
        f'<strong>{esc(run_meta["prefix"])}</strong>'
        f' · TideCluster {esc(run_meta["version"])}'
        f' · generated {esc(run_meta["generated_at"])}'
        f' · {run_meta["stats_line"]}'
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)} — {esc(run_meta["prefix"])}</title>
<link rel="stylesheet" href="{assets_href}datatables/datatables.min.css">
<link rel="stylesheet" href="{assets_href}tidecluster.css">
{extra_head}
</head>
<body class="tc-page-{esc(active)}">
<nav class="tc-nav">
  <span class="tc-brand">TideCluster</span>
  {nav_html}
  <span class="tc-spacer"></span>
  {legacy_link}
</nav>
<div class="tc-run-strip">{run_strip}</div>
<main>
{main_html}
</main>
<script src="{assets_href}datatables/jquery-3.7.1.min.js"></script>
<script src="{assets_href}datatables/datatables.min.js"></script>
<script src="{assets_href}tidecluster.js"></script>
</body>
</html>
"""


def make_run_meta(model):
    stats = model.get("stats", {})
    parts = []
    if "n_trcs_total" in stats:           parts.append(f'{stats["n_trcs_total"]} TRCs')
    if "n_trcs_above_threshold" in stats: parts.append(f'{stats["n_trcs_above_threshold"]} analysed')
    if "n_tras" in stats:                 parts.append(f'{stats["n_tras"]} arrays')
    if "input_sequence_length" in stats:  parts.append(f'{fmt_bp(stats["input_sequence_length"])} input')
    return {
        "prefix":       model["meta"]["prefix"],
        "version":      stats.get("tidecluster_version", "?"),
        "generated_at": model["meta"]["generated_at"],
        "stats_line":   " · ".join(parts),
    }

# ----------------------------------------------------------------------
# GFF3 parsing (minimal, no dependency on tc_utils)
# ----------------------------------------------------------------------

def parse_gff3_attrs(attrs_str):
    """Parse a GFF3 attribute string into a dict."""
    out = {}
    if not attrs_str or attrs_str == ".":
        return out
    for kv in attrs_str.split(";"):
        kv = kv.strip()
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        out[k] = v
    return out


def read_gff3(path):
    """Yield dicts with {seqid,source,type,start,end,score,strand,phase,attrs}."""
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            yield {
                "seqid":  parts[0],
                "source": parts[1],
                "type":   parts[2],
                "start":  int(parts[3]),
                "end":    int(parts[4]),
                "score":  parts[5],
                "strand": parts[6],
                "phase":  parts[7],
                "attrs":  parse_gff3_attrs(parts[8]),
            }


# ----------------------------------------------------------------------
# Source file discovery
# ----------------------------------------------------------------------

def detect_prefix(input_dir: Path) -> str:
    """Derive prefix from filenames like <prefix>_clustering.gff3."""
    matches = sorted(p.name[:-len("_clustering.gff3")]
                     for p in input_dir.glob("*_clustering.gff3"))
    if not matches:
        raise SystemExit(
            f"no *_clustering.gff3 found in {input_dir}; pass --prefix")
    if len(matches) > 1:
        raise SystemExit(
            f"multiple prefixes in {input_dir} ({matches}); pass --prefix")
    return matches[0]


def resolve_paths(input_dir: Path, prefix: str) -> dict:
    """Map every source file to its path (or None if absent)."""
    def first(*names):
        for n in names:
            p = input_dir / n
            if p.exists():
                return p
        return None
    kite_dir   = input_dir / f"{prefix}_kite"
    tarean_dir = input_dir / f"{prefix}_tarean"
    return {
        "cmd_args":              first(f"{prefix}_cmd_args.json"),
        "pipeline_stats":        first(f"{prefix}_pipeline_stats.json"),
        "seqid_lengths":         first(f"{prefix}_seqid_lengths.tsv"),
        "clustering_gff3":       first(f"{prefix}_clustering.gff3"),
        "annotation_gff3":       first(f"{prefix}_annotation.gff3"),
        "tarean_report_tsv":     first(f"{prefix}_tarean_report.tsv"),
        "tarean_report_html":    first(f"{prefix}_tarean_report.html"),
        "kite_top3_csv":         kite_dir / "monomer_size_top3_estimats.csv" if (kite_dir / "monomer_size_top3_estimats.csv").exists() else None,
        "kite_best_csv":         kite_dir / "monomer_size_best_estimate_stat.csv" if (kite_dir / "monomer_size_best_estimate_stat.csv").exists() else None,
        # kitehor >=0.12.0 rescored peaks TSV (24 cols, per-peak with the
        # 15 rescore diagnostics). Read for the per-array Details panel.
        "rescored_peaks_tsv":    kite_dir / "kitehor.rescored.peaks.tsv" if (kite_dir / "kitehor.rescored.peaks.tsv").exists() else None,
        "kite_dir":              kite_dir if kite_dir.exists() else None,
        "tarean_dir":            tarean_dir if tarean_dir.exists() else None,
        "superfamilies_csv":     first(f"{prefix}_trc_superfamilies.csv"),
        "superfamilies_html":    first(f"{prefix}_trc_superfamilies.html"),
        "dotplots_dir":          input_dir / "dotplots" if (input_dir / "dotplots").exists() else None,
        "index_html":            first(f"{prefix}_index.html"),
        "consensus_library":     first(f"{prefix}_consensus_dimer_library.fasta"),
    }


# ----------------------------------------------------------------------
# Fallback: scrape pipeline_stats from index.html (<pre> summary block)
# ----------------------------------------------------------------------

_STATS_KEY_MAP = {
    "Number of TRCs":                  "n_trcs_total",
    "Number of TRCs above threshold":  "n_trcs_above_threshold",
    "Number of SSRs in TRCs":          "n_ssrs",
    "Total length of TRCs":            "total_tr_length",
    "Number of TRAs":                  "n_tras",
    "Input sequence length":           "input_sequence_length",
}


def parse_stats_from_index_html(path: Path) -> dict:
    """Extract key:value numbers from the Summary <pre> block."""
    stats = {}
    txt = path.read_text(errors="replace")
    for pretty, key in _STATS_KEY_MAP.items():
        m = re.search(re.escape(pretty) + r"\s*:\s*([0-9]+)", txt)
        if m:
            stats[key] = int(m.group(1))
    # tidecluster version
    m = re.search(r"TideCluster version\s*:\s*([0-9A-Za-z.\-]+)", txt)
    if m:
        stats["tidecluster_version"] = m.group(1)
    return stats


def load_stats(paths: dict) -> dict:
    if paths["pipeline_stats"]:
        with open(paths["pipeline_stats"]) as f:
            return json.load(f)
    if paths["index_html"]:
        return parse_stats_from_index_html(paths["index_html"])
    return {}


def load_seqid_lengths(paths: dict) -> dict:
    """Return seqid -> length (bp). Prefers <prefix>_seqid_lengths.tsv
    emitted by TideCluster.py tarean(); falls back to a lower-bound
    derived from max(end) across the tidehunter + clustering GFF3s for
    legacy outputs that pre-date the sidecar. A missing seqid is
    fine; consumers use .get() with a per-array-coord fallback."""
    p = paths.get("seqid_lengths")
    out = {}
    if p and p.exists():
        with open(p) as f:
            next(f, None)  # header
            for line in f:
                s, l = line.rstrip("\n").split("\t")[:2]
                try: out[s] = int(l)
                except ValueError: continue
        return out
    # Legacy fallback: scan every GFF3 for max(end) per seqid.
    for key in ("clustering_gff3", "annotation_gff3"):
        p = paths.get(key)
        if not p: continue
        for feat in read_gff3(p):
            sid = feat["seqid"]
            if feat["end"] > out.get(sid, 0):
                out[sid] = feat["end"]
    th = paths.get("input_dir") if isinstance(paths.get("input_dir"), Path) else None
    # Also mine the raw tidehunter GFF3 if available (max(end) there is
    # a tighter lower bound since every detected array contributes).
    prefix = None
    ann = paths.get("clustering_gff3")
    if ann:
        name = ann.name
        if name.endswith("_clustering.gff3"):
            prefix = name[:-len("_clustering.gff3")]
            th_path = ann.parent / f"{prefix}_tidehunter.gff3"
            if th_path.exists():
                for feat in read_gff3(th_path):
                    sid = feat["seqid"]
                    if feat["end"] > out.get(sid, 0):
                        out[sid] = feat["end"]
    return out


# ----------------------------------------------------------------------
# Per-TRA rows from clustering.gff3 / annotation.gff3
# ----------------------------------------------------------------------

def load_tras(paths: dict):
    """Return per-TRA dict list grouped by TRC_ID.

    Uses annotation.gff3 when available (richer attrs), else clustering."""
    gff_path = paths["annotation_gff3"] or paths["clustering_gff3"]
    if not gff_path:
        return {}
    trcs = {}
    for feat in read_gff3(gff_path):
        attrs = feat["attrs"]
        trc_id = attrs.get("Name")
        if not trc_id:
            continue
        entry = {
            "seqid":  feat["seqid"],
            "start":  feat["start"],
            "end":    feat["end"],
            "length": feat["end"] - feat["start"],
            "strand": feat["strand"],
            "ssr":          attrs.get("ssr"),
            "annotation":   attrs.get("annotation"),
            "repeat_type":  attrs.get("repeat_type"),
            "repeat_unit_seq": attrs.get("repeat_unit_seq"),
            "copy_number":  attrs.get("copy_number"),
        }
        trcs.setdefault(trc_id, []).append(entry)
    return trcs


# ----------------------------------------------------------------------
# TSV / CSV loaders (quoted CSV with embedded newlines, etc.)
# ----------------------------------------------------------------------

def _read_rows(path: Path, delim: str):
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter=delim):
            yield {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}


def load_tarean_summary(paths):
    """Per-TRC TAREAN summary from <prefix>_tarean_report.tsv.

    Returns dict keyed by TRC_ID. Some cells contain HTML (img tags,
    <pre>consensus<pre>); we keep the raw strings and extract the PNG
    path for the graph and logo as tidy fields."""
    p = paths["tarean_report_tsv"]
    if not p:
        return {}
    out = {}
    img_src = re.compile(r'src=\\?"([^"\\]+)\\?"')
    for row in _read_rows(p, "\t"):
        trc = (row.get("TRC") or "").strip()
        if not trc:
            continue
        def pull_img(field):
            raw = row.get(field) or ""
            m = img_src.search(raw)
            return m.group(1) if m else None
        def maybe_num(v, to=float):
            try: return to(v)
            except (TypeError, ValueError): return None
        # consensus is wrapped in <pre>...<pre> (note: broken tag in
        # original HTML). Strip the tags, keep just the sequence.
        raw_cons = (row.get("Consensus") or "").strip()
        cons_seq = re.sub(r"<\/?pre>", "", raw_cons).replace("\n", "")
        out[trc] = {
            "index":          maybe_num(row.get("INDEX"), int),
            "kmer":           maybe_num(row.get("kmer"), int),
            "variant":        None,  # not in TSV; pulled from per-TRC summary_table if needed
            "monomer_length": maybe_num(row.get("monomer_length"), int),
            "total_score":    maybe_num(row.get("total_score")),
            "n_gap50":        maybe_num(row.get("n_gap50"), int),
            "tarean_dir":     (row.get("tarean_dir") or "").strip(),
            "graph_png":      pull_img("graph_link"),
            "logo_png":       pull_img("logo_link"),
            "consensus":      cons_seq,
            "annotation":     (row.get("Annotation") or "").strip() or None,
            "ssr_motifs":     (row.get("SSRs") or "").strip() or None,
            "total_size":     maybe_num(row.get("Total_size"), int),
            "n_arrays":       maybe_num(row.get("number_of_arrays"), int),
            "min_array":      maybe_num(row.get("min_array_length"), int),
            "max_array":      maybe_num(row.get("max_array_length"), int),
            "median_array":   maybe_num(row.get("median_array_length")),
            "type":           (row.get("type") or "").strip() or None,
        }
    return out


def load_kite_top3(paths):
    """Per-TRA KITE records keyed by (trc_id, seqid, start, end).

    Four CSV schemas are accepted:
      - 1.11+ (kitehor >=0.12.0): adds `founder_period`,
        `strongest_period`, `multiplicity`, `delta_id_pp`,
        `founder_id_med`, `strongest_id_med`, `founder_fallback` +
        top-2 subrepeat candidates (HIGH / LIKELY tier). The
        combined_class / tv_* columns are gone (the cascade is
        dropped from the pipeline).
      - 1.10  (kitehor 0.10.0): lowercase `hor_*` + `combined_class`
        + `tv_*` (tandem_validate). Surfaces in legacy mode.
      - 1.9.x (legacy R kite.R): capitalised `HOR_status` /
        `HOR_confidence` / `HOR_base_monomer` / `HOR_hor_period` /
        `HOR_n_harmonics`.
      - pre-1.9.0: no HOR columns; recomputed via `compute_hor()`."""
    p = paths["kite_top3_csv"]
    if not p:
        return {}
    out = {}
    # Peek at headers to decide schema.
    with open(p, newline="") as f:
        header = next(csv.reader(f, delimiter="\t"))
    is_v012   = "founder_period" in header
    is_kitehor = (not is_v012) and "hor_status" in header
    is_legacy  = "HOR_status" in header
    for row in _read_rows(p, "\t"):
        trc = (row.get("TRC_ID") or "").strip()
        sid = (row.get("seqid") or "").strip()
        try:
            start = int(row.get("start")); end = int(row.get("end"))
        except (TypeError, ValueError):
            continue
        def num(v):
            try:
                if v in (None, "", "NA"): return None
                return float(v)
            except (TypeError, ValueError): return None
        def ni(v):
            x = num(v)
            return int(x) if x is not None else None
        m1 = ni(row.get("monomer_size"));   s1 = num(row.get("score"))
        m2 = ni(row.get("monomer_size_2")); s2 = num(row.get("score_2"))
        m3 = ni(row.get("monomer_size_3")); s3 = num(row.get("score_3"))
        m4 = ni(row.get("monomer_size_4")); s4 = num(row.get("score_4"))
        m5 = ni(row.get("monomer_size_5")); s5 = num(row.get("score_5"))

        if is_v012:
            hor = {
                # The cascade is gone; keep these keys present (None) so
                # downstream model code that still references them doesn't
                # KeyError. The new fields below carry the meaning.
                "hor_status": None, "hor_confidence": None,
                "hor_base_monomer": None, "hor_hor_period": None,
                "hor_n_harmonics": 0, "hor_multiplicity": None,
                "combined_class": None,
                # v0.12 per-array summary
                "founder_period":   ni(row.get("founder_period")),
                "strongest_period": ni(row.get("strongest_period")),
                "multiplicity":     ni(row.get("multiplicity")) or 1,
                "multiplicity_raw": num(row.get("multiplicity_raw")),
                "irregular_multiplicity":
                                    (row.get("irregular_multiplicity") or "").strip() == "true",
                "delta_id_pp":      num(row.get("delta_id_pp")),
                "founder_id_med":   num(row.get("founder_id_med")),
                "strongest_id_med": num(row.get("strongest_id_med")),
                "founder_fallback": (row.get("founder_fallback") or "").strip() == "true",
                # Top-2 subrepeat candidates (main column)
                "subrepeat_1_period": ni(row.get("subrepeat_1_period")),
                "subrepeat_1_occ":    num(row.get("subrepeat_1_occ")),
                "subrepeat_1_tier":   (row.get("subrepeat_1_tier") or "").strip() or None,
                "subrepeat_2_period": ni(row.get("subrepeat_2_period")),
                "subrepeat_2_occ":    num(row.get("subrepeat_2_occ")),
                "subrepeat_2_tier":   (row.get("subrepeat_2_tier") or "").strip() or None,
                # SSR
                "ssr_flag":           (row.get("ssr_flag") or "").strip() or None,
                "ssr_dominant_motif": (row.get("ssr_dominant_motif") or "").strip() or None,
                "ssr_dominant_motif_coverage_pct": num(row.get("ssr_dominant_motif_coverage_pct")),
                "ssr_total_coverage_pct": num(row.get("ssr_total_coverage_pct")),
                "ssr_top_motifs":     (row.get("ssr_top_motifs") or "").strip() or None,
                "consensus_period_bp": ni(row.get("consensus_period_bp")),
            }
        elif is_kitehor:
            hor = {
                "hor_status":       (row.get("hor_status") or "").strip() or "No HOR",
                "hor_confidence":   num(row.get("hor_confidence")) or 0.0,
                # Map kitehor's hor_founder / hor_tile / hor_multiplicity to
                # the names the downstream model already publishes.
                "hor_base_monomer": ni(row.get("hor_founder")),
                "hor_hor_period":   ni(row.get("hor_tile")),
                "hor_n_harmonics":  ni(row.get("hor_multiplicity")) or 0,
                "hor_multiplicity": ni(row.get("hor_multiplicity")),
                # Structural fields surfaced by kitehor. The tv_* columns
                # (kitehor >=0.10.0 tandem_validate) hold the subrepeat
                # signal: tv_host_period = host monomer, tv_best_candidate_
                # period = subrepeat. The subrepeat_* keys are kept so
                # reports rerendered from 0.9.x output dirs still resolve
                # (they are None on 0.10.0 data).
                "combined_class":   (row.get("combined_class") or "").strip() or None,
                "tv_decision":      (row.get("tv_decision") or "").strip() or None,
                "tv_host_period":   num(row.get("tv_host_period")),
                "tv_best_candidate_period": num(row.get("tv_best_candidate_period")),
                "tv_best_candidate_kind": (row.get("tv_best_candidate_kind") or "").strip() or None,
                "tv_density":       num(row.get("tv_density")),
                "tv_n_windows_total":   ni(row.get("tv_n_windows_total")),
                "tv_n_windows_present": ni(row.get("tv_n_windows_present")),
                "ssr_flag":         (row.get("ssr_flag") or "").strip() or None,
                "ssr_dominant_motif": (row.get("ssr_dominant_motif") or "").strip() or None,
                "ssr_total_coverage_pct": num(row.get("ssr_total_coverage_pct")),
                "ssr_top_motifs":   (row.get("ssr_top_motifs") or "").strip() or None,
                "consensus_period_bp": ni(row.get("consensus_period_bp")),
                # Legacy 0.9.x subrepeat columns (None on 0.10.0 data).
                "subrepeat_period_bp": ni(row.get("subrepeat_period_bp")),
                "subrepeat_coverage_pct": num(row.get("subrepeat_coverage_pct")),
            }
        elif is_legacy:
            hor = {
                "hor_status":       (row.get("HOR_status") or "").strip() or "No HOR",
                "hor_confidence":   num(row.get("HOR_confidence")) or 0.0,
                "hor_base_monomer": ni(row.get("HOR_base_monomer")),
                "hor_hor_period":   ni(row.get("HOR_hor_period")),
                "hor_n_harmonics":  ni(row.get("HOR_n_harmonics")) or 0,
            }
        else:
            hor = compute_hor(m1, m2, m3, s1, s2, s3)

        out[(trc, sid, start, end)] = {
            "m1": m1, "s1": s1, "m2": m2, "s2": s2, "m3": m3, "s3": s3,
            "m4": m4, "s4": s4, "m5": m5, "s5": s5,
            "array_length": ni(row.get("array_length")),
            **hor,
        }
    return out


def load_rescored_peaks(paths):
    """Per-array rescored peaks from kitehor.rescored.peaks.tsv, keyed
    by (trc_id, seqid, start, end). Each value is a list of dicts in
    rank order (one dict per kite peak, with the 24 columns kitehor
    >=0.12.0 emits). Drives the per-array Details child row."""
    p = paths.get("rescored_peaks_tsv")
    if not p:
        return {}
    out = {}
    for row in _read_rows(p, "\t"):
        rec = row.get("case_id") or ""
        if ":" not in rec:
            continue
        trc, rest = rec.split(":", 1)
        parts = rest.rsplit("_", 2)
        if len(parts) != 3:
            continue
        sid, s, e = parts
        try:
            key = (trc, sid, int(s), int(e))
        except ValueError:
            continue
        out.setdefault(key, []).append(row)
    # Sort each per-array list by rank
    for v in out.values():
        v.sort(key=lambda r: int(r.get("rank") or 0))
    return out


# ----------------------------------------------------------------------
# Subrepeat tier classifier — mirrors tc_utils._classify_subrepeat_tier
# so the per-array Details panel can label every non-rejected peak.
# ----------------------------------------------------------------------
def classify_peak_tier(peak, founder_period):
    def _n(v):
        if v in (None, "", "NA"): return None
        try: return float(v)
        except (TypeError, ValueError): return None
    if peak.get("phantom", "").strip() == "true":
        return "REJECT_PHANTOM"
    period = _n(peak.get("period"))
    if period is None or period <= 0:
        return "REJECT_BAD"
    occ = _n(peak.get("scan_occupancy_frac"))
    if founder_period is None or founder_period <= 0:
        return "OBSERVATIONAL" if (occ is not None and occ >= 0.05) else "REJECT_NO_FOUNDER"
    ratio = period / float(founder_period)
    idm    = _n(peak.get("identity_med"))
    covf   = _n(peak.get("coverage_frac"))
    scan_n = _n(peak.get("scan_n_intervals"))
    phaseC = _n(peak.get("kmer_phase_contrast"))
    autoF  = _n(peak.get("kmer_autocorr_founder"))
    srflag = peak.get("subrepeat", "").strip() == "true"
    if (idm is not None and idm >= 0.85
        and covf is not None and covf >= 0.95
        and occ  is not None and occ  >= 0.95):
        return "REJECT_IS_FOUNDER"
    if ratio > 0.33:
        return "REJECT_TOO_LARGE"
    if ratio <= 0.25:
        if (occ is not None and occ >= 0.15
            and (srflag
                 or (phaseC is not None and phaseC >= 0.10)
                 or (autoF  is not None and autoF  >= 0.4))):
            return "HIGH"
        if (occ is not None and occ >= 0.20
            and scan_n is not None and scan_n >= 10):
            return "LIKELY"
        if ((phaseC is not None and phaseC >= 0.10)
            or (autoF is not None and autoF  >= 0.4)):
            return "KMER_SUPPORT"
        return "WEAK"
    if ((occ is not None and occ >= 0.20)
        or (phaseC is not None and phaseC >= 0.10)
        or (autoF  is not None and autoF  >= 0.4)):
        return "AMBIGUOUS"
    return "REJECT_AMBIGUOUS_LOW"


def load_per_trc_variants(paths, trc_id):
    """TAREAN variants table for a single TRC (summary_table.csv)."""
    if not paths["tarean_dir"]:
        return []
    p = paths["tarean_dir"] / f"{trc_id}.fasta_tarean" / "summary_table.csv"
    if not p.exists():
        return []
    rows = []
    for row in _read_rows(p, ","):
        rows.append({
            "kmer":             int(row["kmer"])     if row.get("kmer") else None,
            "variant":          int(row["variant"])  if row.get("variant") else None,
            "total_score":      float(row["total_score"])   if row.get("total_score") else None,
            "monomer_length":   int(row["monomer_length"])  if row.get("monomer_length") else None,
            "consensus_length": int(row["consensus_length"])if row.get("consensus_length") else None,
            "n_gap50":          int(row["n_gap50"])         if row.get("n_gap50") else None,
            "graph_link":       (row.get("graph_link") or "").strip() or None,
            "logo_link":        (row.get("logo_link") or "").strip() or None,
        })
    return rows


def load_superfamilies(paths, input_dir):
    """TRC-to-superfamily mapping + dotplot path (None if no SF analysis).

    Dotplot path is emitted RELATIVE to input_dir so the HTML side only
    needs to prepend "../" to resolve it from report_v2/."""
    p = paths["superfamilies_csv"]
    if not p:
        return []
    # CSV format: "Superfamily","TRC"
    groups = {}
    for row in _read_rows(p, ","):
        sf = row.get("Superfamily")
        trc = row.get("TRC")
        if not sf or not trc: continue
        try: sf_id = int(sf)
        except ValueError: continue
        groups.setdefault(sf_id, []).append(trc.strip())
    out = []
    dotdir = paths["dotplots_dir"]
    for sf_id in sorted(groups):
        members = sorted(groups[sf_id], key=_trc_sort_key)
        dotplot = None
        if dotdir and len(members) > 1:
            indices = sorted(int(x.replace("TRC_", "")) for x in members)
            fname = "TRCS_" + "_".join(str(i) for i in indices) + ".png"
            candidate = dotdir / fname
            if candidate.exists():
                dotplot = str(candidate.relative_to(input_dir))
        out.append({"id": sf_id, "trcs": members, "dotplot": dotplot})
    return out


def _trc_sort_key(trc_id):
    try:
        return int(trc_id.split("_", 1)[1])
    except (ValueError, IndexError):
        return 10**9


# ----------------------------------------------------------------------
# Model assembly
# ----------------------------------------------------------------------

def build_model(input_dir: Path, prefix: str):
    paths = resolve_paths(input_dir, prefix)
    settings = json.loads(paths["cmd_args"].read_text()) if paths["cmd_args"] else {}
    stats = load_stats(paths)
    seqid_lengths    = load_seqid_lengths(paths)
    tras_by_trc      = load_tras(paths)
    tarean_summary   = load_tarean_summary(paths)
    kite_by_array    = load_kite_top3(paths)
    rescored_peaks   = load_rescored_peaks(paths)
    superfams        = load_superfamilies(paths, input_dir)
    sf_by_trc = {trc: sf["id"] for sf in superfams for trc in sf["trcs"]}
    sf_peers  = {trc: [t for t in sf["trcs"] if t != trc]
                 for sf in superfams for trc in sf["trcs"] if len(sf["trcs"]) > 1}
    sf_dotplot_by_trc = {trc: sf["dotplot"]
                         for sf in superfams for trc in sf["trcs"] if sf["dotplot"]}

    kite_dir_rel = paths["kite_dir"].name if paths["kite_dir"] else None

    trcs = []
    for trc_id in sorted(tras_by_trc, key=_trc_sort_key):
        arrays = tras_by_trc[trc_id]
        # Merge KITE per-array fields
        kite_counts = {"no_hor": 0, "hor_weak": 0, "hor_moderate": 0,
                       "hor_strong": 0}
        confidences = []
        m1_values = []
        combined_class_counter = Counter()
        n_subrepeat = 0
        n_ssr = 0
        for arr in arrays:
            key = (trc_id, arr["seqid"], arr["start"], arr["end"])
            k = kite_by_array.get(key)
            if k:
                arr.update(k)
                # Attach the per-array rescored peak list (kitehor v0.12+)
                # for the Details child row. List of dicts (rank-sorted).
                arr["rescored_peaks"] = rescored_peaks.get(key) or []
                status = k.get("hor_status")
                if   status == "HOR strong":    kite_counts["hor_strong"]   += 1
                elif status == "HOR moderate":  kite_counts["hor_moderate"] += 1
                elif status == "HOR weak":      kite_counts["hor_weak"]     += 1
                elif status == "No HOR":        kite_counts["no_hor"]       += 1
                c = k.get("hor_confidence")
                if c is not None: confidences.append(c)
                if k.get("m1") is not None:
                    m1_values.append(k["m1"])
                cc = k.get("combined_class")
                if cc: combined_class_counter[cc] += 1
                # 0.10.0 dropped the standalone subrepeat_flag column;
                # derive the structural counters from combined_class
                # (works for both 0.9.x and 0.10.0 data).
                if cc in {"tr_with_subrepeat", "tr_with_nested_tr"}:
                    n_subrepeat += 1
                if cc in {"pure_ssr", "tr_with_ssr", "hor_with_ssr"} or \
                        (k.get("ssr_flag") or "").lower() in {"yes", "y", "true"}:
                    n_ssr += 1
        # TRC-level aggregates
        repeat_type = next((a["repeat_type"] for a in arrays if a.get("repeat_type")), None)
        annotation  = next((a["annotation"]  for a in arrays if a.get("annotation")),  None)
        ssr_motif   = next((a["ssr"]         for a in arrays if a.get("ssr")),         None)
        lengths = [a["length"] for a in arrays]
        tarean = tarean_summary.get(trc_id)
        has_tarean = tarean is not None and (paths["tarean_dir"]
                                             and (paths["tarean_dir"] / f"{trc_id}.fasta_tarean").exists())
        has_kite = bool(m1_values)
        kite_block = None
        if has_kite and kite_dir_rel:
            # "monomer_primary" = most-frequent m1 across the TRC's arrays
            # (matches original kite_report.R's monomer_size column).
            monomer_primary = Counter(m1_values).most_common(1)[0][0]
            dominant_class = (combined_class_counter.most_common(1)[0][0]
                              if combined_class_counter else None)
            # TRC-level aggregates (v0.12). prevalent_founder = median of
            # per-array founder_period values that cluster within ±10 % of
            # each other (the same anchor used by build_monomer_size_csv's
            # Pass 2 TRC consensus). multiplicity_counts = distribution of
            # the rounded multiplicity for HOR arrays (×k≥2).
            founder_periods = [a["founder_period"] for a in arrays
                               if a.get("founder_period")]
            prev_founder = None
            prev_founder_n = 0
            if len(founder_periods) >= 3:
                best_n, best_members = 0, []
                for anchor in founder_periods:
                    members = [p for p in founder_periods
                               if abs(p - anchor) / max(anchor, 1.0) <= 0.10]
                    if len(members) > best_n:
                        best_n, best_members = len(members), members
                if best_n >= 3:
                    prev_founder   = int(round(statistics.median(best_members)))
                    prev_founder_n = best_n
            mult_counts = Counter(
                int(a["multiplicity"]) for a in arrays
                if (a.get("multiplicity") or 1) > 1)
            n_irreg = sum(1 for a in arrays if a.get("irregular_multiplicity"))
            kite_block = {
                "monomer_primary":  int(monomer_primary),
                "n_no_hor":         kite_counts["no_hor"],
                "n_hor_weak":       kite_counts["hor_weak"],
                "n_hor_moderate":   kite_counts["hor_moderate"],
                "n_hor_strong":     kite_counts["hor_strong"],
                "median_confidence": (round(statistics.median(confidences), 4)
                                      if confidences else None),
                "profile_png":      f"{kite_dir_rel}/profile_plots/profile_{trc_id}.png",
                "profile_top3_png": f"{kite_dir_rel}/profile_plots/profile_top3_{trc_id}.png",
                # New kitehor-derived per-TRC aggregates (None on legacy data).
                "combined_class":   dominant_class,
                "combined_class_counts": dict(combined_class_counter) or None,
                "n_subrepeat":      n_subrepeat,
                "n_ssr":            n_ssr,
                # v0.12 TRC aggregates
                "prevalent_founder":   prev_founder,
                "prevalent_founder_n": prev_founder_n,
                "multiplicity_counts": dict(mult_counts) or None,
                "n_irregular":         n_irreg,
            }
        tarean_block = None
        if has_tarean:
            tarean_block = dict(tarean)
            tarean_block["variants"] = load_per_trc_variants(paths, trc_id)

        trcs.append({
            "id": trc_id,
            "index": _trc_sort_key(trc_id),
            "repeat_type":  repeat_type,
            "annotation":   annotation,
            "ssr_motif":    ssr_motif,
            "n_arrays":     len(arrays),
            "total_size":   sum(lengths),
            "min_array":    min(lengths) if lengths else None,
            "max_array":    max(lengths) if lengths else None,
            "median_array": int(statistics.median(lengths)) if lengths else None,
            "arrays":       arrays,
            "tarean":       tarean_block,
            "kite":         kite_block,
            "superfamily":  sf_by_trc.get(trc_id),
            "sf_peers":     sf_peers.get(trc_id, []),
            "sf_dotplot":   sf_dotplot_by_trc.get(trc_id),
        })

    # kitehor >=0.12.0 produces founder_period per array; the cascade-
    # driven combined_class is gone. Detect once so render functions can
    # swap out the class-driven roll-ups for the v0.12 counters.
    is_v012 = any(a.get("founder_period") is not None
                  for t in trcs for a in t["arrays"])
    return {
        "schema_version": __version__,
        "is_v012":        is_v012,
        "meta": {
            "prefix": prefix,
            "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
            "input_dir": str(input_dir),
        },
        "settings":      settings,
        "stats":         stats,
        "paths":         {k: (str(v.relative_to(input_dir)) if isinstance(v, Path) and v
                              else None)
                          for k, v in paths.items()},
        "seqid_lengths": seqid_lengths,
        "trcs":          trcs,
        "superfamilies": superfams,
    }


# ----------------------------------------------------------------------
# HTML page renderers
# ----------------------------------------------------------------------

def _trc_link(trc_id, text=None, prefix="trc/"):
    """Link to a TRC dashboard. Default prefix matches top-level pages."""
    return f'<a href="{prefix}{esc(trc_id)}.html">{esc(text or trc_id)}</a>'


def _rel_img(src_rel, width=None, alt="", src_prefix="../"):
    """An <img> whose src points back into the original output dir.

    src_prefix defaults to one `../` (right for top-level pages); dashboards
    pass `../../` since they sit one dir deeper."""
    if not src_rel:
        return ""
    w = f' width="{int(width)}"' if width else ""
    return f'<img src="{src_prefix}{esc(src_rel)}" alt="{esc(alt)}"{w}>'


def _tarean_asset_path(trc):
    """Relative path from report_v2/ to the per-TRC TAREAN image directory."""
    tdir = (trc.get("tarean") or {}).get("tarean_dir")
    return f"../{model_tarean_root}/{tdir}" if tdir else None  # not used; explicit per site


def _shell(ctx, title, active, run_meta, main_html):
    return page_shell(title, active, run_meta, main_html,
                      up=ctx["up"], root_href=ctx["root_href"],
                      legacy_href=ctx["legacy_href"],
                      assets_href=ctx["assets_href"])


def _sf_color_map(model):
    """Deterministic distinct colour per superfamily id (golden-angle hue).
    TRCs with no superfamily render grey."""
    sf_ids = sorted({sf["id"] for sf in model.get("superfamilies", [])})
    cmap = {}
    for i, sid in enumerate(sf_ids):
        hue = (i * 137.508) % 360
        cmap[sid] = f"hsl({hue:.0f}, 62%, 52%)"
    return cmap


def _overview_tick(v):
    """Axis tick label for a (positive) decade value."""
    if v >= 1:
        return f"{int(round(v)):,}"
    return f"{v:g}"


# ----------------------------------------------------------------------
# Index-level genome-distribution plot (chromosome ideogram showing
# every TRA across the run, with a click-to-highlight TRC selector on
# the left). Hand-rolled SVG; the highlight + search interactions live
# in tidecluster.js (`initIndexDistribution`).
#
# Distinct from the PER-TRC ideogram (`render_trc_distribution`) which
# only shows arrays of one TRC and lives on each TRC dashboard.
# ----------------------------------------------------------------------

INDEX_DIST_MIN_CONTIG = 1_000_000   # hide scaffolds smaller than this
INDEX_DIST_MAX_CONTIGS = 80         # hard cap so plant assemblies stay legible
INDEX_DIST_SVG_WIDTH   = 760
INDEX_DIST_LABEL_WIDTH = 170
INDEX_DIST_ROW_HEIGHT  = 16
INDEX_DIST_BAR_HEIGHT  = 9


def render_index_distribution(model, ctx):
    """Return a `<section>` HTML block with chromosome ideogram + TRC
    side-selector. Empty string if there are no contigs ≥ the threshold
    or no arrays — both happen on tiny test datasets.

    Each TRA <rect> carries data-trc / data-sf and a per-rectangle CSS
    variable `--sf` set to the superfamily colour. Default rendering is
    uniform muted grey; when the surrounding container's data-active-trc
    matches a rectangle's data-trc the JS adds `tc-idx-active` to it and
    CSS lifts it to its superfamily colour + outline."""
    trcs = model.get("trcs") or []
    seqid_lengths = model.get("seqid_lengths") or {}
    # All arrays, tagged with their TRC id and superfamily colour. We
    # collect by_seqid first so contig filtering is one pass.
    sfcol = _sf_color_map(model)
    by_seqid = {}
    for t in trcs:
        sf_color = sfcol.get(t.get("superfamily"), "#888888")
        for a in t.get("arrays") or []:
            sid = a.get("seqid")
            if not sid:
                continue
            by_seqid.setdefault(sid, []).append((t["id"], sf_color, a))
    if not by_seqid:
        return ""
    # Contig length resolution: prefer seqid_lengths.tsv; fall back to
    # the max array end on that contig so legacy runs without the
    # side-car still get a plot (with slightly under-reported lengths).
    contig_lengths = {sid: int(l) for sid, l in seqid_lengths.items()}
    for sid, items in by_seqid.items():
        if sid not in contig_lengths:
            contig_lengths[sid] = max(a["end"] for _, _, a in items)
    # Keep contigs ≥ threshold, longest first, capped at MAX_CONTIGS.
    kept = sorted(((sid, L) for sid, L in contig_lengths.items()
                   if L >= INDEX_DIST_MIN_CONTIG),
                  key=lambda x: -x[1])[:INDEX_DIST_MAX_CONTIGS]
    if not kept:
        return ""
    n_hidden = sum(1 for L in contig_lengths.values()
                   if L < INDEX_DIST_MIN_CONTIG)
    max_len = max(L for _, L in kept)
    label_w = INDEX_DIST_LABEL_WIDTH
    bar_w   = INDEX_DIST_SVG_WIDTH - label_w - 20
    row_h   = INDEX_DIST_ROW_HEIGHT
    bar_h   = INDEX_DIST_BAR_HEIGHT
    min_rect_w = 1.5
    height = len(kept) * row_h + 44

    parts = [
        f'<svg viewBox="0 0 {INDEX_DIST_SVG_WIDTH} {height}" '
        f'class="tc-idx-svg" style="max-width:100%;height:auto;" '
        f'data-active-trc="">'
    ]
    for i, (sid, length) in enumerate(kept):
        y = i * row_h + 4
        label = sid if len(sid) <= 22 else sid[:20] + "…"
        parts.append(
            f'<text x="4" y="{y + bar_h * 0.85:.0f}" font-size="11" '
            f'fill="currentColor">{esc(label)}</text>')
        parts.append(
            f'<text x="{label_w - 6}" y="{y + bar_h * 0.85:.0f}" '
            f'font-size="10" text-anchor="end" fill="currentColor" '
            f'opacity="0.55">{fmt_bp(length)}</text>')
        contig_px = (length / max_len) * bar_w
        parts.append(
            f'<rect x="{label_w}" y="{y}" width="{contig_px:.2f}" '
            f'height="{bar_h}" fill="#e9e9e9" stroke="none"/>')
        # Draw TRAs ordered by ascending length so the small ones land
        # on top and stay clickable on dense regions.
        items = sorted(by_seqid.get(sid, []),
                       key=lambda it: -(it[2]["end"] - it[2]["start"]))
        for trc_id, sf_color, a in items:
            x = label_w + (a["start"] / max_len) * bar_w
            w = max(min_rect_w,
                    ((a["end"] - a["start"]) / max_len) * bar_w)
            title = (f'{trc_id} · {esc(sid)}:{a["start"]:,}-{a["end"]:,} '
                     f'· {fmt_bp(a["end"] - a["start"])}')
            parts.append(
                f'<rect class="tc-idx-tra" x="{x:.2f}" y="{y:.1f}" '
                f'width="{w:.2f}" height="{bar_h}" '
                f'data-trc="{esc(trc_id)}" '
                f'style="--sf:{sf_color}" '
                f'data-title="{title}"><title>{title}</title></rect>')
    # Scale bar across the bottom.
    scale_y = len(kept) * row_h + 16
    parts.append(
        f'<line x1="{label_w}" y1="{scale_y}" x2="{label_w + bar_w}" '
        f'y2="{scale_y}" stroke="#888" stroke-width="0.5"/>')
    parts.append(
        f'<text x="{label_w}" y="{scale_y + 13}" font-size="10" '
        f'fill="currentColor">0</text>')
    parts.append(
        f'<text x="{label_w + bar_w}" y="{scale_y + 13}" font-size="10" '
        f'text-anchor="end" fill="currentColor">{fmt_bp(max_len)}</text>')
    parts.append("</svg>")
    svg_html = "\n".join(parts)

    # Sidebar: TRCs ordered by total coverage desc (= TRC numbering).
    # Each item carries its data-trc + sf colour swatch.
    up = ctx["up"]
    sf_label = {sf["id"]: sf for sf in model.get("superfamilies", [])}
    trc_items = []
    for t in sorted(trcs, key=lambda t: -(t.get("total_size") or 0)):
        n_tra = len(t.get("arrays") or [])
        if not n_tra:
            continue
        sf = t.get("superfamily")
        sf_html = (f'<span class="tc-idx-sf-chip" '
                   f'style="background:{sfcol[sf]}"></span>'
                   if sf else
                   '<span class="tc-idx-sf-chip tc-idx-sf-chip-empty"></span>')
        cov = fmt_bp(t.get("total_size") or 0)
        trc_items.append(
            f'<li class="tc-idx-trc-item" data-trc="{esc(t["id"])}" '
            f'title="{esc(t["id"])} · {n_tra} TRA · {cov}">'
            f'{sf_html}'
            f'<a class="tc-idx-trc-link" '
            f'href="{up}trc/{esc(t["id"])}.html">{esc(t["id"])}</a>'
            f'<span class="tc-idx-trc-meta">{n_tra}&nbsp;TRA · {cov}</span>'
            f'</li>')

    hidden_note = (f' <span class="tc-idx-hidden-note">'
                   f'{n_hidden} scaffold(s) &lt; '
                   f'{fmt_bp(INDEX_DIST_MIN_CONTIG)} hidden</span>'
                   if n_hidden else '')
    capped_note = (f' <span class="tc-idx-hidden-note">'
                   f'showing {INDEX_DIST_MAX_CONTIGS} longest contigs</span>'
                   if len(kept) == INDEX_DIST_MAX_CONTIGS else '')

    return f"""
    <section class="tc-idx-dist-wrap">
      <h2>TRC distribution across the assembly</h2>
      <p style="font-size:12px;color:var(--fg-muted);margin:4px 0 10px;">
        Each tick is a tandem-repeat array (TRA). Click a TRC on the
        left to highlight just its arrays — others fade. Click the
        active TRC again, or the <em>Show all</em> button, to reset.
        Contigs shorter than {fmt_bp(INDEX_DIST_MIN_CONTIG)} are
        hidden.{hidden_note}{capped_note}
      </p>
      <div class="tc-idx-dist">
        <aside class="tc-idx-sidebar">
          <div class="tc-idx-sidebar-head">
            <input type="search" class="tc-idx-trc-search"
                   placeholder="Filter TRCs…" aria-label="Filter TRCs">
            <button type="button" class="tc-btn tc-idx-trc-reset"
                    title="Clear selection">Show all</button>
          </div>
          <ul class="tc-idx-trc-list">
            {"".join(trc_items)}
          </ul>
        </aside>
        <div class="tc-idx-plot">{svg_html}</div>
      </div>
    </section>"""


def render_cluster_overview(model, ctx):
    """Interactive (hover + click) bubble chart summarising the whole run:
    x = median per-array KITE monomer size (log bp), y = TRC coverage
    (log Mbp), dot radius ∝ number of arrays, dot colour = superfamily.
    Hand-rolled inline SVG; the floating tooltip is the same one the
    genome-distribution ideograms use (see tidecluster.js)."""
    import math
    trc_prefix = ctx["trc_link_prefix"]
    pts, n_excluded = [], 0
    for t in model["trcs"]:
        m1s = [a["m1"] for a in t["arrays"] if a.get("m1")]
        cov_bp = t.get("total_size") or 0
        if not m1s or cov_bp <= 0:
            n_excluded += 1
            continue
        kite = t.get("kite") or {}
        cc = kite.get("combined_class_counts") or {}
        pts.append({
            "id":   t["id"],
            "x":    statistics.median(m1s),
            "y":    cov_bp / 1e6,
            "n":    t["n_arrays"],
            "cov":  cov_bp,
            "hor":  cc.get("hor", 0) + cc.get("hor_with_ssr", 0),
            "cls":  kite.get("combined_class"),
            "sf":   t.get("superfamily"),
            "ann":  t.get("annotation"),
        })
    if not pts:
        return ""

    W, H = 860, 480
    ml, mr, mt, mb = 70, 20, 20, 58
    px0, px1 = ml, W - mr          # x pixel range (left→right)
    py0, py1 = mt, H - mb          # y pixel range (top→bottom)

    def _scale(vals, lo_px, hi_px, invert):
        lo = math.floor(math.log10(min(vals)))
        hi = math.ceil(math.log10(max(vals)))
        if hi <= lo:
            hi = lo + 1
        span = hi - lo

        def f(v):
            tt = (math.log10(v) - lo) / span
            return hi_px - tt * (hi_px - lo_px) if invert else lo_px + tt * (hi_px - lo_px)
        return f, lo, hi

    xf, xlo, xhi = _scale([p["x"] for p in pts], px0, px1, invert=False)
    yf, ylo, yhi = _scale([p["y"] for p in pts], py0, py1, invert=True)
    sfcol = _sf_color_map(model)

    def radius(n):
        return max(3.0, min(22.0, 2.2 * math.sqrt(n)))

    parts = [f'<svg viewBox="0 0 {W} {H}" class="tc-scatter-svg" '
             f'style="max-width:100%;height:auto;" '
             f'role="img" aria-label="Cluster overview scatter">']
    # Gridlines + ticks
    for k in range(xlo, xhi + 1):
        X = xf(10 ** k)
        parts.append(f'<line x1="{X:.1f}" y1="{py0}" x2="{X:.1f}" y2="{py1}" '
                     f'class="tc-grid"/>')
        parts.append(f'<text x="{X:.1f}" y="{py1 + 16}" text-anchor="middle" '
                     f'class="tc-axis-lbl">{_overview_tick(10 ** k)}</text>')
    for k in range(ylo, yhi + 1):
        Y = yf(10 ** k)
        parts.append(f'<line x1="{px0}" y1="{Y:.1f}" x2="{px1}" y2="{Y:.1f}" '
                     f'class="tc-grid"/>')
        parts.append(f'<text x="{px0 - 8}" y="{Y + 3:.1f}" text-anchor="end" '
                     f'class="tc-axis-lbl">{_overview_tick(10 ** k)}</text>')
    # Axis titles
    parts.append(f'<text x="{(px0 + px1) / 2:.0f}" y="{H - 6}" '
                 f'text-anchor="middle" class="tc-axis-title">'
                 f'Median monomer size per TRC (bp, log)</text>')
    parts.append(f'<text x="16" y="{(py0 + py1) / 2:.0f}" '
                 f'text-anchor="middle" class="tc-axis-title" '
                 f'transform="rotate(-90 16 {(py0 + py1) / 2:.0f})">'
                 f'TRC coverage (Mbp, log)</text>')
    # Points — larger dots first so small ones stay clickable on top.
    for p in sorted(pts, key=lambda p: -p["n"]):
        cx, cy, r = xf(p["x"]), yf(p["y"]), radius(p["n"])
        fill = sfcol.get(p["sf"], "#8a8a8a")
        sf_line = (f'<br>superfamily: SF {esc(p["sf"])}' if p["sf"]
                   else '<br>superfamily: unassigned')
        cls_line = f'<br>class: {esc(class_label(p["cls"]))}' if p["cls"] else ''
        ann_line = f'<br>annotation: {esc(p["ann"])}' if p.get("ann") else ''
        title = (f'<strong>{esc(p["id"])}</strong>'
                 f'<br>median monomer: {esc(int(round(p["x"])))} bp'
                 f'<br>coverage: {p["y"]:.3g} Mbp ({p["cov"]:,} bp)'
                 f'<br>arrays (TRA): {esc(p["n"])}'
                 f'<br>HOR arrays: {esc(p["hor"])}'
                 f'{cls_line}{sf_line}{ann_line}')
        plain = f'{p["id"]} · {int(round(p["x"]))} bp · {p["y"]:.3g} Mbp · {p["n"]} TRA'
        parts.append(
            f'<a href="{trc_prefix}{esc(p["id"])}.html">'
            f'<circle class="tc-point" cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
            f'fill="{fill}" data-html="1" data-title="{title}">'
            f'<title>{esc(plain)}</title></circle></a>')
    parts.append('</svg>')

    # Legends: dot-size reference + superfamily swatches.
    size_refs = [n for n in (1, 10, 100) if n <= max(p["n"] for p in pts)] or [1]
    size_legend = "".join(
        f'<span class="tc-size-ref"><svg width="{2*radius(n)+2:.0f}" '
        f'height="{2*radius(n)+2:.0f}" style="vertical-align:middle">'
        f'<circle cx="{radius(n)+1:.1f}" cy="{radius(n)+1:.1f}" r="{radius(n):.1f}" '
        f'class="tc-point" fill="#8a8a8a"/></svg> {n}</span>'
        for n in size_refs)
    sf_present = sorted({p["sf"] for p in pts if p["sf"]})
    sf_swatches = "".join(
        f'<span class="tc-sf-leg"><span class="tc-sf-swatch" '
        f'style="background:{sfcol[s]}"></span>SF {esc(s)}</span>'
        for s in sf_present)
    sf_legend = (f'<span class="tc-sf-leg"><span class="tc-sf-swatch" '
                 f'style="background:#8a8a8a"></span>unassigned</span>{sf_swatches}')

    excl = (f' {n_excluded} TRC(s) without a KITE monomer estimate '
            f'(e.g. below the TAREAN size threshold) are not shown.'
            if n_excluded else "")
    return f"""
    <section>
      <h2>Cluster overview</h2>
      <p style="font-size:12px;color:var(--fg-muted)">
        Each dot is a TRC. <strong>x</strong> = median per-array monomer
        size (KITE m₁, log); <strong>y</strong> = TRC coverage (log Mbp);
        <strong>dot size</strong> ∝ number of arrays (TRA);
        <strong>colour</strong> = superfamily (grey = unassigned). Hover a
        dot for details, click to open the TRC.{excl}
      </p>
      <div class="tc-scatter-wrap">{''.join(parts)}</div>
      <div class="tc-scatter-legend">
        <span class="tc-leg-grp">arrays:&nbsp;{size_legend}</span>
        <span class="tc-leg-grp">{sf_legend}</span>
      </div>
    </section>"""


def render_index(model, out_path, run_meta, ctx):
    settings = model.get("settings", {}) or {}
    stats    = model.get("stats", {}) or {}
    rows_settings = []
    for k in ("command", "fasta", "original_fasta", "prefix", "min_length",
              "min_total_length", "cpu", "no_dust", "tidehunter_arguments", "long",
              "library"):
        if k in settings:
            rows_settings.append(
                f'<dt>{esc(k)}</dt><dd>{esc(settings[k])}</dd>')
    n_tarean    = sum(1 for t in model["trcs"] if t["tarean"])
    n_sf_peers  = sum(1 for t in model["trcs"] if t["superfamily"])
    up = ctx["up"]
    is_v012 = bool(model.get("is_v012"))

    # Per-array aggregates: counts + SSR motif tally weighted by bp.
    # Fallback-founder count is intentionally not surfaced here (see TRC card).
    n_hor = n_subrep = n_ssr_signal = n_ssr_dom = 0
    ssr_bp = Counter()
    total_ssr_bp = 0
    for t in model["trcs"]:
        for a in t["arrays"]:
            if (a.get("multiplicity") or 1) > 1:
                n_hor += 1
            if a.get("subrepeat_1_tier"):
                n_subrep += 1
            motif = (a.get("ssr_dominant_motif") or "").strip()
            pct   = a.get("ssr_total_coverage_pct") or 0
            if motif and motif != "NA":
                n_ssr_signal += 1
                alen = (a.get("array_length")
                        or (a.get("end") or 0) - (a.get("start") or 0))
                bp = int((alen or 0) * float(pct) / 100.0)
                ssr_bp[motif]  += bp
                total_ssr_bp   += bp
                if float(pct) >= 50.0:
                    n_ssr_dom += 1
    top_ssr = ssr_bp.most_common(3)
    # Wrap each motif so the _kv helper's "contains <" detection skips
    # escaping (otherwise &nbsp; entities would be double-encoded).
    top_ssr_html = (", ".join(
        f'<span class="tc-ssr-motif">{esc(m)}&nbsp;({fmt_bp(bp)})</span>'
        for m, bp in top_ssr) if top_ssr else "")

    # ------------------------------------------------------------------
    # RUN SUMMARY — sectioned: Input / Clusters / Arrays. REPORTS tile
    # was dropped (its rows are folded in here or reachable from the
    # bubble chart / per-tab links). Fallback-founder line removed by
    # request; per-TRC card still flags it inline.
    # ------------------------------------------------------------------
    def _kv(rows):
        # Each row: (label, ready-to-render html or plain value). Plain
        # values are esc()'d; explicit HTML strings (links, motifs) pass
        # through. None / "" rows are skipped.
        out = []
        for label, val in rows:
            if val in (None, ""):
                continue
            html = val if (isinstance(val, str) and "<" in val) else esc(val)
            out.append(f'<dt>{esc(label)}</dt><dd>{html}</dd>')
        return "".join(out)

    input_rows = [
        ("sequence length", fmt_bp(stats.get("input_sequence_length"))),
    ]
    cluster_rows = [
        ("total",                  stats.get("n_trcs_total")),
        ("above length threshold", stats.get("n_trcs_above_threshold")),
        ("with TAREAN consensus",
            f'<a href="{up}tarean.html">{n_tarean}</a>'),
        ("total length",           fmt_bp(stats.get("total_tr_length"))),
        ("superfamilies",
            f'<a href="{up}superfamilies.html">{len(model["superfamilies"])}</a>'
            f' ({n_sf_peers} TRCs grouped)'),
    ]
    array_rows = [("total (KITE-analysed)", stats.get("n_tras"))]
    if is_v012:
        array_rows += [
            ("HOR (×k≥2)",           n_hor),
            ("with subrepeat",       n_subrep),
            ("with SSR signal",      n_ssr_signal),
            ("SSR-dominant (≥50 %)", n_ssr_dom),
        ]
    else:
        # Legacy roll-up for pre-v0.12 dirs.
        cls_total = class_totals(model)
        array_rows += [
            ("HOR arrays",
                f'{class_badge("hor")} {cls_total.get("hor", 0)}'),
            ("array classes", class_mix_cell(cls_total)),
        ]
    if total_ssr_bp:
        array_rows += [
            ("SSR total length", fmt_bp(total_ssr_bp)),
            ("top SSR motifs",   top_ssr_html),
        ]

    run_summary_body = (
        f'<div class="tc-kv-section">Input</div>'
        f'<dl class="tc-kv">{_kv(input_rows)}</dl>'
        f'<div class="tc-kv-section">Clusters (TRC)</div>'
        f'<dl class="tc-kv">{_kv(cluster_rows)}</dl>'
        f'<div class="tc-kv-section">Arrays (TRA)</div>'
        f'<dl class="tc-kv">{_kv(array_rows)}</dl>')

    cards = f"""
    <section class="tc-cards tc-cards-2">
      <div class="tc-card"><div class="tc-card-title">Run summary</div>
        {run_summary_body}</div>
      <div class="tc-card"><div class="tc-card-title">Run settings</div>
        <dl class="tc-kv">{"".join(rows_settings)}</dl></div>
    </section>"""
    overview_html = load_section("overview")
    credits_html  = load_section("credits")
    cluster_chart = render_cluster_overview(model, ctx)
    distribution_chart = render_index_distribution(model, ctx)
    main = f"""
    <h1>TideCluster report — {esc(model["meta"]["prefix"])}</h1>
    {cards}
    {cluster_chart}
    {distribution_chart}
    <section class="tc-prose">{overview_html}</section>
    <hr style="margin:28px 0; border:0; border-top:1px solid var(--border);">
    <section class="tc-prose">{credits_html}</section>
    """
    Path(out_path).write_text(_shell(ctx, "Summary", "index", run_meta, main))


def _trc_type_label(t):
    """Human-readable type label for the merged TAREAN table. v0.12+:
    only two values — SSR or TR. (The earlier
    “TR (below threshold)” nuance is now implicit: rows with empty
    TAREAN columns are below threshold; the Type column doesn't repeat
    that.)"""
    if (t.get("repeat_type") or "").upper() == "SSR" or t.get("ssr_motif"):
        return "SSR"
    return "TR"


def _array_size_cell(t):
    """Stacked min/median/max cell, matching the original TAREAN report."""
    mn, md, mx = t.get("min_array"), t.get("median_array"), t.get("max_array")
    if mn is None:
        return ""
    return (f"min: {fmt_bp(mn)}<br>median: {fmt_bp(md)}<br>max: {fmt_bp(mx)}")


def _thumb(src_rel, kind, src_prefix="../", alt=""):
    """Fixed-size thumbnail wrapper; full PNG opens on click.

    kind is one of "graph", "logo", "profile" → maps to a CSS size class."""
    if not src_rel:
        return ""
    klass = {"graph": "tc-thumb-graph",
             "logo":  "tc-thumb-logo",
             "profile": "tc-thumb-profile"}.get(kind, "tc-thumb-graph")
    url = f"{src_prefix}{esc(src_rel)}"
    return (f'<a class="tc-thumb {klass}" href="{url}" target="_blank" '
            f'rel="noopener"><img src="{url}" alt="{esc(alt)}"></a>')


def _render_tarean_row(t, ctx):
    src_prefix = ctx["src_prefix"]
    up = ctx["up"]
    ta = t["tarean"] or {}
    kite = t["kite"] or {}
    # v0.12: the H/S/X structural-signal column was dropped from this
    # table; HOR / subrepeat / SSR counts live on the per-TRC dashboard
    # instead. Legacy dirs (pre-v0.12) still get the class-mix pill row
    # so older reports rerender unchanged.
    if ctx.get("is_v012"):
        class_cell = None
    else:
        class_cell = (class_mix_cell(Counter(kite.get("combined_class_counts") or {}))
                      if kite else "")
    sf_cell = (f'<a href="{up}superfamilies.html#sf-{t["superfamily"]}">SF {t["superfamily"]}</a>'
               if t["superfamily"] else "")
    cons_short = ""
    if ta.get("consensus"):
        cons_short = esc(ta["consensus"][:80])
        if len(ta["consensus"]) > 80:
            cons_short += "…"
    graph_thumb = _thumb(ta.get("graph_png"), "graph", src_prefix=src_prefix, alt="k-mer graph")
    logo_thumb  = _thumb(ta.get("logo_png"),  "logo",  src_prefix=src_prefix, alt="logo")
    return (
        "<tr>"
        f'<td data-order="{t["index"]}">{_trc_link(t["id"], prefix=ctx["trc_link_prefix"])}</td>'
        f'<td>{_trc_type_label(t)}</td>'
        f'<td data-order="{t["total_size"]}">{fmt_bp(t["total_size"])}</td>'
        f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
        f'<td data-order="{t.get("max_array") or 0}">{_array_size_cell(t)}</td>'
        f'<td data-order="{ta.get("monomer_length") or 0}">{esc(ta.get("monomer_length"))}</td>'
        f'<td data-order="{ta.get("total_score") or 0}">'
            f'{esc(round(ta["total_score"], 4)) if ta.get("total_score") is not None else ""}</td>'
        f'<td>{esc((t.get("ssr_motif") or "").replace("%25","%"))}</td>'
        f'<td>{esc(t.get("annotation") or "")}</td>'
        + (f'<td>{class_cell}</td>' if class_cell is not None else '')
        + f'<td>{sf_cell}</td>'
        f'<td>{graph_thumb}</td>'
        f'<td>{logo_thumb}</td>'
        f'<td><code style="font-size:11px;">{cons_short}</code></td>'
        "</tr>")


def render_tarean(model, out_path, run_meta, ctx):
    """Merged TAREAN / All-TRCs table. Every TRC from clustering.gff3
    appears here (TAREAN-analysed, below-threshold, SSR). TAREAN-specific
    columns are empty for rows where TAREAN was not performed."""
    # Thread the v0.12 flag into ctx so _render_tarean_row can pick the
    # right structural-signal cell (no in-place mutation: shallow copy).
    row_ctx = dict(ctx); row_ctx["is_v012"] = model.get("is_v012")
    rows = "".join(_render_tarean_row(t, row_ctx) for t in model["trcs"])
    n_total   = len(model["trcs"])
    n_tarean  = sum(1 for t in model["trcs"] if t["tarean"])
    n_ssr     = sum(1 for t in model["trcs"] if _trc_type_label(t) == "SSR")
    n_below   = n_total - n_tarean - n_ssr
    # v0.12: the H/S/X column is gone (HOR / subrep / SSR live on the
    # per-TRC dashboard). Legacy dirs keep the class-mix column header.
    signal_header = ("" if model.get("is_v012")
                     else "<th>Class mix</th>")
    body = f"""
    <h1>Tandem repeat clusters (TAREAN)</h1>
    <section class="tc-prose">{load_section("tarean")}</section>
    <p>{n_total} TRCs total — {n_tarean} analysed by TAREAN,
       {n_ssr} SSR clusters, {n_below} below the TAREAN size threshold.
       Click a TRC to open its per-TRC dashboard.
       Graph and logo thumbnails open the full image.</p>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="50"
           data-order='[[0,"asc"]]'>
      <thead><tr>
        <th>TRC</th>
        <th>Type</th>
        <th>Total size</th>
        <th>Arrays</th>
        <th>Array sizes<br>(min / med / max)</th>
        <th>Monomer<br>(TAREAN)</th>
        <th>Score<br>(TAREAN)</th>
        <th>SSR motif</th>
        <th>Annotation</th>
        {signal_header}
        <th>Superfamily</th>
        <th>Graph</th>
        <th>Logo</th>
        <th>Consensus<br>(preview)</th>
      </tr></thead>
      <tbody>{rows}</tbody>
    </table>
    </div>
    """
    Path(out_path).write_text(_shell(ctx, "TAREAN", "tarean", run_meta, body))
    # Old (pre-merge) page no longer makes sense; clean it up so the nav
    # stays consistent when rerender runs over an existing report dir.
    old = Path(out_path).parent / "trcs.html"
    if old.exists():
        old.unlink()


def render_kite(model, out_path, run_meta, ctx):
    """KITE overview page. For kitehor >=0.12.0 (no combined_class) the
    per-array distribution bar is replaced with per-TRC structural
    counters: arrays with HOR multiplicity, with subrepeat candidates,
    with SSR, with fallback founders."""
    is_v012 = model.get("is_v012")
    rows = []
    tot_hor = tot_subrep = tot_ssr = tot_fb = 0
    for t in model["trcs"]:
        if not t["kite"]:
            continue
        k = t["kite"]
        profile_thumb = _thumb(k.get("profile_top3_png"), "profile",
                               src_prefix=ctx["src_prefix"], alt="KITE profile")
        tarean_mon = (t["tarean"] or {}).get("monomer_length")
        kite_mon   = k.get("monomer_primary")
        if is_v012:
            n_hor    = sum(1 for a in t["arrays"]
                           if (a.get("multiplicity") or 1) > 1)
            n_subrep = sum(1 for a in t["arrays"] if a.get("subrepeat_1_tier"))
            n_ssr    = sum(1 for a in t["arrays"]
                           if (a.get("ssr_dominant_motif") or "NA") != "NA")
            n_fb     = sum(1 for a in t["arrays"] if a.get("founder_fallback"))
            tot_hor    += n_hor; tot_subrep += n_subrep
            tot_ssr    += n_ssr; tot_fb     += n_fb
            rows.append(
                "<tr>"
                f'<td data-order="{t["index"]}">{_trc_link(t["id"], prefix=ctx["trc_link_prefix"])}</td>'
                f'<td data-order="{kite_mon or 0}">{esc(kite_mon)}</td>'
                f'<td data-order="{tarean_mon or 0}">{esc(tarean_mon)}</td>'
                f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
                f'<td data-order="{t["total_size"]}">{fmt_bp(t["total_size"])}</td>'
                f'<td data-order="{n_hor}">{n_hor or ""}</td>'
                f'<td data-order="{n_subrep}">{n_subrep or ""}</td>'
                f'<td data-order="{n_ssr}">{n_ssr or ""}</td>'
                f'<td data-order="{n_fb}">{("*" * min(n_fb,3)) if n_fb else ""}</td>'
                f'<td>{profile_thumb}</td>'
                "</tr>")
        else:
            med_conf  = k.get("median_confidence")
            med_str   = f"{med_conf:.3f}" if med_conf is not None else ""
            cc_counts = k.get("combined_class_counts") or {}
            dominant  = k.get("combined_class")
            rows.append(
                "<tr>"
                f'<td data-order="{t["index"]}">{_trc_link(t["id"], prefix=ctx["trc_link_prefix"])}</td>'
                f'<td data-order="{kite_mon or 0}">{esc(kite_mon)}</td>'
                f'<td data-order="{tarean_mon or 0}">{esc(tarean_mon)}</td>'
                f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
                f'<td data-order="{t["total_size"]}">{fmt_bp(t["total_size"])}</td>'
                f'<td data-order="{esc(dominant or "")}">{class_badge(dominant)}</td>'
                f'<td>{class_mix_cell(Counter(cc_counts))}</td>'
                f'<td data-order="{med_conf or 0}">{med_str}</td>'
                f'<td>{profile_thumb}</td>'
                "</tr>")

    if is_v012:
        head = ('<tr>'
                '<th>TRC</th>'
                '<th>Monomer<br>(KITE&nbsp;primary)</th>'
                '<th>Monomer<br>(TAREAN)</th>'
                '<th>Arrays</th>'
                '<th>Total size</th>'
                '<th>HOR<br>(×k≥2)</th>'
                '<th>Subrep</th>'
                '<th>SSR</th>'
                '<th>Fallback</th>'
                '<th>Profile</th>'
                '</tr>')
        summary = (
            '<h3>Structural signal counts across all analysed arrays</h3>'
            '<p style="font-size:12px;color:var(--fg-muted)">'
            f'HOR (×k≥2): <strong>{tot_hor}</strong> &nbsp;·&nbsp; '
            f'subrepeat: <strong>{tot_subrep}</strong> &nbsp;·&nbsp; '
            f'SSR: <strong>{tot_ssr}</strong> &nbsp;·&nbsp; '
            f'fallback founder: <strong>{tot_fb}</strong>. '
            'Subrepeat counts use the strict HIGH/LIKELY tiers '
            '(see <em>Subrepeat</em> column legend on any TRC page).'
            '</p>'
            '<h2>Per-TRC structural summary</h2>'
            '<p style="font-size:12px;color:var(--fg-muted)">'
            '<strong>Monomer (KITE)</strong> is the most-frequent primary '
            'rescore founder period across the arrays of a given TRC. '
            '<strong>HOR (×k≥2)</strong> counts arrays whose strongest '
            'period is an integer multiple of the founder; '
            '<strong>Subrep</strong> counts arrays with at least one '
            'HIGH/LIKELY subrepeat candidate; <strong>SSR</strong> counts '
            'arrays with a dominant short-motif tandem; '
            '<strong>Fallback</strong> marks arrays where rescore could '
            'not assign a founder and TideCluster fell back to the top-'
            'scored kite peak.'
            '</p>')
    else:
        cls_total = class_totals(model)
        total = sum(cls_total.values()) or 1
        bar = ""; legend_bits = []
        for raw in CANONICAL_CLASS_ORDER:
            n = cls_total.get(raw, 0)
            if not n: continue
            pct = 100 * n / total
            css = COMBINED_CLASS_CSS[raw]
            bar += (f'<div class="tc-bar-seg class-bar-{css}" '
                    f'style="width:{pct:.2f}%" '
                    f'title="{esc(class_label(raw))}: {n}">'
                    f'{n if pct > 3 else ""}</div>')
            legend_bits.append(f'{class_badge(raw)} {n}')
        head = ('<tr>'
                '<th>TRC</th><th>Monomer<br>(KITE&nbsp;primary)</th>'
                '<th>Monomer<br>(TAREAN)</th><th>Arrays</th>'
                '<th>Total size</th><th>Cluster<br>class</th>'
                '<th>Class mix</th><th>Median<br>conf.</th><th>Profile</th>'
                '</tr>')
        summary = (
            '<h3>Array class distribution across all analysed arrays</h3>'
            f'<div class="tc-bar" title="Total arrays: {total}">{bar}</div>'
            '<p style="font-size:12px;color:var(--fg-muted)">'
            + " · ".join(legend_bits) +
            f' (total {total} arrays).</p>'
            '<h2>Per-TRC class summary</h2>')

    body = f"""
    <h1>K-mer interval tandem repeat estimation (KITE)</h1>
    <section class="tc-prose">{load_section("kite")}</section>
    {summary}
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[0,"asc"]]'>
      <thead>{head}</thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
    </div>
    """
    Path(out_path).write_text(_shell(ctx, "KITE", "kite", run_meta, body))


def render_superfamilies(model, out_path, run_meta, ctx):
    out_path = Path(out_path)
    if not model["superfamilies"]:
        body = (f'<h1>TRC Superfamilies</h1>'
                f'<section class="tc-prose">{load_section("superfamilies")}</section>'
                f'<p class="tc-callout">No TRC superfamilies were identified '
                f'(either insufficient TRCs reached the clustering step or none shared '
                f'significant consensus similarity).</p>')
        out_path.write_text(_shell(ctx, "Superfamilies", "sf", run_meta, body))
        return
    src_prefix = ctx["src_prefix"]
    sections = []
    for sf in model["superfamilies"]:
        peers = ", ".join(_trc_link(t, prefix=ctx["trc_link_prefix"]) for t in sf["trcs"])
        img = (f'<div class="tc-fig"><a href="{src_prefix}{esc(sf["dotplot"])}">'
               f'<img src="{src_prefix}{esc(sf["dotplot"])}" width="500" alt="dotplot"></a>'
               f'<div class="tc-fig-caption">pairwise dotplot of superfamily '
               f'{sf["id"]} consensus sequences</div></div>'
               if sf["dotplot"] else "")
        sections.append(
            f'<section id="sf-{sf["id"]}">'
            f'<h2>Superfamily {sf["id"]} ({len(sf["trcs"])} TRCs)</h2>'
            f'<p>{peers}</p>{img}</section>')
    body = (f'<h1>TRC superfamilies</h1>'
            f'<section class="tc-prose">{load_section("superfamilies")}</section>'
            f'<p>{len(model["superfamilies"])} superfamilies identified.</p>'
            + "".join(sections))
    out_path.write_text(_shell(ctx, "Superfamilies", "sf", run_meta, body))


# ----------------------------------------------------------------------
# Per-TRC genome-distribution visualisation (ideogram + minor-contig
# table). Pure inline SVG; no JS dep beyond the existing DataTables
# enhancement on the minor table. See the "TRC distribution" section
# of the docs for the design rationale.
# ----------------------------------------------------------------------

TRC_DIST_MAJOR_MAX_COUNT  = 50          # all contigs kept as major if <= this
TRC_DIST_MAJOR_MIN_LENGTH = 1_000_000   # fallback length threshold above that
TRC_DIST_SVG_WIDTH        = 920
TRC_DIST_LABEL_WIDTH      = 190
TRC_DIST_ROW_HEIGHT       = 22
TRC_DIST_BAR_HEIGHT       = 12
TRC_DIST_MINI_WIDTH       = 180
TRC_DIST_MINI_HEIGHT      = 12

# SVG fill palette — saturated variants of the pastel badge colours so
# small rectangles (down to ~2 px) stay legible without needing a
# stroke. The badge CSS keeps its softer palette for text legibility.
# Keyed by kitehor combined_class. Unresolved gets a flag colour.
_CLASS_FILL = {
    "hor":               "#2e8b2e",  # forest green
    "hor_with_ssr":      "#2f9e78",  # green-teal (HOR + SSR)
    "tr":                "#8a8a8a",  # neutral grey (baseline)
    "tr_with_subrepeat": "#2bb3a3",  # teal
    "tr_with_ssr":       "#e88f00",  # amber
    "pure_ssr":          "#8e5bd0",  # purple
    "unresolved":        "#c0392b",  # muted red — needs attention
    "tr_with_nested_tr": "#3a7bd5",  # legacy 0.9.x — blue
}


def _class_fill(arr):
    cc = arr.get("combined_class")
    if cc:
        return _CLASS_FILL.get(cc, "#cccccc")
    # v0.12 fallback: colour by the dominant structural signal so the
    # ideogram still carries some semantic information without the
    # combined_class cascade.
    if (arr.get("multiplicity") or 1) > 1:
        return _CLASS_FILL["hor"]                # ×k≥2 → HOR green
    if arr.get("subrepeat_1_tier"):
        return _CLASS_FILL["tr_with_subrepeat"]  # teal
    if (arr.get("ssr_dominant_motif") or "NA") != "NA":
        return _CLASS_FILL["pure_ssr"]           # purple
    if arr.get("founder_fallback"):
        return _CLASS_FILL["unresolved"]         # muted red — needs attention
    return _CLASS_FILL["tr"]                     # neutral grey baseline


def _array_title(arr):
    cc = arr.get("combined_class")
    if cc:
        label = class_label(cc)
    else:
        # v0.12 label assembled from structural signals.
        parts = []
        mult = arr.get("multiplicity") or 1
        if mult > 1: parts.append(f"HOR ×{mult}")
        if arr.get("subrepeat_1_tier"): parts.append("subrepeat")
        if (arr.get("ssr_dominant_motif") or "NA") != "NA":
            parts.append("SSR")
        if arr.get("founder_fallback"): parts.append("fallback")
        label = " · ".join(parts) if parts else (arr.get("hor_status") or "—")
    span = f"{arr['seqid']}:{arr['start']:,}-{arr['end']:,}"
    length = f"{fmt_bp(arr['end'] - arr['start'])}"
    return f"{span} · {length} · {label}"


def _render_ideogram(majors, by_seqid):
    """Return inline SVG for a linear ideogram stack of major contigs.

    `majors` is a list of (seqid, length). Each TRA in by_seqid[seqid]
    is drawn as a coloured rectangle whose width is the array length
    scaled to the shared x-axis (= max major length)."""
    if not majors: return ""
    max_len = max(l for _, l in majors)
    label_w = TRC_DIST_LABEL_WIDTH
    bar_w   = TRC_DIST_SVG_WIDTH - label_w - 20
    row_h   = TRC_DIST_ROW_HEIGHT
    bar_h   = TRC_DIST_BAR_HEIGHT
    min_rect_w = 2.0
    height = len(majors) * row_h + 50
    parts = [
        f'<svg viewBox="0 0 {TRC_DIST_SVG_WIDTH} {height}" '
        f'class="tc-ideogram-svg" style="max-width:100%;height:auto;">'
    ]
    for i, (sid, length) in enumerate(majors):
        y = i * row_h + 6
        # Label (seqid + length, truncated if the name is long)
        label = sid if len(sid) <= 22 else sid[:20] + "…"
        parts.append(
            f'<text x="4" y="{y + bar_h * 0.75:.0f}" font-size="11" '
            f'fill="currentColor">{esc(label)}</text>'
        )
        parts.append(
            f'<text x="{label_w - 6}" y="{y + bar_h * 0.75:.0f}" font-size="10" '
            f'text-anchor="end" fill="currentColor" opacity="0.65">{fmt_bp(length)}</text>'
        )
        # Backing contig bar
        contig_px = (length / max_len) * bar_w
        parts.append(
            f'<rect x="{label_w}" y="{y}" width="{contig_px:.2f}" '
            f'height="{bar_h}" fill="#e6e6e6" stroke="none"/>'
        )
        # Arrays — sorted by ascending confidence so stronger ones draw on top.
        # No stroke: the fill has been darkened so small rectangles stay
        # prominent without a border muting their colour. Empty list for
        # major scaffolds that carry no arrays of this TRC — the backing
        # contig bar stays visible (absence is informative).
        arrs = sorted(by_seqid.get(sid, []),
                      key=lambda a: (a.get("hor_confidence") or 0))
        for arr in arrs:
            x_start = label_w + (arr["start"] / max_len) * bar_w
            w = max(min_rect_w, ((arr["end"] - arr["start"]) / max_len) * bar_w)
            parts.append(
                f'<rect class="tc-tra" x="{x_start:.2f}" y="{y}" '
                f'width="{w:.2f}" height="{bar_h}" '
                f'fill="{_class_fill(arr)}" stroke="none" '
                f'data-title="{esc(_array_title(arr))}">'
                f'<title>{esc(_array_title(arr))}</title></rect>'
            )
    # Scale bar
    scale_y = len(majors) * row_h + 18
    parts.append(
        f'<line x1="{label_w}" y1="{scale_y}" x2="{label_w + bar_w}" '
        f'y2="{scale_y}" stroke="#888" stroke-width="0.5"/>')
    parts.append(
        f'<text x="{label_w}" y="{scale_y + 14}" font-size="10" '
        f'fill="currentColor">0</text>')
    parts.append(
        f'<text x="{label_w + bar_w}" y="{scale_y + 14}" font-size="10" '
        f'text-anchor="end" fill="currentColor">{fmt_bp(max_len)}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def _render_minor_table(minors, by_seqid):
    """DataTables row per minor contig with inline mini-SVG + HOR mix."""
    if not minors: return ""
    rows = []
    for sid, length in minors:
        arrs = by_seqid[sid]
        n = len(arrs)
        total_arr_len = sum(a["end"] - a["start"] for a in arrs)
        class_counts = Counter(a.get("combined_class") for a in arrs
                               if a.get("combined_class"))
        mix_cells = class_mix_cell(class_counts)
        mini = [f'<svg viewBox="0 0 {TRC_DIST_MINI_WIDTH} {TRC_DIST_MINI_HEIGHT}" '
                f'class="tc-mini-svg" style="width:{TRC_DIST_MINI_WIDTH}px;'
                f'height:{TRC_DIST_MINI_HEIGHT}px;">']
        mini.append(
            f'<rect x="0" y="4" width="{TRC_DIST_MINI_WIDTH}" height="4" '
            f'fill="#e6e6e6" stroke="none"/>')
        for arr in sorted(arrs, key=lambda a: (a.get("hor_confidence") or 0)):
            x_start = (arr["start"] / max(length, 1)) * TRC_DIST_MINI_WIDTH
            w = max(1.5, ((arr["end"] - arr["start"]) / max(length, 1))
                    * TRC_DIST_MINI_WIDTH)
            mini.append(
                f'<rect class="tc-tra" x="{x_start:.2f}" y="2" '
                f'width="{w:.2f}" height="8" '
                f'fill="{_class_fill(arr)}" stroke="none" '
                f'data-title="{esc(_array_title(arr))}">'
                f'<title>{esc(_array_title(arr))}</title></rect>')
        mini.append("</svg>")
        rows.append(
            "<tr>"
            f'<td>{esc(sid)}</td>'
            f'<td data-order="{length}">{fmt_bp(length)}</td>'
            f'<td data-order="{n}">{n}</td>'
            f'<td data-order="{total_arr_len}">{fmt_bp(total_arr_len)}</td>'
            f'<td>{mix_cells}</td>'
            f'<td>{"".join(mini)}</td>'
            "</tr>"
        )
    return f"""
    <h3>Minor contigs carrying arrays of this TRC ({len(minors)})</h3>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="10"
           data-order='[[1,"desc"]]'>
      <thead><tr>
        <th>Contig</th>
        <th>Length</th>
        <th>Arrays</th>
        <th>Total TRA length</th>
        <th>Class mix</th>
        <th>Positions</th>
      </tr></thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
    </div>
    """


def render_trc_distribution(trc, seqid_lengths):
    """Return the `<h2> + summary + ideogram + minor-table` HTML block,
    or an empty string if the TRC has no arrays.

    The ideogram always draws the SAME set of major contigs for every
    TRC (the longest scaffolds in the assembly) even if this TRC has
    no arrays on them — the absence of signal on a chromosome is itself
    informative. The minor-contig table below only lists minor contigs
    that *do* carry arrays of this TRC, since listing the long tail of
    empty small contigs would drown out the actual data."""
    by_seqid = {}
    for arr in trc["arrays"]:
        by_seqid.setdefault(arr["seqid"], []).append(arr)
    if not by_seqid:
        return ""

    # Build a dict of all contigs we're aware of. The side-car gives
    # the assembly-wide set; any seqid that appears in arrays but is
    # missing from the side-car (legacy runs without seqid_lengths.tsv)
    # gets a lower-bound length from its max(end).
    contig_lengths = {sid: int(l) for sid, l in seqid_lengths.items()
                      if l and l > 0}
    for sid, arrs in by_seqid.items():
        if sid not in contig_lengths:
            contig_lengths[sid] = max(a["end"] for a in arrs)
    all_contigs = sorted(contig_lengths.items(),
                         key=lambda kv: (-kv[1], kv[0]))

    # Partition major/minor on the full assembly set (not just contigs
    # with arrays): for small assemblies every contig is major, for
    # fragmented ones the 1 Mb filter + 50-cap kicks in.
    if len(all_contigs) <= TRC_DIST_MAJOR_MAX_COUNT:
        majors = all_contigs
    else:
        majors = [c for c in all_contigs if c[1] >= TRC_DIST_MAJOR_MIN_LENGTH]
        if len(majors) > TRC_DIST_MAJOR_MAX_COUNT:
            majors = majors[:TRC_DIST_MAJOR_MAX_COUNT]
    major_ids = {s for s, _ in majors}
    # Minor table only lists minor contigs that actually carry arrays
    # of this TRC; empty minor contigs would just be noise.
    minors_with_arrays = [(s, l) for s, l in all_contigs
                          if s not in major_ids and s in by_seqid]

    n_majors_with_arrays = sum(1 for s, _ in majors if s in by_seqid)
    n_major_arrays       = sum(len(by_seqid[s]) for s, _ in majors if s in by_seqid)
    n_minor_arrays       = sum(len(by_seqid[s]) for s, _ in minors_with_arrays)

    summary_parts = []
    summary_parts.append(
        f"{len(majors)} major scaffold{'' if len(majors) == 1 else 's'} "
        f"shown ({n_majors_with_arrays} with arrays, "
        f"{n_major_arrays} array{'' if n_major_arrays == 1 else 's'})")
    if minors_with_arrays:
        summary_parts.append(
            f"{len(minors_with_arrays)} minor contig"
            f"{'' if len(minors_with_arrays) == 1 else 's'} with arrays "
            f"({n_minor_arrays} array{'' if n_minor_arrays == 1 else 's'})")
    summary = (
        f'<p style="font-size:12px;color:var(--fg-muted)">'
        f'Genome distribution of this TRC&#39;s arrays across the '
        f'assembly: '
        + " · ".join(summary_parts) +
        f'. Every major scaffold is drawn even when it carries no '
        f'arrays of this TRC &mdash; absence is informative too. '
        f'Array colour encodes the kitehor class &mdash; '
        f'{class_legend_inline()}.'
        f' Hover a rectangle to see coordinates and class.'
        f'</p>'
    )
    return (
        f'<h2>TRA genome distribution</h2>'
        f'{summary}'
        f'<div class="tc-ideogram">{_render_ideogram(majors, by_seqid)}</div>'
        f'{_render_minor_table(minors_with_arrays, by_seqid)}'
    )


# ----------------------------------------------------------------------
# Per-TRC dashboard
# ----------------------------------------------------------------------

def _is_ssr(trc):
    return (trc.get("repeat_type") or "").upper() == "SSR" or trc.get("ssr_motif")


_TIER_SHORT = {
    "HIGH":                 ("HIGH",   "high"),
    "LIKELY":               ("LIKELY", "likely"),
    "KMER_SUPPORT":         ("k-mer",  "kmer"),
    "AMBIGUOUS":            ("ambig",  "ambig"),
    "WEAK":                 ("weak",   "weak"),
    "OBSERVATIONAL":        ("obs",    "obs"),
    "REJECT_PHANTOM":       ("phantom",     "rej"),
    "REJECT_TOO_LARGE":     ("&gt;1/3",     "rej"),
    "REJECT_NO_FOUNDER":    ("no fdr",      "rej"),
    "REJECT_IS_FOUNDER":    ("is fdr",      "rej"),
    "REJECT_AMBIGUOUS_LOW": ("ambig low",   "rej"),
    "REJECT_BAD":           ("bad",         "rej"),
}


def _tier_pill(tier):
    """Coloured pill for a subrepeat tier value."""
    if not tier:
        return ""
    label, css = _TIER_SHORT.get(tier, ("?", "rej"))
    return f'<span class="tc-tier tc-tier-{css}">{label}</span>'


def _other_periods_cell(a):
    """Peaks with identity_med ≥ 0.7 excluding founder & strongest, plus
    short peaks (period < min_period so rescore returned NA identity)
    whose kite rank is better than the founder's rank. Sorted by
    identity desc (NAs by rank)."""
    peaks = a.get("rescored_peaks") or []
    if not peaks:
        return "—"
    founder = a.get("founder_period")
    strongest = a.get("strongest_period")
    founder_rank = None
    if founder is not None:
        for p in peaks:
            try:
                if int(round(float(p.get("period") or 0))) == founder:
                    founder_rank = int(p.get("rank") or 0)
                    break
            except (TypeError, ValueError):
                continue

    def _num(v):
        if v in (None, "", "NA"): return None
        try: return float(v)
        except (TypeError, ValueError): return None

    high_id = []     # (identity_med desc, period)
    short_high_rank = []  # (rank asc, period)
    for p in peaks:
        try:
            period = int(round(float(p.get("period") or 0)))
            rank = int(p.get("rank") or 0)
        except (TypeError, ValueError):
            continue
        if period in (founder, strongest):
            continue
        idm = _num(p.get("identity_med"))
        if idm is not None and idm >= 0.7:
            high_id.append((idm, period))
        elif idm is None and founder_rank is not None and rank < founder_rank:
            short_high_rank.append((rank, period))
    high_id.sort(key=lambda x: -x[0])
    short_high_rank.sort()

    def _id_str(i):
        return f".{int(round(i * 100)):02d}"
    bits = [f'{p}<span class="tc-score">({_id_str(i)})</span>'
            for i, p in high_id[:6]]
    bits += [f'{p}<span class="tc-score">(rank&nbsp;{r})</span>'
             for r, p in short_high_rank[:2]]
    return " ".join(bits) if bits else "—"


def _subrepeat_cell(a):
    """Up to two HIGH/LIKELY subrepeat candidates, period + occ + tier."""
    bits = []
    for n in (1, 2):
        period = a.get(f"subrepeat_{n}_period")
        tier   = a.get(f"subrepeat_{n}_tier")
        if not period or not tier:
            continue
        occ = a.get(f"subrepeat_{n}_occ")
        occ_str = f"occ&nbsp;{occ:.2f}" if occ is not None else ""
        bits.append(
            f'{esc(period)} {_tier_pill(tier)}'
            f'<span class="tc-meta">({occ_str})</span>')
    return "<br>".join(bits) if bits else "—"


def _ssr_cell(a):
    """Compact SSR summary: dominant motif + coverage% + optional top motifs."""
    motif = a.get("ssr_dominant_motif")
    if not motif or motif == "NA":
        return "—"
    pct = a.get("ssr_dominant_motif_coverage_pct")
    s = f'{esc(motif)} {pct:.0f}%' if pct is not None else esc(motif)
    top = a.get("ssr_top_motifs")
    if top and top not in (motif, f"{motif}:{pct}%"):
        # ssr_top_motifs uses "M:pct%;M2:pct%" — strip the dominant if echoed.
        parts = _fmt_motifs(top)
        if parts:
            s += f'<br><span class="tc-meta">top: {parts}</span>'
    return s


def _details_html(a):
    """HTML for the per-array Details child row (DataTables row.child()).
    Header line + SSR scan summary + table of every rescored peak,
    tier-labelled. Stored in `data-details` on the main row."""
    peaks = a.get("rescored_peaks") or []
    fp   = a.get("founder_period")
    sp   = a.get("strongest_period")
    mult = a.get("multiplicity") or 1
    delta = a.get("delta_id_pp")
    f_id = a.get("founder_id_med")
    s_id = a.get("strongest_id_med")
    fb   = a.get("founder_fallback")

    def _idstr(v):
        return f"id_med {v:.3f}" if v is not None else ""

    head_bits = []
    if fp is not None:
        head_bits.append(
            f'<strong>Founder:</strong> {esc(fp)}'
            f'{("<sup title=\"rescore fallback to top-scored kite peak\">*</sup>" if fb else "")}'
            f' <span class="tc-dim">{_idstr(f_id)}</span>')
    else:
        head_bits.append('<strong>Founder:</strong> NA')
    if sp is not None:
        head_bits.append(
            f'<strong>Strongest:</strong> {esc(sp)} '
            f'<span class="tc-dim">{_idstr(s_id)}</span>')
    irreg     = a.get("irregular_multiplicity")
    mult_raw  = a.get("multiplicity_raw")
    if irreg and mult_raw is not None:
        head_bits.append(
            f'<strong>&times;{mult}</strong> '
            f'<span class="tc-irreg" title="irregular multiplicity">irreg</span>'
            f' <span class="tc-dim">(k = {mult_raw:.3f})</span>')
    else:
        head_bits.append(f'<strong>&times;{mult}</strong>')
    if delta is not None:
        head_bits.append(f'<span class="tc-dim">&Delta;id {delta:+.2f}&nbsp;pp</span>')
    header = ('<div class="tc-details-summary">'
              + ' &nbsp;·&nbsp; '.join(head_bits)
              + '</div>')

    # SSR
    ssr_html = ""
    motif = a.get("ssr_dominant_motif")
    if motif and motif != "NA":
        pct = a.get("ssr_dominant_motif_coverage_pct")
        tot = a.get("ssr_total_coverage_pct")
        top = a.get("ssr_top_motifs")
        ssr_html = '<div class="tc-details-ssr">SSR scan: '
        ssr_html += f'<strong>{esc(motif)}</strong>'
        if pct is not None: ssr_html += f' {pct:.1f}% dominant'
        if tot is not None and tot != pct: ssr_html += f' &nbsp;·&nbsp; {tot:.1f}% total'
        if top and top != motif:
            ssr_html += f' &nbsp;·&nbsp; top: {_fmt_motifs(top)}'
        ssr_html += '</div>'

    # All peaks
    if not peaks:
        return header + ssr_html + (
            '<div class="tc-dim" style="margin-top:6px">'
            'no rescored peaks for this array</div>')

    fp_f = float(fp) if fp is not None else None
    def _f(v, n=3):
        if v in (None, "", "NA"): return ""
        try: return f"{float(v):.{n}f}"
        except (TypeError, ValueError): return esc(v)
    def _yn(v):
        s = (v or "").strip()
        return "✓" if s == "true" else ""
    body_rows = []
    for p in peaks:
        tier = classify_peak_tier(p, fp_f)
        body_rows.append(
            "<tr>"
            f"<td>{esc(p.get('rank',''))}</td>"
            f"<td>{esc(p.get('period',''))}</td>"
            f"<td>{_f(p.get('score'), 4)}</td>"
            f"<td>{_f(p.get('identity_med'))}</td>"
            f"<td>{_f(p.get('identity_iqr'))}</td>"
            f"<td>{_f(p.get('scan_occupancy_frac'))}</td>"
            f"<td>{esc(p.get('scan_n_intervals',''))}</td>"
            f"<td>{_f(p.get('coverage_frac'))}</td>"
            f"<td>{_f(p.get('spatial_contrast'))}</td>"
            f"<td>{_f(p.get('kmer_autocorr_founder'))}</td>"
            f"<td>{_f(p.get('kmer_phase_contrast'))}</td>"
            f"<td>{_yn(p.get('subrepeat'))}</td>"
            f"<td>{_yn(p.get('phantom'))}</td>"
            f"<td>{_tier_pill(tier)}</td>"
            "</tr>")
    table = (
        '<div class="tc-details-table-wrap">'
        '<table class="tc-details-table">'
        '<thead><tr>'
        '<th>rank</th><th>period</th><th>score</th>'
        '<th>id_med</th><th>id_iqr</th>'
        '<th>occ</th><th>scan_n</th><th>cov</th>'
        '<th>spat</th><th>autoF</th><th>phaseC</th>'
        '<th>sub</th><th>phtm</th><th>tier</th>'
        '</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody></table></div>')
    return header + ssr_html + table


def _arrays_table(arrays, include_hor=True):
    """Render the per-array DataTable body for a TRC dashboard.

    kitehor >=0.12.0 (include_hor=True) — columns:
       ▸  # · seqid · start · end · length · Founder · Δid · Strongest · ×k
          · Other periods · Subrepeat · SSR
    Each row carries a `data-details` attribute holding the HTML for the
    Details child row (founder/strongest summary, SSR scan, full table of
    rescored peaks with tier badges); a click on the ▸ control toggles
    it via DataTables' row().child() API (see tidecluster.js).
    The else branch keeps the legacy minimal layout for SSR-typed
    archival dirs without rescored peaks."""
    rows = []
    for i, a in enumerate(arrays, 1):
        if include_hor:
            fp     = a.get("founder_period")
            sp     = a.get("strongest_period")
            mult   = a.get("multiplicity") or 1
            delta  = a.get("delta_id_pp")
            fb     = a.get("founder_fallback")
            fp_html = (esc(fp) + ('<span class="tc-fallback">*</span>' if fb else "")
                       if fp is not None else "—")
            sp_html = esc(sp) if sp is not None else "—"
            delta_html = (f"{delta:+.1f}<span class='tc-meta'>&nbsp;pp</span>"
                          if delta is not None else "—")
            irreg     = a.get("irregular_multiplicity")
            mult_raw  = a.get("multiplicity_raw")
            mult_html = f"&times;{mult}"
            if irreg and mult_raw is not None:
                mult_html += (' <span class="tc-irreg" '
                              f'title="irregular multiplicity (k = {mult_raw:.2f})">~</span>')
            details_attr = html.escape(_details_html(a), quote=True)
            rows.append(
                f'<tr data-details="{details_attr}">'
                '<td class="tc-details-control"></td>'
                f'<td>{i}</td>'
                f'<td>{esc(a.get("seqid"))}</td>'
                f'<td data-order="{a.get("start") or 0}">{esc(a.get("start"))}</td>'
                f'<td data-order="{a.get("end") or 0}">{esc(a.get("end"))}</td>'
                f'<td data-order="{a.get("length") or 0}">{fmt_bp(a.get("length"))}</td>'
                f'<td data-order="{fp or 0}">{fp_html}</td>'
                f'<td data-order="{delta if delta is not None else -9999}">{delta_html}</td>'
                f'<td data-order="{sp or 0}">{sp_html}</td>'
                f'<td data-order="{mult_raw or mult}">{mult_html}</td>'
                f'<td>{_other_periods_cell(a)}</td>'
                f'<td>{_subrepeat_cell(a)}</td>'
                f'<td>{_ssr_cell(a)}</td>'
                '</tr>')
        else:
            rows.append(
                "<tr>"
                f'<td>{i}</td>'
                f'<td>{esc(a.get("seqid"))}</td>'
                f'<td data-order="{a.get("start") or 0}">{esc(a.get("start"))}</td>'
                f'<td data-order="{a.get("end") or 0}">{esc(a.get("end"))}</td>'
                f'<td data-order="{a.get("length") or 0}">{fmt_bp(a.get("length"))}</td>'
                f'<td>{esc(a.get("ssr") or "")}</td>'
                "</tr>")
    return "".join(rows)


def arrays_legend():
    """Trimmed per-TRA column legend. Wrapped in a <details> so it stays
    collapsed by default; click the summary to unfold. Structured as a
    proper <dl> (previous markup put <dt>/<dd> directly under a <div>,
    which DataTables sometimes mis-styled)."""
    return """
<details class="tc-legend">
  <summary>Column legend</summary>
  <dl class="tc-legend-list">
    <dt>Founder / Strongest</dt><dd>
        <strong>Strongest</strong> = highest-identity peak (rescore's
        <code>founder_period</code>). <strong>Founder</strong> = smallest
        divisor P of Strongest with k ∈ [2, 30] and id_med ≥ 0.7; falls
        back to Strongest. Red <code>*</code> = NA from rescore, fell
        back to top kite peak.</dd>
    <dt>&Delta;id</dt><dd>id(Strongest) − id(Founder), pp. NA when
        Founder = Strongest.</dd>
    <dt>&times;k</dt><dd>round(Strongest / Founder); 1 if none. Amber
        <code>~</code> pill = irregular (non-integer); hover for raw k.</dd>
    <dt>Other periods</dt><dd>peaks with id_med ≥ 0.7 other than
        Founder and Strongest, plus short peaks (period &lt;
        rescore's min-period) whose kite rank is better than the
        Founder's. Sorted by identity desc.</dd>
    <dt>Subrepeat</dt><dd>≤ 2 candidates ≤ Founder/3 with partial
        occupancy, ranked by <code>scan_occupancy_frac</code> desc.
        Tiered {hi} (scan + alignment / k-mer support agree) or
        {li} (scan only). Click
        <span class="tc-details-control tc-details-control-static"></span>
        on a row for every peak with its tier and full diagnostics.</dd>
    <dt>SSR</dt><dd>dominant short-motif tandem repeat from kitehor's
        ssr-scan + coverage %; top motifs in the expanded child row.</dd>
  </dl>
</details>
""".format(hi=_tier_pill("HIGH"), li=_tier_pill("LIKELY"))


def _variants_table(variants, tarean_dir, show=20, src_prefix="../"):
    """Top-N TAREAN variants table with graph + logo thumbs."""
    if not variants or not tarean_dir:
        return ""
    rows = []
    for v in variants[:show]:
        graph = (_thumb(f'{tarean_dir}/{v["graph_link"]}', "graph",
                        src_prefix=src_prefix, alt="k-mer graph")
                 if v.get("graph_link") else "")
        logo  = (_thumb(f'{tarean_dir}/{v["logo_link"]}',  "logo",
                        src_prefix=src_prefix, alt="logo")
                 if v.get("logo_link") else "")
        rows.append(
            "<tr>"
            f'<td>{esc(v["kmer"])}</td>'
            f'<td>{esc(v["variant"])}</td>'
            f'<td data-order="{v.get("total_score") or 0}">'
              f'{esc(round(v["total_score"], 4)) if v.get("total_score") is not None else ""}</td>'
            f'<td>{esc(v.get("consensus_length"))}</td>'
            f'<td>{esc(v.get("n_gap50"))}</td>'
            f'<td>{graph}</td>'
            f'<td>{logo}</td>'
            "</tr>"
        )
    more = ""
    if len(variants) > show:
        more = (f'<p class="tc-callout">Showing top {show} of '
                f'{len(variants)} variants by score. Full list in the original '
                f'<a href="{src_prefix}{esc(tarean_dir)}/report.html">TAREAN TRC report</a>.</p>')
    return f"""
    <details>
      <summary><strong>TAREAN variants tested</strong> ({len(variants)})</summary>
      <div class="tc-table-wrap">
      <table class="tc-table" data-order='[[2,"desc"]]'>
        <thead><tr>
          <th>k-mer length</th>
          <th>Variant index</th>
          <th>k-mer coverage score</th>
          <th>Consensus length</th>
          <th>n_gap50</th>
          <th>k-mer based graph</th>
          <th>Sequence logo</th>
        </tr></thead>
        <tbody>{"".join(rows)}</tbody>
      </table>
      </div>
      {more}
    </details>"""


def _trc_jumper(ordered_ids, current_id):
    opts = "".join(
        f'<option value="{esc(tid)}"{" selected" if tid == current_id else ""}>{esc(tid)}</option>'
        for tid in ordered_ids)
    return (f'<select class="tc-trc-jumper" data-base="">'
            f'{opts}</select>')


def render_trc_dashboard(trc, model, out_dir, ordered_ids, idx, run_meta, ctx):
    SITE = ctx["up"]           # prefix for sibling report pages
    SRC  = ctx["src_prefix"]   # prefix for original-output data files
    prev_id = ordered_ids[idx - 1] if idx > 0 else None
    next_id = ordered_ids[idx + 1] if idx + 1 < len(ordered_ids) else None
    prev_html = (f'<a class="tc-btn" href="{esc(prev_id)}.html">← {esc(prev_id)}</a>'
                 if prev_id else '<span class="tc-btn" style="opacity:.4">← prev</span>')
    next_html = (f'<a class="tc-btn" href="{esc(next_id)}.html">{esc(next_id)} →</a>'
                 if next_id else '<span class="tc-btn" style="opacity:.4">next →</span>')
    jump = _trc_jumper(ordered_ids, trc["id"])

    def trc_link(tid, text=None):
        return _trc_link(tid, text, prefix="")  # sibling in same dir
    def img_src(src_rel, width=None, alt=""):
        return _rel_img(src_rel, width=width, alt=alt, src_prefix=SRC)

    # Stats card
    stats_kv = [
        ("Type",         trc.get("repeat_type") or "—"),
        ("Arrays",       trc["n_arrays"]),
        ("Total size",   fmt_bp(trc["total_size"])),
        ("Array size",
         f'min {fmt_bp(trc.get("min_array"))} · max {fmt_bp(trc.get("max_array"))}'),
        ("Annotation",   trc.get("annotation") or "—"),
    ]
    stats_html = "".join(f'<dt>{esc(k)}</dt><dd>{esc(v)}</dd>' for k, v in stats_kv)

    # Classification card. For kitehor >=0.12.0 (combined_class is gone)
    # we surface the v0.12-shaped counters: HOR calls (multiplicity≥2),
    # subrepeats (HIGH/LIKELY in main column), SSR, fallback founders.
    # Old-dir rendering (with combined_class present) keeps the legacy
    # "Cluster class / Class breakdown" rows for archival reports.
    kite = trc.get("kite") or {}
    class_rows = []
    is_v012 = any(a.get("founder_period") is not None for a in trc["arrays"])
    if kite and not is_v012:
        cc_counts = kite.get("combined_class_counts") or {}
        dominant = kite.get("combined_class")
        if cc_counts or dominant:
            class_rows.append(("Cluster class",
                               class_badge(dominant) if dominant else "—"))
            class_rows.append(("Class breakdown",
                               class_mix_cell(Counter(cc_counts))))
            if cc_counts.get("hor") and kite.get("median_confidence") is not None:
                class_rows.append(("Median HOR confidence",
                                   f'{kite["median_confidence"]:.3f}'))
    if is_v012:
        n_hor    = sum(1 for a in trc["arrays"]
                       if (a.get("multiplicity") or 1) > 1)
        n_subrep = sum(1 for a in trc["arrays"] if a.get("subrepeat_1_tier"))
        n_ssr    = sum(1 for a in trc["arrays"]
                       if (a.get("ssr_dominant_motif") or "NA") != "NA")
        n_fb     = sum(1 for a in trc["arrays"] if a.get("founder_fallback"))
        # TRC-level aggregates from build_model: prevalent founder and
        # multiplicity distribution across the TRC's HOR arrays.
        prev_f   = kite.get("prevalent_founder")
        prev_n   = kite.get("prevalent_founder_n", 0)
        mc       = kite.get("multiplicity_counts") or {}
        n_irreg  = kite.get("n_irregular", 0)
        if prev_f:
            class_rows.append(
                ("Prevalent founder",
                 f'{prev_f} bp <span class="tc-dim">'
                 f'(in {prev_n}/{len(trc["arrays"])} arrays)</span>'))
        if mc:
            # Unicode ×/NBSP keep the line readable in the dl row even
            # though the class_html builder escapes plain values (only
            # values containing '<' bypass esc()).
            mc_str = ", ".join(f'×{k} ({v})'
                               for k, v in sorted(mc.items()))
            class_rows.append(("HOR multiplicities", mc_str))
        class_rows.append(("Arrays with HOR call (×k≥2)", str(n_hor)))
        if n_irreg:
            class_rows.append(
                ("Arrays with irregular multiplicity",
                 f'{n_irreg} <span class="tc-dim">(non-integer k, '
                 f'founder taken from TRC consensus)</span>'))
        class_rows.append(("Arrays with subrepeat",       str(n_subrep)))
        class_rows.append(("Arrays with SSR",             str(n_ssr)))
        if n_fb:
            class_rows.append(("Arrays with fallback founder", str(n_fb)))
    if trc.get("superfamily"):
        class_rows.append(
            ("Superfamily",
             f'<a href="{SITE}superfamilies.html#sf-{trc["superfamily"]}">SF {trc["superfamily"]}</a>'
             f' ({len(trc.get("sf_peers", []))} peers)'))
    else:
        class_rows.append(("Superfamily", "—"))
    class_html = "".join(
        f'<dt>{esc(k)}</dt><dd>{v if "<" in str(v) else esc(v)}</dd>'
        for k, v in class_rows)

    # Classification details (SSR motif / HOR / below-threshold callouts)
    callouts = []
    if _is_ssr(trc):
        detail = ("TAREAN is not performed on SSR TRCs by design; the "
                  "per-array table below carries the kitehor SSR classification "
                  "and motif coverage." if kite else
                  "TAREAN and KITE are not performed on SSR TRCs by design.")
        callouts.append(
            f'<div class="tc-callout"><strong>Simple Sequence Repeat.</strong> '
            f'{detail} '
            f'Motif(s): <code>{esc((trc.get("ssr_motif") or "").replace("%25","%"))}</code></div>')
    elif not trc.get("tarean"):
        callouts.append(
            '<div class="tc-callout"><strong>Below TAREAN threshold.</strong> '
            f'Total TRC length {fmt_bp(trc["total_size"])} is less than the '
            f'run&#39;s minimum total length; TAREAN and KITE were not '
            'performed for this cluster.</div>')

    # Consensus + graph/logo
    consensus_block = ""
    ta = trc.get("tarean") or {}
    if ta.get("consensus"):
        cons = ta["consensus"]
        pre_id = f"cons-{trc['id']}"
        length = len(cons)
        gc = 0
        if cons:
            gc = 100.0 * sum(1 for c in cons.upper() if c in "GC") / max(1, len(cons))
        graph = img_src(ta.get("graph_png"), width=160, alt="k-mer graph")
        logo  = img_src(ta.get("logo_png"),  width=320, alt="logo")
        consensus_block = f"""
        <h2>TAREAN consensus</h2>
        <div class="tc-consensus-toolbar">
          <span>monomer {ta.get("monomer_length")} bp · k-mer {ta.get("kmer")} ·
                score {ta.get("total_score"):.4f} ·
                length {length} bp · GC {gc:.1f}%</span>
          <button class="tc-btn" data-copy-target="#{pre_id}">copy</button>
          <a class="tc-btn" data-download-source="#{pre_id}"
             data-download="{esc(trc["id"])}_consensus.fasta">download .fa</a>
        </div>
        <pre id="{pre_id}" class="tc-consensus">&gt;{esc(trc["id"])}
{esc(cons)}</pre>
        <div style="display:flex; gap:16px; flex-wrap:wrap;">
          {f'<div class="tc-fig">{graph}<div class="tc-fig-caption">k-mer graph</div></div>' if graph else ""}
          {f'<div class="tc-fig">{logo}<div class="tc-fig-caption">sequence logo</div></div>' if logo else ""}
        </div>"""

    # Arrays table — use the kitehor Option-B layout whenever kitehor
    # produced per-array data, including SSR-typed TRCs (their arrays
    # carry the pure_ssr / tr_with_ssr structure cell). SSR TRCs without
    # kitehor data fall through to the minimal motif-only table.
    include_hor = any(a.get("m1") is not None for a in trc["arrays"])
    if include_hor:
        head = ('<th class="tc-details-col"></th>'
                "<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>Founder</th>"
                "<th>&Delta;id</th>"
                "<th>Strongest</th>"
                "<th>&times;k</th>"
                "<th>Other periods</th>"
                "<th>Subrepeat</th>"
                "<th>SSR</th>")
        legend = arrays_legend()
    elif _is_ssr(trc):
        head = ("<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>SSR motif</th>")
        legend = ""
    else:
        head = ("<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>—</th>")
        legend = ""
    # Start column moves one slot right when the ▸ details control is
    # prepended (include_hor=True) — keep ascending order on it.
    order_col = 3 if include_hor else 2
    arrays_section = f"""
    <h2>Tandem repeat arrays ({trc["n_arrays"]})</h2>
    {legend}
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[{order_col},"asc"]]'>
      <thead><tr>{head}</tr></thead>
      <tbody>{_arrays_table(trc["arrays"], include_hor=include_hor)}</tbody>
    </table>
    </div>"""

    # Genome distribution — every TRC gets this section regardless of
    # type (arrays are always present; HOR colouring gracefully falls
    # back to grey for arrays without a computed HOR status).
    distribution_section = render_trc_distribution(trc, model.get("seqid_lengths", {}))

    # KITE profile heatmap — suppressed for SSR (not biologically meaningful)
    kite_section = ""
    if kite and kite.get("profile_png") and not _is_ssr(trc):
        kite_section = f"""
        <h2>KITE array profile</h2>
        <div class="tc-fig">
          <img src="{SRC}{esc(kite["profile_png"])}"
               alt="per-array monomer-size heatmap">
          <div class="tc-fig-caption">per-array k-mer-interval monomer-size
          heatmap; rows are arrays, columns are candidate monomer lengths.</div>
        </div>"""

    # TAREAN variants drill-down — suppress for SSR
    variants_section = ""
    if ta.get("variants") and ta.get("tarean_dir") and not _is_ssr(trc):
        tarean_root = Path(model["paths"]["tarean_dir"] or "").as_posix()
        variants_section = _variants_table(
            ta["variants"], f'{tarean_root}/{ta["tarean_dir"]}',
            show=20, src_prefix=SRC)

    # Superfamily context
    sf_section = ""
    if trc.get("superfamily"):
        peers = ", ".join(trc_link(t) for t in trc.get("sf_peers", []))
        dotplot = ""
        if trc.get("sf_dotplot"):
            dotplot = (f'<div class="tc-fig"><img src="{SRC}{esc(trc["sf_dotplot"])}" '
                       f'width="500" alt="superfamily dotplot">'
                       f'<div class="tc-fig-caption">pairwise dotplot of '
                       f'superfamily {trc["superfamily"]} members</div></div>')
        sf_section = f"""
        <h2>Superfamily context</h2>
        <p>This TRC belongs to <a href="{SITE}superfamilies.html#sf-{trc["superfamily"]}">
        Superfamily {trc["superfamily"]}</a> with peers: {peers}.</p>
        {dotplot}"""

    # Parity links
    paths = model["paths"]
    tarean_href = (f'{SRC}{paths["tarean_dir"]}/{trc["id"]}.fasta_tarean/report.html'
                   if paths.get("tarean_dir") and ta else None)
    parity_bits = []
    if tarean_href:
        parity_bits.append(f'<a href="{esc(tarean_href)}">original TAREAN TRC page</a>')
    if paths.get("clustering_gff3"):
        parity_bits.append(f'<a href="{SRC}{esc(paths["clustering_gff3"])}">clustering.gff3</a>')
    if paths.get("annotation_gff3"):
        parity_bits.append(f'<a href="{SRC}{esc(paths["annotation_gff3"])}">annotation.gff3</a>')
    parity_html = (f'<p class="tc-callout">Parity links: {" · ".join(parity_bits)}</p>'
                   if parity_bits else "")

    body = f"""
    <div class="tc-trc-nav">
      {prev_html}
      {jump}
      {next_html}
      <span class="tc-spacer" style="flex:1"></span>
    </div>
    <h1>{esc(trc["id"])}</h1>
    {"".join(callouts)}
    <section class="tc-cards">
      <div class="tc-card"><div class="tc-card-title">Cluster stats</div>
        <dl class="tc-kv">{stats_html}</dl></div>
      <div class="tc-card"><div class="tc-card-title">Classification</div>
        <dl class="tc-kv">{class_html}</dl></div>
    </section>
    {consensus_block}
    {arrays_section}
    {distribution_section}
    {kite_section}
    {variants_section}
    {sf_section}
    {parity_html}
    """
    (out_dir / f"{trc['id']}.html").write_text(
        _shell(ctx, trc["id"], "trcs", run_meta, body))


def render_all_trc_dashboards(model, out_dir, run_meta, ctx):
    out_dir.mkdir(parents=True, exist_ok=True)
    ordered = [t["id"] for t in model["trcs"]]  # already sorted numerically
    for i, trc in enumerate(model["trcs"]):
        render_trc_dashboard(trc, model, out_dir, ordered, i, run_meta, ctx)
    return len(ordered)


# ----------------------------------------------------------------------
# Asset copy
# ----------------------------------------------------------------------

def copy_assets(dest: Path, script_path: Path):
    src = script_path.parent / "tarean" / "assets"
    if not src.exists():
        raise SystemExit(f"asset source missing: {src}")
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def build_report(input_dir, prefix=None, output_dir=None, *, quiet=False):
    """Importable entry point for the report-v2 build.

    Layout (since v1.9.0):
      <input_dir>/<prefix>_index.html         v2 landing (Summary)
      <input_dir>/<prefix>_report/            v2 subpages + assets
      <input_dir>/<prefix>_report_legacy/     v1 HTML moved here by
                                              TideCluster.py._move_v1_to_legacy

    Returns the subtree output directory as a Path."""
    input_dir = Path(input_dir).resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"input-dir is not a directory: {input_dir}")
    prefix  = prefix or detect_prefix(input_dir)

    # One-off cleanup: earlier 1.9.0-dev snapshots wrote to
    # <prefix>_report_v2/. If a stale one is present, remove it so it
    # doesn't clutter the output root next to the new <prefix>_report/.
    stale = input_dir / f"{prefix}_report_v2"
    if stale.exists():
        shutil.rmtree(stale)
        if not quiet:
            print(f"removed stale {stale.name}/ from the previous 1.9.0-dev layout")

    out_dir = Path(output_dir) if output_dir else (input_dir / f"{prefix}_report")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "data").mkdir(exist_ok=True)

    if not quiet:
        print(f"input:   {input_dir}")
        print(f"prefix:  {prefix}")
        print(f"root:    {input_dir}/{prefix}_index.html")
        print(f"subtree: {out_dir}")

    model = build_model(input_dir, prefix)

    trc_count = len(model["trcs"])
    analysed  = sum(1 for t in model["trcs"] if t.get("tarean"))
    sf_count  = len(model["superfamilies"])
    if not quiet:
        if model.get("is_v012"):
            n_hor    = sum(1 for t in model["trcs"] for a in t["arrays"]
                           if (a.get("multiplicity") or 1) > 1)
            n_subrep = sum(1 for t in model["trcs"] for a in t["arrays"]
                           if a.get("subrepeat_1_tier"))
            n_ssr    = sum(1 for t in model["trcs"] for a in t["arrays"]
                           if (a.get("ssr_dominant_motif") or "NA") != "NA")
            n_fb     = sum(1 for t in model["trcs"] for a in t["arrays"]
                           if a.get("founder_fallback"))
            signal_str = (f"HOR(×k≥2): {n_hor}, subrep: {n_subrep}, "
                          f"SSR: {n_ssr}, fallback: {n_fb}")
        else:
            cls_total = class_totals(model)
            signal_str = ", ".join(f"{class_label(raw)}: {cls_total[raw]}"
                                   for raw in CANONICAL_CLASS_ORDER
                                   if cls_total.get(raw)) or "n/a"
        print(f"TRCs: {trc_count}  TAREAN-analysed: {analysed}  "
              f"superfamilies: {sf_count}  signals: {signal_str}")

    copy_assets(out_dir / "assets", Path(__file__).resolve())
    copy_report_content_source(out_dir)
    with (out_dir / "data" / "report.json").open("w") as f:
        json.dump(model, f, indent=2, default=str)
    if not quiet:
        print(f"wrote: {out_dir / 'data' / 'report.json'}")

    run_meta = make_run_meta(model)
    # Three depth contexts share a schema. See each field in page_shell.
    # *_href fields are full relative hrefs from the page being rendered;
    # *_prefix fields are prepended to local refs in the page body.
    ctx_root = {                  # <input_dir>/<prefix>_index.html
        "up":              f"{prefix}_report/",
        "root_href":       f"{prefix}_index.html",
        "legacy_href":     f"{prefix}_report_legacy/{prefix}_index.html",
        "assets_href":     f"{prefix}_report/assets/",
        "src_prefix":      "",
        "trc_link_prefix": f"{prefix}_report/trc/",
    }
    ctx_top = {                   # <out_dir>/{tarean,kite,superfamilies}.html
        "up":              "",
        "root_href":       f"../{prefix}_index.html",
        "legacy_href":     f"../{prefix}_report_legacy/{prefix}_index.html",
        "assets_href":     "assets/",
        "src_prefix":      "../",
        "trc_link_prefix": "trc/",
    }
    ctx_dash = {                  # <out_dir>/trc/TRC_N.html
        "up":              "../",
        "root_href":       f"../../{prefix}_index.html",
        "legacy_href":     f"../../{prefix}_report_legacy/{prefix}_index.html",
        "assets_href":     "../assets/",
        "src_prefix":      "../../",
        "trc_link_prefix": "",
    }

    render_index(model, input_dir / f"{prefix}_index.html", run_meta, ctx_root)
    render_tarean(model, out_dir / "tarean.html", run_meta, ctx_top)
    render_kite(model, out_dir / "kite.html", run_meta, ctx_top)
    render_superfamilies(model, out_dir / "superfamilies.html", run_meta, ctx_top)
    n_dash = render_all_trc_dashboards(model, out_dir / "trc", run_meta, ctx_dash)
    if not quiet:
        print(f"rendered: {prefix}_index.html + tarean / kite / "
              f"superfamilies + {n_dash} TRC dashboards")
    return out_dir


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Directory containing <prefix>_*.gff3 etc.")
    p.add_argument("--prefix",
                   help="Prefix used when the pipeline ran (auto-detect if only one).")
    p.add_argument("--output-dir", type=Path,
                   help="Output directory; defaults to <input-dir>/<prefix>_report_v2/.")
    args = p.parse_args(argv)
    build_report(args.input_dir, prefix=args.prefix, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
