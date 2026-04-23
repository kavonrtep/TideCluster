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
<body>
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

    Both schemas are accepted:
      - New (post-1.9.0): reads HOR_status / HOR_confidence / HOR_base_monomer
        / HOR_hor_period / HOR_n_harmonics directly.
      - Legacy (pre-1.9.0, incl. HOR_multiple column): the HOR columns are
        recomputed from (m1..m3, s1..s3) using compute_hor(). This lets
        rerender produce up-to-date HOR calls on old output directories
        without re-running the KITE R script."""
    p = paths["kite_top3_csv"]
    if not p:
        return {}
    out = {}
    # Peek at headers to decide legacy vs new.
    with open(p, newline="") as f:
        header = next(csv.reader(f, delimiter="\t"))
    is_new = "HOR_confidence" in header
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

        if is_new:
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
            "array_length": ni(row.get("array_length")),
            **hor,
        }
    return out


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
        for arr in arrays:
            key = (trc_id, arr["seqid"], arr["start"], arr["end"])
            k = kite_by_array.get(key)
            if k:
                arr.update(k)
                status = k.get("hor_status")
                if   status == "HOR strong":    kite_counts["hor_strong"]   += 1
                elif status == "HOR moderate":  kite_counts["hor_moderate"] += 1
                elif status == "HOR weak":      kite_counts["hor_weak"]     += 1
                elif status == "No HOR":        kite_counts["no_hor"]       += 1
                c = k.get("hor_confidence")
                if c is not None: confidences.append(c)
                if k.get("m1") is not None:
                    m1_values.append(k["m1"])
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

    return {
        "schema_version": __version__,
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
    rows_stats = [
        ("Number of TRCs",                   stats.get("n_trcs_total")),
        ("TRCs above threshold",             stats.get("n_trcs_above_threshold")),
        ("Number of SSRs in TRCs",           stats.get("n_ssrs")),
        ("Total length of TRCs",             fmt_bp(stats.get("total_tr_length"))),
        ("Number of TRAs",                   stats.get("n_tras")),
        ("Input sequence length",            fmt_bp(stats.get("input_sequence_length"))),
    ]
    stats_html = "".join(
        f'<dt>{esc(label)}</dt><dd>{esc(val) if val != "" else ""}</dd>'
        for label, val in rows_stats)
    hor_total = {"strong": 0, "mod": 0, "weak": 0, "none": 0}
    for t in model["trcs"]:
        if t["kite"]:
            hor_total["strong"] += t["kite"]["n_hor_strong"]
            hor_total["mod"]    += t["kite"]["n_hor_moderate"]
            hor_total["weak"]   += t["kite"]["n_hor_weak"]
            hor_total["none"]   += t["kite"]["n_no_hor"]
    n_tarean    = sum(1 for t in model["trcs"] if t["tarean"])
    n_sf_peers  = sum(1 for t in model["trcs"] if t["superfamily"])
    up = ctx["up"]
    cards = f"""
    <section class="tc-cards">
      <div class="tc-card"><div class="tc-card-title">Run summary</div>
        <dl class="tc-kv">{stats_html}</dl></div>
      <div class="tc-card"><div class="tc-card-title">Reports</div>
        <dl class="tc-kv">
          <dt>TRCs total</dt>     <dd><a href="{up}tarean.html">{len(model["trcs"])}</a></dd>
          <dt>TAREAN-analysed</dt><dd><a href="{up}tarean.html">{n_tarean}</a></dd>
          <dt>HOR calls</dt>      <dd>{hor_badge("HOR strong")} {hor_total["strong"]} ·
                                      {hor_badge("HOR moderate")} {hor_total["mod"]} ·
                                      {hor_badge("HOR weak")} {hor_total["weak"]}</dd>
          <dt>Superfamilies</dt>  <dd><a href="{up}superfamilies.html">{len(model["superfamilies"])}</a>,
                                      {n_sf_peers} TRCs grouped</dd>
        </dl></div>
      <div class="tc-card"><div class="tc-card-title">Run settings</div>
        <dl class="tc-kv">{"".join(rows_settings)}</dl></div>
    </section>"""
    main = f"""
    <h1>TideCluster report — {esc(model["meta"]["prefix"])}</h1>
    {cards}
    <p class="tc-callout">This is report v2 — the default TideCluster
    report since 1.9.0. The original v1 report is preserved under
    <code>{esc(model["meta"]["prefix"])}_report_legacy/</code>; reach it
    via the <strong>Legacy report&nbsp;&#x2197;</strong> link at the
    top-right of the navigation bar.</p>
    """
    Path(out_path).write_text(_shell(ctx, "Summary", "index", run_meta, main))


def _trc_type_label(t):
    """Human-readable type label for the merged TAREAN table."""
    if (t.get("repeat_type") or "").upper() == "SSR" or t.get("ssr_motif"):
        return "SSR"
    if t.get("tarean"):
        return "TR"
    return "TR&nbsp;(below&nbsp;threshold)"


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
    hor_cell = (f'{hor_count_cell(kite["n_hor_strong"],   "strong")}'
                f'{hor_count_cell(kite["n_hor_moderate"], "mod")}'
                f'{hor_count_cell(kite["n_hor_weak"],     "weak")}'
                f'{hor_count_cell(kite["n_no_hor"],       "none")}') if kite else ""
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
        f'<td>{hor_cell}</td>'
        f'<td>{sf_cell}</td>'
        f'<td>{graph_thumb}</td>'
        f'<td>{logo_thumb}</td>'
        f'<td><code style="font-size:11px;">{cons_short}</code></td>'
        "</tr>")


def render_tarean(model, out_path, run_meta, ctx):
    """Merged TAREAN / All-TRCs table. Every TRC from clustering.gff3
    appears here (TAREAN-analysed, below-threshold, SSR). TAREAN-specific
    columns are empty for rows where TAREAN was not performed."""
    rows = "".join(_render_tarean_row(t, ctx) for t in model["trcs"])
    n_total   = len(model["trcs"])
    n_tarean  = sum(1 for t in model["trcs"] if t["tarean"])
    n_ssr     = sum(1 for t in model["trcs"] if _trc_type_label(t) == "SSR")
    n_below   = n_total - n_tarean - n_ssr
    body = f"""
    <h1>Tandem repeat clusters (TAREAN)</h1>
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
        <th>HOR<br>(strong / mod / weak / none)</th>
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
    hor_total = {"strong": 0, "mod": 0, "weak": 0, "none": 0}
    for t in model["trcs"]:
        if t["kite"]:
            hor_total["strong"] += t["kite"]["n_hor_strong"]
            hor_total["mod"]    += t["kite"]["n_hor_moderate"]
            hor_total["weak"]   += t["kite"]["n_hor_weak"]
            hor_total["none"]   += t["kite"]["n_no_hor"]
    total = sum(hor_total.values()) or 1
    bar = ""
    for kind, label in (("strong", "HOR strong"), ("mod", "HOR moderate"),
                        ("weak",   "HOR weak"),   ("none", "No HOR")):
        pct = 100 * hor_total[kind] / total
        if pct < 1 and hor_total[kind] == 0:
            continue
        bar += (f'<div class="tc-bar-seg hor-{kind}" '
                f'style="width:{pct:.2f}%" title="{label}: {hor_total[kind]}">'
                f'{hor_total[kind] if pct > 3 else ""}</div>')
    rows = []
    for t in model["trcs"]:
        if not t["kite"]:
            continue
        k = t["kite"]
        profile_thumb = _thumb(k.get("profile_top3_png"), "profile",
                               src_prefix=ctx["src_prefix"], alt="KITE profile")
        tarean_mon = (t["tarean"] or {}).get("monomer_length")
        kite_mon   = k.get("monomer_primary")
        med_conf   = k.get("median_confidence")
        med_str    = f"{med_conf:.3f}" if med_conf is not None else ""
        rows.append(
            "<tr>"
            f'<td data-order="{t["index"]}">{_trc_link(t["id"], prefix=ctx["trc_link_prefix"])}</td>'
            f'<td data-order="{kite_mon or 0}">{esc(kite_mon)}</td>'
            f'<td data-order="{tarean_mon or 0}">{esc(tarean_mon)}</td>'
            f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
            f'<td data-order="{t["total_size"]}">{fmt_bp(t["total_size"])}</td>'
            f'<td data-order="{k["n_hor_strong"]}">{hor_count_cell(k["n_hor_strong"],   "strong")}</td>'
            f'<td data-order="{k["n_hor_moderate"]}">{hor_count_cell(k["n_hor_moderate"], "mod")}</td>'
            f'<td data-order="{k["n_hor_weak"]}">{hor_count_cell(k["n_hor_weak"],     "weak")}</td>'
            f'<td data-order="{k["n_no_hor"]}">{hor_count_cell(k["n_no_hor"],       "none")}</td>'
            f'<td data-order="{med_conf or 0}">{med_str}</td>'
            f'<td>{profile_thumb}</td>'
            "</tr>")
    body = f"""
    <h1>K-mer interval tandem repeat estimation (KITE)</h1>
    <h3>HOR distribution across all analysed arrays</h3>
    <div class="tc-bar" title="Total arrays: {total}">{bar}</div>
    <p style="font-size:12px;color:var(--fg-muted)">
      {hor_badge("HOR strong")} {hor_total["strong"]}
      · {hor_badge("HOR moderate")} {hor_total["mod"]}
      · {hor_badge("HOR weak")} {hor_total["weak"]}
      · {hor_badge("No HOR")} {hor_total["none"]}
      (total {total} arrays). See <code>docs/hor_classification.md</code>
      for the scoring formula and category thresholds.
    </p>
    <h2>Per-TRC HOR summary</h2>
    <p style="font-size:12px;color:var(--fg-muted)">
      <strong>Monomer (KITE)</strong> is the most-frequent primary
      k-mer-interval estimate across the arrays of a given TRC.
      <strong>Monomer (TAREAN)</strong> is the TAREAN de&nbsp;Bruijn
      consensus length. <strong>Median conf.</strong> is the median
      HOR confidence of the TRC's arrays.
    </p>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[0,"asc"]]'>
      <thead><tr>
        <th>TRC</th>
        <th>Monomer<br>(KITE&nbsp;primary)</th>
        <th>Monomer<br>(TAREAN)</th>
        <th>Arrays</th>
        <th>Total size</th>
        <th>HOR<br>strong</th>
        <th>HOR<br>moderate</th>
        <th>HOR<br>weak</th>
        <th>No&nbsp;HOR</th>
        <th>Median<br>conf.</th>
        <th>Profile</th>
      </tr></thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
    </div>
    """
    Path(out_path).write_text(_shell(ctx, "KITE", "kite", run_meta, body))


def render_superfamilies(model, out_path, run_meta, ctx):
    out_path = Path(out_path)
    if not model["superfamilies"]:
        body = ('<h1>TRC Superfamilies</h1>'
                '<p class="tc-callout">No TRC superfamilies were identified '
                '(either insufficient TRCs reached the clustering step or none shared '
                'significant consensus similarity).</p>')
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
_HOR_FILL = {
    "HOR strong":   "#2e8b2e",  # forest green
    "HOR moderate": "#e88f00",  # saturated amber
    "HOR weak":     "#d4a017",  # saturated gold
    "No HOR":       "#909090",  # medium grey
}


def _hor_fill(arr):
    return _HOR_FILL.get(arr.get("hor_status"), "#cccccc")


def _array_title(arr):
    conf = arr.get("hor_confidence")
    status = arr.get("hor_status") or "—"
    span = f"{arr['seqid']}:{arr['start']:,}-{arr['end']:,}"
    length = f"{fmt_bp(arr['end'] - arr['start'])}"
    conf_s = f"{conf:.3f}" if conf is not None else "—"
    return f"{span} · {length} · {status} · conf {conf_s}"


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
                f'fill="{_hor_fill(arr)}" stroke="none" '
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
        statuses = Counter(a.get("hor_status", "No HOR") for a in arrs)
        mix_cells = (
            f'{hor_count_cell(statuses.get("HOR strong", 0),   "strong")}'
            f'{hor_count_cell(statuses.get("HOR moderate", 0), "mod")}'
            f'{hor_count_cell(statuses.get("HOR weak", 0),     "weak")}'
            f'{hor_count_cell(statuses.get("No HOR", 0),       "none")}'
        )
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
                f'fill="{_hor_fill(arr)}" stroke="none" '
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
        <th>HOR mix<br>(strong / mod / weak / none)</th>
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
        f'Array colour encodes HOR confidence &mdash; '
        f'{hor_badge("HOR strong")} '
        f'{hor_badge("HOR moderate")} '
        f'{hor_badge("HOR weak")} '
        f'{hor_badge("No HOR")}.'
        f' Hover a rectangle to see coordinates and confidence.'
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


def _fmt_score(x):
    if x is None: return ""
    try: return f"{float(x):.3f}"
    except (TypeError, ValueError): return esc(x)


def _arrays_table(arrays, include_hor=True):
    """Render the per-array DataTable body for a TRC dashboard.

    When include_hor is True the table carries m1..m3, s1..s3,
    HOR status, confidence, base monomer, and HOR period columns.
    SSR dashboards instead get an `SSR motif` column since neither
    TAREAN nor KITE apply to them meaningfully."""
    rows = []
    for i, a in enumerate(arrays, 1):
        if include_hor:
            hor_bcell = hor_badge(a.get("hor_status"))
            conf = a.get("hor_confidence")
            base = a.get("hor_base_monomer")
            perd = a.get("hor_hor_period")
            conf_str = f"{conf:.3f}" if conf is not None else ""
            rows.append(
                "<tr>"
                f'<td>{i}</td>'
                f'<td>{esc(a.get("seqid"))}</td>'
                f'<td data-order="{a.get("start") or 0}">{esc(a.get("start"))}</td>'
                f'<td data-order="{a.get("end") or 0}">{esc(a.get("end"))}</td>'
                f'<td data-order="{a.get("length") or 0}">{fmt_bp(a.get("length"))}</td>'
                f'<td data-order="{a.get("m1") or 0}">{esc(a.get("m1"))}</td>'
                f'<td data-order="{a.get("s1") or 0}">{_fmt_score(a.get("s1"))}</td>'
                f'<td data-order="{a.get("m2") or 0}">{esc(a.get("m2"))}</td>'
                f'<td data-order="{a.get("s2") or 0}">{_fmt_score(a.get("s2"))}</td>'
                f'<td data-order="{a.get("m3") or 0}">{esc(a.get("m3"))}</td>'
                f'<td data-order="{a.get("s3") or 0}">{_fmt_score(a.get("s3"))}</td>'
                f'<td>{hor_bcell}</td>'
                f'<td data-order="{conf or 0}">{conf_str}</td>'
                f'<td data-order="{base or 0}">{esc(base) if base else ""}</td>'
                f'<td data-order="{perd or 0}">{esc(perd) if perd else ""}</td>'
                "</tr>")
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


ARRAYS_LEGEND = """
<div class="tc-legend">
  <dt>m<sub>1</sub></dt><dd>Monomer size — primary estimate</dd>
  <dt>m<sub>2</sub>, m<sub>3</sub></dt><dd>Alternative estimates (2nd and 3rd KITE peaks)</dd>
  <dt>s<sub>1</sub>, s<sub>2</sub>, s<sub>3</sub></dt><dd>k-mer-interval scores for each peak</dd>
  <dt>HOR status</dt><dd>category derived from the continuous confidence:
      <span class="hor-badge hor-none">No HOR</span>
      <span class="hor-badge hor-weak">HOR weak</span>
      <span class="hor-badge hor-mod">HOR moderate</span>
      <span class="hor-badge hor-strong">HOR strong</span></dd>
  <dt>Confidence</dt><dd>continuous HOR score (base&times;harmonic geometric-mean,
      bonus for multiple harmonics); 0 when the array has no base+harmonic structure.</dd>
  <dt>Base (bp) / HOR period (bp)</dt><dd>fitted base monomer <em>m*</em> and
      the largest supported harmonic <em>k<sub>max</sub>&middot;m*</em>.</dd>
</div>
"""


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

    # Classification / HOR card
    kite = trc.get("kite") or {}
    class_rows = []
    if kite:
        class_rows.append(("HOR strong arrays",   kite["n_hor_strong"]))
        class_rows.append(("HOR moderate arrays", kite["n_hor_moderate"]))
        class_rows.append(("HOR weak arrays",     kite["n_hor_weak"]))
        class_rows.append(("No HOR detected",     kite["n_no_hor"]))
        if kite.get("median_confidence") is not None:
            class_rows.append(("Median HOR confidence",
                               f'{kite["median_confidence"]:.3f}'))
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
        callouts.append(
            f'<div class="tc-callout"><strong>Simple Sequence Repeat.</strong> '
            f'TAREAN and KITE are not performed on SSR TRCs by design. '
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

    # Arrays table
    include_hor = not _is_ssr(trc) and any(a.get("m1") is not None for a in trc["arrays"])
    if include_hor:
        head = ("<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>m<sub>1</sub></th><th>s<sub>1</sub></th>"
                "<th>m<sub>2</sub></th><th>s<sub>2</sub></th>"
                "<th>m<sub>3</sub></th><th>s<sub>3</sub></th>"
                "<th>HOR status</th>"
                "<th>Confidence</th>"
                "<th>Base (bp)</th>"
                "<th>HOR period (bp)</th>")
        legend = ARRAYS_LEGEND
    elif _is_ssr(trc):
        head = ("<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>SSR motif</th>")
        legend = ""
    else:
        head = ("<th>#</th><th>seqid</th><th>start</th><th>end</th><th>length</th>"
                "<th>—</th>")
        legend = ""
    arrays_section = f"""
    <h2>Tandem repeat arrays ({trc["n_arrays"]})</h2>
    {legend}
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[2,"asc"]]'>
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
    hor_s = sum(t["kite"]["n_hor_strong"]   for t in model["trcs"] if t.get("kite"))
    hor_m = sum(t["kite"]["n_hor_moderate"] for t in model["trcs"] if t.get("kite"))
    hor_w = sum(t["kite"]["n_hor_weak"]     for t in model["trcs"] if t.get("kite"))
    if not quiet:
        print(f"TRCs: {trc_count}  TAREAN-analysed: {analysed}  "
              f"superfamilies: {sf_count}  "
              f"HOR strong/moderate/weak: {hor_s}/{hor_m}/{hor_w}")

    copy_assets(out_dir / "assets", Path(__file__).resolve())
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
