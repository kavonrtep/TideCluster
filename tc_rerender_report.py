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
import sys
from pathlib import Path

__version__ = "1"  # schema version for report.json

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
    cls = {"HOR-dominant": "hor-dom",
           "HOR-visible":  "hor-vis",
           "No HOR detected": "hor-none"}.get(status, "hor-none")
    return f'<span class="hor-badge {cls}">{esc(status)}</span>'


def hor_count_cell(n, kind):
    cls = {"dom": "hor-cnt-dom", "vis": "hor-cnt-vis", "none": "hor-cnt-none"}[kind]
    return f'<span class="{cls}">{int(n)}</span>'


def page_shell(title, active, run_meta, main_html, extra_head=""):
    """Return a full HTML document string wrapping `main_html`."""
    nav_items = [
        ("index.html",          "Summary",       "index"),
        ("trcs.html",           "All TRCs",      "trcs"),
        ("tarean.html",         "TAREAN",        "tarean"),
        ("kite.html",           "KITE",          "kite"),
        ("superfamilies.html",  "Superfamilies", "sf"),
    ]
    def _nav_link(href, label, key):
        cls = ' class="active"' if key == active else ""
        return f'<a href="{href}"{cls}>{label}</a>'
    nav_html = "".join(_nav_link(*item) for item in nav_items)
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
<link rel="stylesheet" href="assets/datatables/datatables.min.css">
<link rel="stylesheet" href="assets/tidecluster.css">
{extra_head}
</head>
<body>
<nav class="tc-nav">
  <span class="tc-brand">TideCluster</span>
  {nav_html}
  <span class="tc-spacer"></span>
</nav>
<div class="tc-run-strip">{run_strip}</div>
<main>
{main_html}
</main>
<script src="assets/datatables/jquery-3.7.1.min.js"></script>
<script src="assets/datatables/datatables.min.js"></script>
<script src="assets/tidecluster.js"></script>
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
    """Per-TRA KITE records keyed by (trc_id, seqid, start, end)."""
    p = paths["kite_top3_csv"]
    if not p:
        return {}
    out = {}
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
        hor_mult = row.get("HOR_multiple")
        out[(trc, sid, start, end)] = {
            "m1": ni(row.get("monomer_size")),
            "s1": num(row.get("score")),
            "m2": ni(row.get("monomer_size_2")),
            "s2": num(row.get("score_2")),
            "m3": ni(row.get("monomer_size_3")),
            "s3": num(row.get("score_3")),
            "array_length": ni(row.get("array_length")),
            "hor_status":   (row.get("HOR_status") or "").strip() or None,
            "hor_multiple": None if hor_mult in (None, "", "NA") else int(float(hor_mult)),
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
        kite_counts = {"no_hor": 0, "hor_visible": 0, "hor_dominant": 0}
        for arr in arrays:
            key = (trc_id, arr["seqid"], arr["start"], arr["end"])
            k = kite_by_array.get(key)
            if k:
                arr.update(k)
                status = k.get("hor_status")
                if status == "HOR-dominant": kite_counts["hor_dominant"] += 1
                elif status == "HOR-visible": kite_counts["hor_visible"] += 1
                elif status == "No HOR detected": kite_counts["no_hor"] += 1
        # TRC-level aggregates
        repeat_type = next((a["repeat_type"] for a in arrays if a.get("repeat_type")), None)
        annotation  = next((a["annotation"]  for a in arrays if a.get("annotation")),  None)
        ssr_motif   = next((a["ssr"]         for a in arrays if a.get("ssr")),         None)
        lengths = [a["length"] for a in arrays]
        tarean = tarean_summary.get(trc_id)
        has_tarean = tarean is not None and (paths["tarean_dir"]
                                             and (paths["tarean_dir"] / f"{trc_id}.fasta_tarean").exists())
        has_kite = any(k for k in [arrays] if any("m1" in a for a in arrays))
        kite_block = None
        if has_kite and kite_dir_rel:
            kite_block = {
                "n_no_hor":       kite_counts["no_hor"],
                "n_hor_visible":  kite_counts["hor_visible"],
                "n_hor_dominant": kite_counts["hor_dominant"],
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
        "trcs":          trcs,
        "superfamilies": superfams,
    }


# ----------------------------------------------------------------------
# HTML page renderers
# ----------------------------------------------------------------------

def _trc_link(trc_id, text=None):
    return f'<a href="trc/{esc(trc_id)}.html">{esc(text or trc_id)}</a>'


def _rel_img(src_rel, width=None, alt=""):
    """An <img> whose src points back into the original output dir."""
    if not src_rel:
        return ""
    w = f' width="{int(width)}"' if width else ""
    return f'<img src="../{esc(src_rel)}" alt="{esc(alt)}"{w}>'


def _tarean_asset_path(trc):
    """Relative path from report_v2/ to the per-TRC TAREAN image directory."""
    tdir = (trc.get("tarean") or {}).get("tarean_dir")
    return f"../{model_tarean_root}/{tdir}" if tdir else None  # not used; explicit per site


def render_index(model, out_dir, run_meta):
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
    hor_total = {"dom": 0, "vis": 0, "none": 0}
    for t in model["trcs"]:
        if t["kite"]:
            hor_total["dom"]  += t["kite"]["n_hor_dominant"]
            hor_total["vis"]  += t["kite"]["n_hor_visible"]
            hor_total["none"] += t["kite"]["n_no_hor"]
    n_tarean    = sum(1 for t in model["trcs"] if t["tarean"])
    n_sf_peers  = sum(1 for t in model["trcs"] if t["superfamily"])
    cards = f"""
    <section class="tc-cards">
      <div class="tc-card"><div class="tc-card-title">Run summary</div>
        <dl class="tc-kv">{stats_html}</dl></div>
      <div class="tc-card"><div class="tc-card-title">Reports</div>
        <dl class="tc-kv">
          <dt>All TRCs</dt>       <dd><a href="trcs.html">{len(model["trcs"])} total</a></dd>
          <dt>TAREAN-analysed</dt><dd><a href="tarean.html">{n_tarean}</a></dd>
          <dt>HOR calls</dt>      <dd>{hor_badge("HOR-dominant")} {hor_total["dom"]} ·
                                      {hor_badge("HOR-visible")} {hor_total["vis"]}</dd>
          <dt>Superfamilies</dt>  <dd><a href="superfamilies.html">{len(model["superfamilies"])}</a>,
                                      {n_sf_peers} TRCs grouped</dd>
        </dl></div>
      <div class="tc-card"><div class="tc-card-title">Run settings</div>
        <dl class="tc-kv">{"".join(rows_settings)}</dl></div>
    </section>"""
    main = f"""
    <h1>TideCluster report — {esc(model["meta"]["prefix"])}</h1>
    {cards}
    <p class="tc-callout">This is report v2. The original reports are
    still available alongside: <a href="../{esc(model["paths"].get("index_html") or "")}">index.html</a>,
    <a href="../{esc(model["paths"].get("tarean_report_html") or "")}">TAREAN</a>,
    and KITE page. Nothing was modified.</p>
    """
    (out_dir / "index.html").write_text(page_shell("Summary", "index", run_meta, main))


def _all_trcs_rows(model):
    rows = []
    for t in model["trcs"]:
        tarean = t["tarean"] or {}
        kite = t["kite"] or {}
        hor_cell = ""
        if kite:
            hor_cell = (f'{hor_count_cell(kite["n_hor_dominant"], "dom")}'
                        f'{hor_count_cell(kite["n_hor_visible"], "vis")}'
                        f'{hor_count_cell(kite["n_no_hor"], "none")}')
        sf_cell = (f'<a href="superfamilies.html#sf-{t["superfamily"]}">SF {t["superfamily"]}</a>'
                   if t["superfamily"] else "")
        rows.append(
            "<tr>"
            f'<td>{_trc_link(t["id"])}</td>'
            f'<td>{esc(t["repeat_type"] or "")}</td>'
            f'<td>{esc(t["annotation"] or "")}</td>'
            f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
            f'<td data-order="{t["total_size"]}">{fmt_bp(t["total_size"])}</td>'
            f'<td data-order="{t.get("min_array") or 0}">{fmt_bp(t.get("min_array"))}</td>'
            f'<td data-order="{t.get("max_array") or 0}">{fmt_bp(t.get("max_array"))}</td>'
            f'<td data-order="{tarean.get("monomer_length") or 0}">'
              f'{esc(tarean.get("monomer_length")) if tarean.get("monomer_length") else ""}</td>'
            f'<td>{hor_cell}</td>'
            f'<td>{sf_cell}</td>'
            "</tr>")
    return "".join(rows)


def render_trcs(model, out_dir, run_meta):
    body = f"""
    <h1>All tandem repeat clusters</h1>
    <p>{len(model["trcs"])} TRCs. Click a TRC to open its dashboard.</p>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="50"
           data-order='[[3,"desc"]]'>
      <thead><tr>
        <th>TRC</th><th>Type</th><th>Annotation</th>
        <th>Arrays</th><th>Total size</th><th>Min array</th><th>Max array</th>
        <th>Monomer</th><th>HOR&nbsp;(dom/vis/none)</th><th>Superfamily</th>
      </tr></thead>
      <tbody>{_all_trcs_rows(model)}</tbody>
    </table>
    </div>
    """
    (out_dir / "trcs.html").write_text(page_shell("All TRCs", "trcs", run_meta, body))


def render_tarean(model, out_dir, run_meta):
    rows = []
    for t in model["trcs"]:
        if not t["tarean"]:
            continue
        ta = t["tarean"]
        graph = _rel_img(ta.get("graph_png"), width=100, alt="k-mer graph")
        logo  = _rel_img(ta.get("logo_png"),  width=200, alt="logo")
        cons_short = esc((ta.get("consensus") or "")[:80])
        if ta.get("consensus") and len(ta["consensus"]) > 80:
            cons_short += "…"
        rows.append(
            "<tr>"
            f'<td>{_trc_link(t["id"])}</td>'
            f'<td data-order="{ta.get("monomer_length") or 0}">{esc(ta.get("monomer_length"))}</td>'
            f'<td data-order="{ta.get("total_score") or 0}">{esc(round(ta["total_score"], 4)) if ta.get("total_score") is not None else ""}</td>'
            f'<td data-order="{ta.get("total_size") or 0}">{fmt_bp(ta.get("total_size"))}</td>'
            f'<td>{esc(ta.get("annotation") or "")}</td>'
            f'<td data-order="{ta.get("n_arrays") or 0}">{esc(ta.get("n_arrays"))}</td>'
            f'<td><code>{cons_short}</code></td>'
            f'<td>{graph}</td>'
            f'<td>{logo}</td>'
            "</tr>")
    body = f"""
    <h1>TAREAN-analysed TRCs</h1>
    <p>{sum(1 for t in model["trcs"] if t["tarean"])} TRCs passed the
    minimum total length threshold and were processed by TAREAN.</p>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[3,"desc"]]'>
      <thead><tr>
        <th>TRC</th><th>Monomer</th><th>Score</th><th>Total size</th>
        <th>Annotation</th><th>Arrays</th><th>Consensus (preview)</th>
        <th>Graph</th><th>Logo</th>
      </tr></thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
    </div>
    """
    (out_dir / "tarean.html").write_text(page_shell("TAREAN", "tarean", run_meta, body))


def render_kite(model, out_dir, run_meta):
    hor_total = {"dom": 0, "vis": 0, "none": 0}
    for t in model["trcs"]:
        if t["kite"]:
            hor_total["dom"]  += t["kite"]["n_hor_dominant"]
            hor_total["vis"]  += t["kite"]["n_hor_visible"]
            hor_total["none"] += t["kite"]["n_no_hor"]
    total = sum(hor_total.values()) or 1
    bar = ""
    for kind, label in (("dom", "HOR-dominant"), ("vis", "HOR-visible"), ("none", "No HOR")):
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
        profile = _rel_img(k.get("profile_top3_png"), width=300, alt="profile")
        monomer = (t["tarean"] or {}).get("monomer_length")
        rows.append(
            "<tr>"
            f'<td>{_trc_link(t["id"])}</td>'
            f'<td data-order="{monomer or 0}">{esc(monomer)}</td>'
            f'<td data-order="{t["n_arrays"]}">{t["n_arrays"]}</td>'
            f'<td data-order="{k["n_hor_dominant"]}">{hor_count_cell(k["n_hor_dominant"],"dom")}</td>'
            f'<td data-order="{k["n_hor_visible"]}">{hor_count_cell(k["n_hor_visible"],"vis")}</td>'
            f'<td data-order="{k["n_no_hor"]}">{hor_count_cell(k["n_no_hor"],"none")}</td>'
            f'<td>{profile}</td>'
            "</tr>")
    body = f"""
    <h1>K-mer interval tandem repeat estimation (KITE)</h1>
    <h3>HOR distribution across all analysed arrays</h3>
    <div class="tc-bar" title="Total arrays: {total}">{bar}</div>
    <p style="font-size:12px;color:var(--fg-muted)">
      {hor_badge("HOR-dominant")} {hor_total["dom"]}
      · {hor_badge("HOR-visible")} {hor_total["vis"]}
      · {hor_badge("No HOR detected")} {hor_total["none"]}
      (total {total} arrays).
    </p>
    <h2>Per-TRC HOR summary</h2>
    <div class="tc-table-wrap">
    <table class="tc-table tc-datatable" data-page-length="25" data-order='[[3,"desc"]]'>
      <thead><tr>
        <th>TRC</th><th>Monomer<br>(TAREAN)</th><th>Arrays</th>
        <th>HOR-dom</th><th>HOR-vis</th><th>No&nbsp;HOR</th>
        <th>Profile</th>
      </tr></thead>
      <tbody>{"".join(rows)}</tbody>
    </table>
    </div>
    """
    (out_dir / "kite.html").write_text(page_shell("KITE", "kite", run_meta, body))


def render_superfamilies(model, out_dir, run_meta):
    if not model["superfamilies"]:
        body = ('<h1>TRC Superfamilies</h1>'
                '<p class="tc-callout">No TRC superfamilies were identified '
                '(either insufficient TRCs reached the clustering step or none shared '
                'significant consensus similarity).</p>')
        (out_dir / "superfamilies.html").write_text(
            page_shell("Superfamilies", "sf", run_meta, body))
        return
    sections = []
    for sf in model["superfamilies"]:
        peers = ", ".join(_trc_link(t) for t in sf["trcs"])
        img = (f'<div class="tc-fig"><a href="../{esc(sf["dotplot"])}">'
               f'<img src="../{esc(sf["dotplot"])}" width="500" alt="dotplot"></a>'
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
    (out_dir / "superfamilies.html").write_text(
        page_shell("Superfamilies", "sf", run_meta, body))


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

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Directory containing <prefix>_*.gff3 etc.")
    p.add_argument("--prefix",
                   help="Prefix used when the pipeline ran (auto-detect if only one).")
    p.add_argument("--output-dir", type=Path,
                   help="Output directory; defaults to <input-dir>/<prefix>_report_v2/.")
    args = p.parse_args(argv)

    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        p.error(f"input-dir is not a directory: {input_dir}")
    prefix = args.prefix or detect_prefix(input_dir)
    out_dir = args.output_dir or (input_dir / f"{prefix}_report_v2")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "data").mkdir(exist_ok=True)

    print(f"input:  {input_dir}")
    print(f"prefix: {prefix}")
    print(f"output: {out_dir}")

    model = build_model(input_dir, prefix)

    # Stats one-liner for sanity.
    trc_count = len(model["trcs"])
    analysed  = sum(1 for t in model["trcs"] if t.get("tarean"))
    sf_count  = len(model["superfamilies"])
    hor_dom   = sum((t["kite"]["n_hor_dominant"] if t.get("kite") else 0) for t in model["trcs"])
    hor_vis   = sum((t["kite"]["n_hor_visible"]  if t.get("kite") else 0) for t in model["trcs"])
    print(f"TRCs: {trc_count}  TAREAN-analysed: {analysed}  "
          f"superfamilies: {sf_count}  HOR dom/vis: {hor_dom}/{hor_vis}")

    script_path = Path(__file__).resolve()
    copy_assets(out_dir / "assets", script_path)

    with (out_dir / "data" / "report.json").open("w") as f:
        json.dump(model, f, indent=2, default=str)
    print(f"wrote: {out_dir / 'data' / 'report.json'}")

    run_meta = make_run_meta(model)
    render_index(model, out_dir, run_meta)
    render_trcs(model, out_dir, run_meta)
    render_tarean(model, out_dir, run_meta)
    render_kite(model, out_dir, run_meta)
    render_superfamilies(model, out_dir, run_meta)
    print(f"rendered: index, trcs, tarean, kite, superfamilies")


if __name__ == "__main__":
    main()
