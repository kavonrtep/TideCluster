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
import json
import os
import re
import shutil
import sys
from pathlib import Path

__version__ = "1"  # schema version for report.json

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


if __name__ == "__main__":
    main()
