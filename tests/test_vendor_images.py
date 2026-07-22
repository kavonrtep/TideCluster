#!/usr/bin/env python3
"""Unit test for self-contained report image vendoring.

The report v2 used to reference PNGs by reaching outside <prefix>_report/ into
the kite / tarean / dotplots scratch trees via `../` paths, so deleting those
(multi-GB) trees broke the report's images. build_report now copies exactly the
referenced images into <prefix>_report/img/ under a flattened layout and rewrites
the refs, so the report survives a purge of the source trees.

This test covers the three moving parts without a full render:
  1. _vendor_map  — the input-dir-relative source path -> flattened dest mapping
     for each source tree (tarean / kite / dotplots) plus the misc/ fallback.
  2. _img_url     — records src_rel -> dest in the shared collector and returns
     the per-page href (img_href + dest).
  3. copy_vendored_images — copies each referenced image to img/<dest>, flattens
     correctly, rebuilds from scratch, and tolerates a missing source.

Run: python3 tests/test_vendor_images.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import tc_rerender_report as r  # noqa: E402
from pathlib import Path  # noqa: E402


def test_vendor_map():
    cases = {
        # TAREAN: keyed by the <TRC>.fasta_tarean segment, namespaced per TRC.
        "tc_tarean/TRC_1.fasta_tarean/img/graph_11mer_1.png":
            "tarean/TRC_1/graph_11mer_1.png",
        "tc_tarean/TRC_12.fasta_tarean/img/logo_23mer_2.png":
            "tarean/TRC_12/logo_23mer_2.png",
        # KITE profile plots: flattened to kite/<basename> (TRC id is in the name).
        "tc_kite/profile_plots/profile_TRC_1.png": "kite/profile_TRC_1.png",
        "tc_kite/profile_plots/profile_top3_TRC_9.png": "kite/profile_top3_TRC_9.png",
        # Superfamily dotplots.
        "dotplots/TRCS_1_4_7_15.png": "dotplots/TRCS_1_4_7_15.png",
        # Unknown shape -> stays self-contained under misc/.
        "some/other/thing.png": "misc/thing.png",
    }
    for src, want in cases.items():
        got = r._vendor_map(src)
        assert got == want, f"_vendor_map({src!r}) = {got!r}, want {want!r}"
    # A non-default kite dir name still routes via the /profile_plots/ marker.
    assert r._vendor_map("KITE_out/profile_plots/profile_TRC_2.png") == \
        "kite/profile_TRC_2.png"
    print("  _vendor_map: tarean / kite / dotplots / misc mappings OK")


def test_img_url_registers_and_builds_href():
    for img_href in ("tc_report/img/", "img/", "../img/"):
        ctx = {"img_href": img_href, "img_map": {}}
        url = r._img_url(ctx, "tc_kite/profile_plots/profile_TRC_1.png")
        assert url == f"{img_href}kite/profile_TRC_1.png", url
        # src_rel -> flattened dest recorded for the copy step.
        assert ctx["img_map"] == {
            "tc_kite/profile_plots/profile_TRC_1.png": "kite/profile_TRC_1.png"
        }
    # Empty / falsy src_rel -> no URL, nothing recorded.
    ctx = {"img_href": "img/", "img_map": {}}
    assert r._img_url(ctx, None) is None
    assert r._img_url(ctx, "") is None
    assert ctx["img_map"] == {}
    print("  _img_url: records collector entry + builds per-page href OK")


def test_copy_vendored_images():
    tmp = Path(tempfile.mkdtemp(prefix="tc_vendor_test_"))
    input_dir = tmp / "run"
    # Create three source images mirroring real tree layouts + one missing.
    srcs = {
        "tc_kite/profile_plots/profile_TRC_1.png": b"KITE",
        "tc_tarean/TRC_1.fasta_tarean/img/graph_11mer_1.png": b"GRAPH",
        "dotplots/TRCS_2_16.png": b"DOTPLOT",
    }
    for rel, data in srcs.items():
        p = input_dir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)

    # The collector, as build_report would hand it over (incl. one missing src).
    img_map = {rel: r._vendor_map(rel) for rel in srcs}
    img_map["tc_kite/profile_plots/profile_TRC_404.png"] = "kite/profile_TRC_404.png"

    img_dir = tmp / "run" / "tc_report" / "img"
    # Pre-seed a stale file to prove the dir is rebuilt from scratch each run.
    img_dir.mkdir(parents=True, exist_ok=True)
    (img_dir / "STALE.png").write_bytes(b"old")

    stats = r.copy_vendored_images(input_dir, img_dir, img_map, quiet=True)

    assert stats["copied"] == 3, stats
    assert stats["missing"] == 1, stats
    # Flattened copies exist with intact content.
    assert (img_dir / "kite/profile_TRC_1.png").read_bytes() == b"KITE"
    assert (img_dir / "tarean/TRC_1/graph_11mer_1.png").read_bytes() == b"GRAPH"
    assert (img_dir / "dotplots/TRCS_2_16.png").read_bytes() == b"DOTPLOT"
    # Missing source is skipped, not copied.
    assert not (img_dir / "kite/profile_TRC_404.png").exists()
    # Stale pre-existing file was cleared by the rmtree rebuild.
    assert not (img_dir / "STALE.png").exists()
    print("  copy_vendored_images: flatten + rebuild + missing-tolerant OK")


def test_copy_tarean_drilldowns():
    tmp = Path(tempfile.mkdtemp(prefix="tc_dd_test_"))
    input_dir = tmp / "run"
    tdir = input_dir / "tc_tarean" / "TRC_1.fasta_tarean"
    (tdir / "img").mkdir(parents=True)
    (tdir / "report.html").write_bytes(b"<html>report</html>")
    (tdir / "img" / "graph_7mer_1.png").write_bytes(b"G")
    (tdir / "ppm_7mer_1.csv").write_bytes(b"ppm")
    # Heavy intermediates that must NOT be copied.
    (tdir / "TRC_1.fasta").write_bytes(b"X" * 1000)
    (tdir / "TRC_1.fasta_7.kmers").write_bytes(b"Y" * 1000)
    # A TAREAN TRC with no report.html must be skipped, not crash.
    (input_dir / "tc_tarean" / "TRC_2.fasta_tarean").mkdir(parents=True)

    model = {
        "paths": {"tarean_dir": "tc_tarean"},
        "trcs": [
            {"id": "TRC_1", "tarean": {"has_report_html": True,
                                       "tarean_dir": "TRC_1.fasta_tarean"}},
            {"id": "TRC_2", "tarean": {"has_report_html": False,
                                       "tarean_dir": "TRC_2.fasta_tarean"}},
        ],
    }
    out_dir = tmp / "run" / "tc_report"
    stats = r.copy_tarean_drilldowns(input_dir, out_dir, model, quiet=True)
    assert stats["copied"] == 1, stats
    dd = out_dir / "tarean" / "TRC_1"
    assert (dd / "report.html").read_bytes() == b"<html>report</html>"
    assert (dd / "img" / "graph_7mer_1.png").read_bytes() == b"G"
    assert (dd / "ppm_7mer_1.csv").read_bytes() == b"ppm"
    # Heavy files left behind.
    assert not (dd / "TRC_1.fasta").exists()
    assert not (dd / "TRC_1.fasta_7.kmers").exists()
    # The report.html-less TRC produced no drill-down dir.
    assert not (out_dir / "tarean" / "TRC_2").exists()
    print("  copy_tarean_drilldowns: light-only copy + skip-missing OK")


def test_vendor_legacy_report():
    tmp = Path(tempfile.mkdtemp(prefix="tc_legacy_test_"))
    input_dir = tmp / "run"
    legacy = input_dir / "tc_report_legacy"
    legacy.mkdir(parents=True)
    (legacy / "tc_index.html").write_text(
        '<a href="tc_kite_report.html"><strong>Report</strong></a>'
        '<a href="tc_tarean_report.html">TAREAN</a>')
    # Mirror the real v1 hwriter output: mixed quotes, and _move_v1_to_legacy
    # only prepended ../ to DOUBLE-quoted refs, so single-quoted ones (the
    # per-TRC report.html drill-down anchors) arrive WITHOUT ../.
    (legacy / "tc_tarean_report.html").write_text(
        '<img src="../tc_tarean/TRC_1.fasta_tarean/img/graph_7mer_1.png">'   # double, ../
        "<a href='tc_tarean/TRC_1.fasta_tarean/report.html'>t1</a>"          # single, no ../
        '<a href="../tc_tarean/TRC_2.fasta_tarean/report.html">t2</a>'       # double, ../
        '<img src="../dotplots/TRCS_1_2.png">'
        '<a href="https://example.org/x.png">ext</a>')

    img_map = {}
    stats = r._vendor_legacy_report(input_dir, "tc", img_map, quiet=True)
    assert stats["rewritten"] == 2, stats
    tr = (legacy / "tc_tarean_report.html").read_text()
    # Image refs -> shared vendored img/, and added to img_map for copying.
    assert 'src="../tc_report/img/tarean/TRC_1/graph_7mer_1.png"' in tr
    assert 'src="../tc_report/img/dotplots/TRCS_1_2.png"' in tr
    assert img_map["tc_tarean/TRC_1.fasta_tarean/img/graph_7mer_1.png"] == \
        "tarean/TRC_1/graph_7mer_1.png"
    assert img_map["dotplots/TRCS_1_2.png"] == "dotplots/TRCS_1_2.png"
    # Per-TRC report.html -> vendored drill-down, BOTH the single-quoted/no-../
    # anchor (the bug) and the double-quoted/../ one.
    assert "href='../tc_report/tarean/TRC_1/report.html'" in tr
    assert 'href="../tc_report/tarean/TRC_2/report.html"' in tr
    assert "tc_tarean/" not in tr.replace("tc_report", "")  # no raw-tree ref left
    # External URL untouched.
    assert 'href="https://example.org/x.png"' in tr
    # Dangling sibling link unwrapped (label kept); valid sibling link intact.
    idx = (legacy / "tc_index.html").read_text()
    assert 'href="tc_kite_report.html"' not in idx and "<strong>Report</strong>" in idx
    assert 'href="tc_tarean_report.html"' in idx
    print("  _vendor_legacy_report: repoint imgs/report + drop dangling link OK")


if __name__ == "__main__":
    test_vendor_map()
    test_img_url_registers_and_builds_href()
    test_copy_vendored_images()
    test_copy_tarean_drilldowns()
    test_vendor_legacy_report()
    print("test_vendor_images: PASSED")
