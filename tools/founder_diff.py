#!/usr/bin/env python3
"""tools/founder_diff.py <run_dir> [run_dir ...]

Founder-invariance check for a pipeline change: rebuild the KITE monomer CSV
with the *working-tree* tc_utils.py and with the *git HEAD* tc_utils.py, then
diff the founder-bearing columns. Prints, per run, the number of changed cells
and the added/removed columns — so a "this change is additive / display-only"
claim can be verified in one command.

  tools/founder_diff.py test_data/Solanum_lycopersicum
  tools/founder_diff.py test_data/Solanum_lycopersicum test_data/Arabidopsis_thaliana/run_e2e

Exit 0 if every run has zero founder-column diffs, 1 otherwise. Run it before
committing a tc_utils change to confirm you only touched what you meant to.
"""
import csv
import glob
import importlib.util
import os
import subprocess
import sys
import tempfile

KEY = ["TRC_ID", "seqid", "start", "end", "founder_period", "strongest_period",
       "multiplicity", "irregular_multiplicity", "founder_method",
       "hor_order_confidence"]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _loadmod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.argv = ["x"]                       # tc_utils parses argv at import on some paths
    spec.loader.exec_module(mod)
    return mod


def _gen(mod, run_dir):
    kd = os.path.join(run_dir, "tc_kite")
    cl = glob.glob(os.path.join(run_dir, "*_clustering.gff3"))
    out = tempfile.mktemp(suffix=".tsv")
    opt = lambda p: p if os.path.exists(p) else None
    mod.build_monomer_size_csv(
        kite_tsv=opt(f"{kd}/kitehor.kite.tsv"),
        ssr_tsv=opt(f"{kd}/kitehor.ssr.tsv"),
        rescored_peaks_tsv=f"{kd}/kitehor.rescored.peaks.tsv",
        out_csv=out,
        tandem_validate_tsv=opt(f"{kd}/kitehor.tandem_validate.tsv"),
        trc_repeat_type=mod.parse_trc_ssr_motif_len(cl[0]) if cl else None,
    )
    return out


def main(run_dirs):
    head_src = os.path.join(tempfile.mkdtemp(), "tc_utils_head.py")
    with open(head_src, "w") as fh:
        fh.write(subprocess.check_output(
            ["git", "-C", ROOT, "show", "HEAD:tc_utils.py"], text=True))
    head = _loadmod(head_src, "tc_head")
    work = _loadmod(os.path.join(ROOT, "tc_utils.py"), "tc_work")

    bad = 0
    for rd in run_dirs:
        if not os.path.exists(os.path.join(rd, "tc_kite",
                                           "kitehor.rescored.peaks.tsv")):
            print(f"{rd}: SKIP (no tc_kite/kitehor.rescored.peaks.tsv)")
            continue
        hr = list(csv.DictReader(open(_gen(head, rd)), delimiter="\t"))
        wr = list(csv.DictReader(open(_gen(work, rd)), delimiter="\t"))
        diffs = sum(1 for a, b in zip(hr, wr) for k in KEY if a.get(k) != b.get(k))
        new = [c for c in wr[0] if c not in hr[0]]
        gone = [c for c in hr[0] if c not in wr[0]]
        flag = "" if diffs == 0 else "  <-- FOUNDER COLUMNS CHANGED"
        print(f"{rd}: founder-col diffs={diffs}  +cols={new}  -cols={gone}{flag}")
        bad += (diffs != 0)
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
