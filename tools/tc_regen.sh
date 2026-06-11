#!/bin/bash
# tools/tc_regen.sh <run_dir> [prefix]
#
# One-shot: regenerate the KITE monomer-size CSV from a run's existing kitehor
# output, then re-render the HTML report. Collapses the regenerate-CSV ->
# re-render loop (the most-repeated dev action) into a single command, with the
# kitehor PATH-shadow baked in.
#
#   tools/tc_regen.sh test_data/Solanum_lycopersicum
#   tools/tc_regen.sh test_data/Arabidopsis_thaliana/run_e2e tc
#
# It locates <run_dir>/<prefix>_kite (the kitehor.*.tsv) and
# <run_dir>/<prefix>_clustering.gff3, rebuilds <prefix>_kite/
# monomer_size_top3_estimats.csv via tc_utils.build_monomer_size_csv, and runs
# tc_rerender_report.py --input-dir <run_dir> --prefix <prefix>.
#
# ── Updating kitehor ──────────────────────────────────────────────────────────
# KITE binaries come from a dedicated isolated env, PATH-shadowed below so the
# (often-lagging) dev `tidecluster` env is not mutated. When a new kitehor ships:
#   1) build the isolated env once (host or via the allowed `mamba create -p`):
#        mamba create -y -p /envs/conda/envs/kitehor0XYZ -c conda-forge \
#            -c bioconda -c petrnovak kitehor=0.XY.Z
#   2) bump KITEHOR_ENV below (or export KITEHOR_ENV=/envs/conda/envs/kitehor0XYZ).
# Set KITEHOR_ENV="" to use whatever kitehor is already on PATH (e.g. once the
# dev `tidecluster` env itself is updated — then no shadow is needed).
set -euo pipefail

KITEHOR_ENV="${KITEHOR_ENV-/envs/conda/envs/kitehor0132}"   # current pinned KITE env

RUN_DIR="${1:?usage: tools/tc_regen.sh <run_dir> [prefix]}"
PREFIX="${2:-tc}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Activate the dev env (for R / mmseqs / etc. the renderer may touch) unless one
# is already active, then PATH-shadow the pinned kitehor.
if [ "${CONDA_DEFAULT_ENV:-}" = "" ]; then
  eval "$(conda shell.bash hook)" && conda activate tidecluster
fi
if [ -n "$KITEHOR_ENV" ] && [ -x "$KITEHOR_ENV/bin/kitehor" ]; then
  export PATH="$KITEHOR_ENV/bin:$PATH"
fi
echo "[tc_regen] kitehor -> $(command -v kitehor || echo none) ($(kitehor --version 2>/dev/null || echo n/a))"

KITE_DIR="$RUN_DIR/${PREFIX}_kite"
GFF3="$RUN_DIR/${PREFIX}_clustering.gff3"
[ -d "$KITE_DIR" ] || { echo "[tc_regen] no kite dir: $KITE_DIR" >&2; exit 1; }

echo "[tc_regen] rebuilding $KITE_DIR/monomer_size_top3_estimats.csv"
python3 - "$KITE_DIR" "$GFF3" <<'PY'
import os, sys, tc_utils as tc
kd, gff3 = sys.argv[1], sys.argv[2]
def opt(p): return p if os.path.exists(p) else None
n = tc.build_monomer_size_csv(
    kite_tsv=opt(f"{kd}/kitehor.kite.tsv"),
    ssr_tsv=opt(f"{kd}/kitehor.ssr.tsv"),
    rescored_peaks_tsv=f"{kd}/kitehor.rescored.peaks.tsv",
    out_csv=f"{kd}/monomer_size_top3_estimats.csv",
    tandem_validate_tsv=opt(f"{kd}/kitehor.tandem_validate.tsv"),
    trc_repeat_type=tc.parse_trc_ssr_motif_len(gff3) if os.path.exists(gff3) else None,
)
print(f"[tc_regen] wrote {n} array rows")
PY

echo "[tc_regen] re-rendering report"
python3 ./tc_rerender_report.py --input-dir "$RUN_DIR" --prefix "$PREFIX"
echo "[tc_regen] done: $RUN_DIR/${PREFIX}_report/"
