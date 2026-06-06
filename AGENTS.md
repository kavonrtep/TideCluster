# Repository Guidelines

## Project Structure & Module Organization
`TideCluster.py` is the main CLI entry point for the tandem-repeat pipeline (`tidehunter`, `clustering`, `annotation`, `tarean`, `run_all`). Shared Python helpers live in `tc_utils.py`; small maintenance utilities include `tc_update_gff3.py`, `tc_reannotate.py`, and `tc_merge_annotations.py`. R-based analysis and reporting live in `tarean/`. Report assets are stored in `html/`, `html_template/`, and `js/`. Use `test_data/`, `micro_pea.fasta`, and `tiny_pea.fasta` for local validation; write scratch outputs under `tmp/`.

## Build, Test, and Development Commands
Dependencies are pinned in `conda-deps.txt` (intentionally **not**
`requirements.txt`, so Dependabot ignores it) and installed in the
`tidecluster` conda env; the external tools must be on `PATH`.

```bash
python TideCluster.py --help
python TideCluster.py run_all -c 4 -pr tmp/demo -f micro_pea.fasta
python TideCluster.py clustering -c 4 -pr tmp/demo -f micro_pea.fasta -m 5000
# Comparative analysis: -i is a TAB-delimited table with columns
# input_dir / sample_code / tidecluster_prefix (one row per sample).
Rscript tc_comparative_analysis.R -i input_table.tsv -o tmp/comparative -c 4
python -m py_compile TideCluster.py tc_utils.py tc_update_gff3.py tidecluster_viz.py
```

Use `run_all` for an end-to-end smoke test, stage-specific commands for focused debugging, and `py_compile` for a fast syntax check before review.

## Coding Style & Naming Conventions
Match the surrounding style instead of reformatting unrelated code. Python uses 4-space indentation, `snake_case` for functions, and `CamelCase` for classes such as `GRange` and `Gff3Feature`. R scripts also favor `snake_case` and `optparse`-based CLIs. Keep filename and output naming prefix-based, for example `sample_clustering.gff3`, `sample_tarean/`, and `sample_consensus/`.

## Testing Guidelines
Run `./tests.sh {smoke|short|long|rerender|determinism|all}` (it auto-activates the `tidecluster` env; `NCPU`/positional arg sets threads). Tiered Python checks live in `tests/` and run as scripts (e.g. `python3 tests/test_strongest_by_identity.py`), not under `pytest`. Note `short` (not `smoke`) exercises the KITE path. Also validate changes with small bundled inputs and confirm the expected artifacts for the affected step, such as `*_tidehunter.gff3`, `*_clustering.gff3`, `*_index.html`, or comparative-analysis tables under the chosen output directory.

## Commit & Pull Request Guidelines
Recent history favors short, focused subjects such as `version changed`, `some refactoring`, and `readme updated`; keep commit titles concise and specific to one change. Pull requests should summarize pipeline impact, list the commands you ran, note any required external tools or reference data, and include screenshots when HTML reports or plots change. Avoid committing large generated outputs or container artifacts unless the change explicitly requires them.
