# TideCluster — Claude Code Project Notes

TideCluster (`version.py` → 1.9.3) is a Python + R pipeline that wraps
TideHunter, mmseqs2, BLAST, RepeatMasker and TAREAN to detect, cluster
and characterise tandem repeats in genome assemblies.

## Sandbox: Hermit (Singularity)

This project is developed **inside a Hermit Singularity sandbox**.
The sandbox itself lives in `./hermit/` and is **not part of the
project source tree** — do not edit anything under `hermit/` as part
of TideCluster changes.

Key sandbox facts (see `hermit/CLAUDE.md` and `hermit/envs/hermit.env`):

- Project root `/home/petr/PycharmProjects/TideCluster` is bound
  read-write inside the container (it is the hermit `DATA_OUTPUT_DIR`).
- `$AGENT_CPUS` / `$AGENT_MEMORY` advertise the available budget;
  use `$((AGENT_CPUS - 1))` for `-c` / thread flags.
- Bare `conda install`, `mamba install` and `pip install` are
  **blocked**. Use `htool <name>` for single tools, or
  `mamba create -p /envs/conda/envs/<name>` for multi-package envs.
- Persistent conda envs live under `/envs/conda/envs/`; pip cache
  under `/envs/pip`; npm under `/envs/npm_global`. Anything you
  install there survives sandbox restarts.

## Runtime environment for development

Dependencies are pinned in `conda-deps.txt` (intentionally not named
`requirements.txt` so Dependabot ignores it). They are installed in
a single env named `tidecluster`:

```bash
# Already created in this sandbox at /envs/conda/envs/tidecluster
conda activate tidecluster        # name resolves via $CONDA_ENVS_PATH
# or explicitly: conda activate /envs/conda/envs/tidecluster
```

The env supplies: TideHunter 1.4.3, mmseqs2, NCBI BLAST,
RepeatMasker 4.1.2.p1, mafft, bedtools, networkx, kitehor 0.13.2
(`-c petrnovak`), and the R/Bioconductor stack (Biostrings,
GenomicRanges, rtracklayer, GenomeInfoDbData, igraph, optparse,
plyr/dplyr, reshape2, scales, hwriter, R2HTML, jsonlite).

To add a dependency: edit `conda-deps.txt` (the source of truth), then
update the env. Update `TideCluster.def` if the change should also reach
the Singularity image.

**Reality check on installing into the env (learned the hard way):**
the hermit hook **blocks `mamba install -n tidecluster ...`** despite
older notes calling it "the documented exception". Two ways through:

- **Don't touch the shared env — build an isolated one** (the allowed
  `mamba create -p` path) and PATH-shadow it for the run/test:
  ```bash
  mamba create -y -p /envs/conda/envs/<iso> -c conda-forge -c bioconda \
      -c petrnovak --file conda-deps.txt        # or a single tool
  conda activate tidecluster
  export PATH=/envs/conda/envs/<iso>/bin:$PATH  # the iso binary wins
  ```
  `tests.sh` skips its own `conda activate` when `CONDA_DEFAULT_ENV` is
  already set, so this PATH-shadow survives into the test run. This is
  how to test a new kitehor/tool version without mutating the dev env.
- **Or** have the user run the install themselves outside the hook
  (`! mamba install -n tidecluster ... <pkg>`).

The dev `tidecluster` env can lag the `conda-deps.txt` pin (it has, e.g.
kitehor 0.9.3 while the pin moved to 0.12/0.13) — check
`kitehor --version` before assuming the pinned version is what runs.

## Entry points

| Script                              | Purpose                                                  |
|-------------------------------------|----------------------------------------------------------|
| `TideCluster.py`                    | Main CLI: `tidehunter`, `clustering`, `annotation`, `tarean`, `run_all` |
| `tc_reannotate.py`                  | RepeatMasker reannotation using a TRC library            |
| `tc_per_tra_consensus.py`           | Per-TRA consensus + grading post-processing              |
| `tc_merge_annotations.py`           | Merge annotation TSVs across samples                     |
| `tc_update_gff3.py`                 | Rename TRCs in a GFF3 via conversion table               |
| `tc_rerender_report.py`             | Regenerate HTML report from an existing run              |
| `tc_comparative_analysis.R`         | Cross-sample TRC → satellite family clustering           |
| `tc_summarize_comparative_analysis.R` | HTML report for comparative analysis                   |
| `tc_utils.py`                       | Shared library (FASTA chunking, GFF3 IO, clustering glue)|

R helpers live in `tarean/` (TAREAN, KITE, MSA consensus, plotting).

## Tests

`./tests.sh {smoke|short|long|rerender|determinism|all}` —
auto-activates the `tidecluster` conda env. `NCPU` (or positional arg)
controls thread count; default is 2. Use `smoke` for quick sanity
checks; `long` runs `run_all` on the bundled
`test_data/CEN6_ver_220406.fasta`.

- The **kite path** (KITE/tarean) is exercised by `short` (and above),
  **not** `smoke`. To test it against a non-default kitehor, PATH-shadow
  an isolated env (see "Reality check" above) before calling `tests.sh`.
- `determinism` covers the comparative-analysis reproducibility
  regression; it needs a multi-sample fixture and **skips** if none is
  present (looks for `$TC_COMPARATIVE_FIXTURE`, `tests/data/comparative`,
  then `test_data/analysis_1.10.5`).

### Gotchas

- `<prefix>_kite/monomer_size_top3_estimats.csv` is **TAB-delimited**
  despite the `.csv` name — read it with `delimiter='\t'` /
  `sep="\t"`, not comma.
- The sandbox has **no `/usr/bin/time` and no `bc`**. Time things with
  `date +%s%N` (integer-ns deltas) and do arithmetic in Python, not `bc`.
- Founder / strongest / multiplicity / HOR / subrepeat have precise,
  easy-to-conflate meanings. Before touching
  `tc_utils.build_monomer_size_csv` (the 5-pass + Pass-6 kh_deeper
  founder logic) or the subrepeat code, read the canonical definitions:
  the README "Concepts" / KITE section, and the persistent memory note
  on founder/strongest/HOR semantics.

## Releasing

Releases are **annotated git tags** named with a bare version
(`1.13.1`); pushing one fires `release.yml` (SIF → GHCR + GitHub
Release) and `conda-release.yml` (→ anaconda.org, which **asserts the
tag name equals `version.py`**). The convention:

1. Bump `version.py`, prepend a `## X.Y.Z (date)` block to
   `changelog.md`, bump the README install/pull version strings
   (`tidecluster=…`, `sif:…`).
2. Commit `release: X.Y.Z — <summary>` **directly on `main`** (release
   commits go on main; feature work lands first via fast-forward, no
   merge commits).
3. `git tag -a X.Y.Z -m "…" HEAD`, with the tag name == `version.py`.
4. **Push from the host, not the sandbox** (no SSH trust here):
   `git push origin main && git push origin X.Y.Z`.

The `/release` skill automates steps 1–3 and prints the push commands.

## Singularity image

`TideCluster.def` builds `TideCluster.sif` from the same
`conda-deps.txt`. Build with:

```bash
sudo singularity build TideCluster.sif TideCluster.def
```

Released images are published to
`oras://ghcr.io/kavonrtep/tidecluster/sif:<tag>`.

## Don'ts

- Don't edit `hermit/` as part of project work.
- Don't run bare `conda install` / `pip install` (hook-blocked). Keep
  `conda-deps.txt` the source of truth and update envs via
  `mamba create -p` isolated envs (see "Reality check on installing"),
  or have the user run the install outside the hook.
- Don't push git tags/commits from the sandbox (SSH host-key fails) —
  give the user the `git push` commands to run on the host.
- Don't rename `conda-deps.txt` to `requirements.txt`.
