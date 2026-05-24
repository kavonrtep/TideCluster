# KITE → kitehor integration plan

Branch: `feat/kite-kitehor` (off `main`).

This document is the implementation plan for replacing the in-tree R
KITE analysis (`tarean/kite.R`) with [kavonrtep/kitehor](https://github.com/kavonrtep/kitehor),
a Rust reimplementation that follows the same k-mer-interval principle
but emits substantially richer per-array structural information (HOR
verdict + multiplicity + founder/tile, subrepeat scan, SSR scan,
per-window classification, combined class label).

## Decisions taken upfront

| Question                                              | Decision                                                        |
|-------------------------------------------------------|-----------------------------------------------------------------|
| Per-array heatmaps (`profile_plots/profile_*.png`)    | Keep; render in **R** (`tarean/kite_heatmaps.R`, base R only) from kitehor's `--periodogram` bundle. No new conda deps. |
| Old `tarean/kite.R` analysis                          | Remove. kitehor is the only analysis backend                    |
| `tc_rerender_report.py` Python HOR port               | Remove. Read kitehor `.summary.tsv` directly                    |
| Conda distribution                                    | Add `kitehor=0.9.3` from `-c petrnovak` channel                 |
| Output column casing                                  | **Switch to kitehor's lowercase `hor_*`** in the joined CSV (clean break, all consumers updated in this branch) |
| `peaks_list.RDS`                                      | Drop. Profiles persisted as the kitehor `.periodogram` FASTA-like bundle. |
| Branch                                                | `feat/kite-kitehor`                                             |

## What kitehor produces (recap)

Command: `kitehor analyze <in.fa> -o <out_prefix> --threads N`.
Inputs: one multi-FASTA, **one record per tandem-repeat array**.
Outputs (sibling files at `<out_prefix>.*`):

| File                       | Role                                                                                  |
|----------------------------|---------------------------------------------------------------------------------------|
| `.kite.tsv`                | Per-array top-3 monomer-size + score (replaces old `monomer_size_top3_estimats.csv`)  |
| `.kite.peaks.tsv`          | Per-array ranked peaks (period, height, score, score2_norm, background)               |
| `.verdicts.tsv`            | HOR verdict + founder/tile/multiplicity/confidence/reason                             |
| `.subrepeat.tsv`           | Subrepeat / nested-TR flag with block coverage                                        |
| `.windows.tsv`             | Per-window class (sliding classifier output)                                          |
| `.ssr.tsv`                 | SSR scan summary per array                                                            |
| `.ssr.regions.tsv`         | SSR per-region detail                                                                 |
| `.hor_within_tile.tsv`     | HOR-within-tile diagnostics (founder_density, phase_contrast, ratios)                 |
| `.summary.tsv` (32 cols)   | **Canonical merged output** including `combined_class` ∈ {hor, hor_with_ssr, tr, tr_with_ssr, tr_with_nested_tr, tr_with_subrepeat, pure_ssr, unresolved} |

Parallelism: internal `rayon`; honors `--threads`. Pure Rust deps, no
system libs beyond optional `image` (behind `viz` feature, off by
default for our build).

## Input adaptation

The current pipeline materialises one fasta per TRC under
`{prefix}_tarean/fasta/TRC_<N>.fasta`, with sequences stored as
**dimers** (sequence concatenated with itself) and headers shaped like
`TRC_<N>:<seqid>_<start>_<end>`. `kite.R` halves the dimer via
`subseq(x, start=1, end=nchar(x)/2)` before analysis.

For kitehor we will:

1. Build **one multi-FASTA** `{prefix}_kite_input.fasta` by reading
   every `{prefix}_tarean/fasta/TRC_*.fasta`, halving each record back
   to the single-array sequence, and rewriting the header to
   `<TRC>:<seqid>_<start>_<end>` (unchanged from today).
2. Call `kitehor analyze {prefix}_kite_input.fasta -o {prefix}_kite/kitehor --periodogram {prefix}_kite/kitehor.periodogram --threads {cpu}`.
   Single invocation; emits the 9 per-stage TSVs **plus** the
   FASTA-like periodogram bundle containing per-array `H[d]` and
   `bg[d]` lines. `--periodogram` is a new flag added in kitehor 0.9.3.
3. Remove the temp multi-FASTA — kitehor's `record_id` keeps the
   TRC/seqid/start/end recoverable.

`record_id` from kitehor will equal the input FASTA header verbatim,
so the existing `TRC:seqid_start_end` split in `tc_rerender_report.py`
(`peaks_best$TRC_ID`, `$seqid`, `$start`, `$end` derivation) ports
trivially.

**Historical note**: The first internal branches (against kitehor
0.9.2) had to run kitehor twice — `analyze` for the structural TSVs,
then `kite-periodicity --dump-profile DIR` for the per-array spectrum
vectors. kitehor 0.9.3 forwards a single `--periodogram <PATH>` flag
through `analyze`, replacing both that two-call pattern and the
per-array TSV directory with one bundled output file.

**GLIBC compatibility**: Resolved in kitehor 0.9.3. The 0.9.2 conda
binary required glibc ≥2.34 and aborted on Ubuntu 20.04 / CentOS 8
hosts; 0.9.3 is rebuilt against an older target and installs cleanly
via `-c petrnovak` on those distributions.

## Output mapping

We keep the on-disk filename `monomer_size_top3_estimats.csv` (it is
referenced by `tc_per_tra_consensus.py:55` and by
`tc_rerender_report.py:489`) but the column casing **switches to
lowercase** to match kitehor's native schema, and the file is now
joined from kitehor's `.kite.tsv` + `.summary.tsv`:

| Column                    | Source                                                          |
|---------------------------|-----------------------------------------------------------------|
| `TRC_ID`                  | parsed from `record_id`                                         |
| `seqid`, `start`, `end`   | parsed from `record_id`                                         |
| `monomer_size`            | `.kite.tsv:monomer_size`                                        |
| `score`                   | `.kite.tsv:score`                                               |
| `array_length`            | `.kite.tsv:array_length`                                        |
| `monomer_size_2`, `score_2`, `monomer_size_3`, `score_3` | `.kite.tsv` cols of same name        |
| `hor_status`              | mapped from `.summary.tsv:hor_verdict` (see mapping below)      |
| `hor_confidence`          | `.summary.tsv:hor_confidence`                                   |
| `hor_founder`             | `.summary.tsv:hor_founder` (bp; was `HOR_base_monomer`)         |
| `hor_tile`                | `.summary.tsv:hor_tile` (bp; was `HOR_hor_period`)              |
| `hor_multiplicity`        | `.summary.tsv:hor_multiplicity` (k; was `HOR_n_harmonics`)      |
| `combined_class`          | `.summary.tsv:combined_class`                                   |
| `subrepeat_flag`          | `.summary.tsv:subrepeat_flag`                                   |
| `subrepeat_period_bp`     | `.summary.tsv:subrepeat_period_bp`                              |
| `ssr_flag`                | `.summary.tsv:ssr_flag`                                         |
| `ssr_dominant_motif`      | `.summary.tsv:ssr_dominant_motif`                               |
| `founder_density`         | `.summary.tsv:founder_density`                                  |
| `phase_contrast`          | `.summary.tsv:phase_contrast`                                   |

The full unmodified kitehor TSVs are also written under
`{prefix}_kite/kitehor.*.tsv` for downstream tools that want the raw
data. Per-array spectrum vectors land under
`{prefix}_kite/profiles/<record_id>.kite.tsv` (sparse `d \t H \t bg`).

### HOR status mapping

`hor_status` is the 4-bin string the v1 HTML report and v2 report
already use. kitehor's `hor_verdict` (one of `hor`, `simple_tr`,
`no_hor`, `unresolved`) is mapped:

| kitehor verdict / confidence                          | hor_status (legacy bin) |
|-------------------------------------------------------|-------------------------|
| `hor_verdict ∈ {no_hor, simple_tr, unresolved}`       | `No HOR`                |
| `hor_verdict = hor` & `hor_confidence < 0.20`          | `HOR weak`              |
| `hor_verdict = hor` & `0.20 ≤ hor_confidence < 0.40`   | `HOR moderate`          |
| `hor_verdict = hor` & `hor_confidence ≥ 0.40`          | `HOR strong`            |

Thresholds match `HOR_BIN_WEAK/MODERATE/STRONG` in current `kite.R` so
report colours, badges, and per-TRC roll-up counts stay numerically
comparable when graphs/tables are read across versions. If kitehor's
confidence calibration differs in practice from the R version, the
mapping table is the single tunable.

## Heatmap rendering (R, from kitehor's `--periodogram` bundle)

`tarean/kite_heatmaps.R` (new, ~140 lines, base R + `optparse`):

- Input: kitehor's `{prefix}_kite/kitehor.periodogram` (one
  `>case_id|H` line per array followed by a single dense-numeric
  line, alternating with `>case_id|bg` lines) plus the joined
  `monomer_size_top3_estimats.csv` for per-TRC top-3 peak overlay.
- Renders the same two PNGs per TRC the old `kite.R` produced:
  - `profile_plots/profile_{TRC}.png` — per-array `H[d]` heatmap
    (rows = arrays, columns = distance d, white→black intensity).
  - `profile_plots/profile_top3_{TRC}.png` — per-TRC summed profile
    on a log-x axis.
- Pure base R; no `numpy` / `matplotlib` dependency added.
  TideCluster's existing R stack (Biostrings, plyr/dplyr, scales,
  R2HTML, hwriter, etc.) already serves TAREAN, comparative analysis,
  and karyotype plots, so this just adds one more script alongside
  them.

Net change vs. 1.9.x: `tarean/kite.R` (736 lines, analysis + plotting +
HTML) → `tarean/kite_heatmaps.R` (~140 lines, plotting only) +
`tc_utils.build_kite_multifasta` / `build_monomer_size_csv` (~80
Python lines). All analysis logic moves out of TideCluster into the
external Rust binary.

## Integration points

### `TideCluster.py` (tarean step)

Replace lines 129–133:

```python
# RUN kite
cmd = (F"{script_path}/tarean/kite.R -d {prefix}_tarean/fasta"
       F" -p {prefix} -c {cpu}")
print("Making tandem repeat array profiles.")
tc.run_cmd(cmd)
```

with:

```python
# RUN kitehor (single analyze --periodogram call) + R heatmap renderer
kite_dir = F"{prefix}_kite"
os.makedirs(kite_dir, exist_ok=True)
multi_fa = F"{kite_dir}/_kite_input.fasta"
tc.build_kite_multifasta(F"{prefix}_tarean/fasta", multi_fa)
tc.run_cmd(F"kitehor analyze {multi_fa}"
           F" -o {kite_dir}/kitehor"
           F" --periodogram {kite_dir}/kitehor.periodogram"
           F" --threads {cpu}")
tc.build_monomer_size_csv(
    kite_tsv=F"{kite_dir}/kitehor.kite.tsv",
    summary_tsv=F"{kite_dir}/kitehor.summary.tsv",
    out_csv=F"{kite_dir}/monomer_size_top3_estimats.csv")
tc.run_cmd(F"{script_path}/tarean/kite_heatmaps.R"
           F" --periodogram {kite_dir}/kitehor.periodogram"
           F" --top3-csv {kite_dir}/monomer_size_top3_estimats.csv"
           F" --out-dir {kite_dir}/profile_plots")
os.remove(multi_fa)
```

`build_kite_multifasta()` and `build_monomer_size_csv()` are new
helpers in `tc_utils.py`. The former halves dimers; the latter joins
kitehor outputs and applies the hor_status binning mapping above.

The legacy v1 `_kite_report.html` is no longer produced by the
analysis step. To keep `_move_v1_to_legacy()` idempotent we either
(a) keep an HTML stub describing where the new data lives, or
(b) drop the legacy mover for the kite report (the function already
skips missing files). Option (b) is cleaner and matches the
removal-of-kite.R direction.

### `tc_rerender_report.py`

- Delete `_score_m_star`, `_hor_candidates`, `compute_hor` (lines 41–122).
- `load_kite_top3()` (line 673) updated to read lowercase
  `hor_status` / `hor_confidence` / `combined_class` / etc. columns
  directly (no re-derivation).
- Add `combined_class`, `subrepeat_flag`, `ssr_flag` to the per-TRC
  card and to the v2 `kite.html` aggregate.
- Backwards-compatibility for old runs (CSV with capitalised `HOR_*`
  columns and no `combined_class`): detect on column-name presence and
  upcase the headers in-memory before parsing. No separate legacy
  module — schema is close enough that one branch in `load_kite_top3`
  suffices.

### `tc_per_tra_consensus.py`

No change required. `monomer_size_top3_estimats.csv` keeps the same
primary columns; the new HOR / structural columns are appended at the
end, ignored by this consumer.

### `tarean/tarean_report.R` & per-TRC v1 HTML

The v1 per-TRC pages from `kite.R` (`trc_TRC_<N>.html`) are no longer
generated. The v2 report (`{prefix}_report_v2/trc/TRC_<N>.html`) is
already the documented landing page since 1.9.x; this just makes it
the only one. Update the v1 fallback message in `_move_v1_to_legacy`
accordingly.

### `batch_visualizer*.py` / `tidecluster_viz.*`

Quick grep showed no direct KITE references in these visualisers, but
sanity-check during implementation: if they consume any kite artefact,
update the schema reader.

## Conda + Singularity

`conda-deps.txt` — append:

```
kitehor=0.9.2          # add `-c petrnovak` to install command
```

CI / install instructions (`README.md`) update the `mamba create`
example:

```bash
mamba create -n tidecluster \
    -c conda-forge -c bioconda -c petrnovak \
    tidecluster kitehor
```

`TideCluster.def` (`%post`) — same channel addition for the SIF build.
This is the same `petrnovak` channel TideCluster itself is already
distributed from, so no new trust boundary.

A `mamba install` smoke check in the existing `tests/smoke.sh` should
add a `command -v kitehor && kitehor --help` probe.

If a bioconda recipe lands upstream later, dropping `-c petrnovak`
becomes a one-line change.

## Tests

| Test               | Change                                                                                          |
|--------------------|-------------------------------------------------------------------------------------------------|
| `tests/smoke.sh`   | Add `kitehor --help` probe in the CLI sanity section.                                           |
| `tests/short.sh`   | After `run_all`, assert `kitehor.summary.tsv` and `monomer_size_top3_estimats.csv` exist & have the new HOR + `combined_class` columns. |
| `tests/long.sh`    | Same as short, plus per-TRC heatmap PNG sanity (non-empty, reasonable size).                    |
| `tests/scripts/validate_kite_new_on_test3.R` | Already validates the R kite on a synthetic HOR ground-truth set; port the validation to read kitehor's `combined_class` / `hor_verdict` + `hor_multiplicity` columns and assert recall/exact-N against the test3 truth. |

Add a new regression test fixture if kitehor's classifier returns
substantively different categories on the canonical
`test_data/CEN6_ver_220406.fasta` carve — golden file checked in.

## Documentation

| File                                                      | Update                                                                                      |
|-----------------------------------------------------------|---------------------------------------------------------------------------------------------|
| `README.md`                                               | KITE section: "Powered by [kitehor](https://github.com/kavonrtep/kitehor)"; new output cols |
| `docs/hor_detection_design_notes.md`                      | Note migration to kitehor; legacy R algorithm history preserved                             |
| `docs/hor_classification.md` (referenced from `kite.R`)   | Replace with pointer to kitehor's classifier description                                    |
| `CLAUDE.md`                                               | Add kitehor to runtime deps section, note `-c petrnovak` channel                            |
| `changelog.md`                                            | New entry: "KITE step now uses kitehor; new per-array columns: combined_class, …"           |

## Phased implementation

1. **Plumbing.** Add `kitehor` to `conda-deps.txt`; `mamba install -n tidecluster -c petrnovak kitehor`. Verify `kitehor --version` runs in this sandbox. Write a 30-line standalone script that takes one `{prefix}_tarean/fasta/` dir and produces `monomer_size_top3_estimats.csv` of the new shape, to validate the input/output adaptation in isolation against `test_data/example/`.
2. **Helper move-out.** Extract the profile/heatmap code from `kite.R` into `tarean/kite_profiles.R`, with all empty-peaks guards preserved. Verify against `test_data/example/tc_kite/profile_plots/` (byte-identical PNGs ideal but not required; same dimensions and content).
3. **Wire into pipeline.** Implement `tc_utils.build_kite_multifasta`, `build_monomer_size_csv`, swap `TideCluster.py:129–133` block. Delete `tarean/kite.R`. Confirm `tests/smoke.sh` + `tests/short.sh` pass.
4. **Rerender consumer.** Remove `compute_hor` from `tc_rerender_report.py`; read HOR columns from the new CSV. Add `combined_class`, subrepeat, SSR cells to the v2 kite card. Add legacy fallback module for old outputs.
5. **Surface new outputs.** Expose `.windows.tsv` / `.verdicts.tsv` / `.hor_within_tile.tsv` in the v2 report's per-TRC dashboard (links + small table). Optionally a per-array `combined_class` chip beside the HOR badge.
6. **Validation.** Run `tests/scripts/validate_kite_new_on_test3.R` (ported) against `tmp/hor_test3/`; record recall / exact-N / FPR_null in `docs/hor_detector_v2_kite_defensive_guards.md` as a baseline for the new backend.
7. **Container + release.** Update `TideCluster.def`; rebuild SIF; bump `version.py` to `1.10.0` (breaks the kite output schema slightly → minor bump per existing convention); update `changelog.md`; update `README.md`.

## Open / deferred questions

- **kitehor bioconda recipe**: deferred — opportunistic contribution
  later; would let us drop `-c petrnovak`.

(The kitehor 0.9.2 glibc issue and the `--dump-profile` two-call
pattern are both resolved by upgrading to kitehor 0.9.3 + the new
`analyze --periodogram` flag.)
