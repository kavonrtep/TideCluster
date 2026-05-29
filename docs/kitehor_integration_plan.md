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

## Addendum: kitehor 0.10.0 (shipped in TideCluster 1.10.0)

The integration was first built against kitehor 0.9.3 (8-class
schema). The 1.10.0 release pins **kitehor 0.10.0**, which is a
breaking change to the classifier and the `summary.tsv` schema. The
`analyze --periodogram` invocation is unchanged; only the joined CSV
and the report consumers were updated.

**Categories (7, was 8).** `tr_with_nested_tr` is removed (those
arrays now fall to `tr`, `tr_with_subrepeat`, or `unresolved`);
`hor_with_ssr` is added. Final set: `hor`, `hor_with_ssr`, `tr`,
`tr_with_ssr`, `tr_with_subrepeat`, `pure_ssr`, `unresolved`. Decision
rules (per `kavonrtep/kitehor` `docs/rule_proto.md`): `pure_ssr` when
`ssr_dominant_motif_coverage_pct ≥ 80`; `tr_with_subrepeat` when
`tv_decision = localized_subrepeat`; `hor_with_ssr` / `hor` from
`hor_verdict=hor` (± `ssr_flag`); `tr_with_ssr` / `tr` from
`hor_verdict=simple_tr` (± `ssr_flag`); else `unresolved`.

**`summary.tsv` (still 32 cols, reshuffled).** Dropped: `length_bp`,
all `subrepeat_*`, `density_hint`, `founder_density`, `phase_contrast`,
`density_n_windows`. Added 10 `tandem_validate` columns: `tv_decision`,
`tv_host_period`, `tv_best_candidate_period`, `tv_best_candidate_kind`,
`tv_density`, `tv_spatial_contrast`, `tv_phase_contrast`,
`tv_n_windows_total`, `tv_n_windows_present`, `tv_reason`. The
subrepeat signal is now: host monomer = `tv_host_period`, subrepeat
period = `tv_best_candidate_period`. `length_bp` is recovered by
joining on `record_id` against `.kite.tsv` (`array_length`).

**`analyze` outputs (6 TSVs, was 9).** `.tandem_validate.tsv` replaces
`.subrepeat.tsv` + `.hor_within_tile.tsv`; `.windows.tsv` is gone.
Removed subcommands: `subrepeat-scan`, `hor-validate`.

**TideCluster code touched:** `tc_utils.build_monomer_size_csv`
(lift `tv_*` instead of `subrepeat_*`/`founder_density`/`phase_contrast`),
`tc_rerender_report.py` (class set + palette + `structure_cell` +
`load_kite_top3`), `tarean/assets/tidecluster.css` (add
`hor_with_ssr`). The 0.9.x `subrepeat_*` / `tr_with_nested_tr` paths
are kept as legacy fallbacks so old output dirs still rerender.

## Addendum: kitehor 0.12.0 (TideCluster 1.10.0 shipped)

After the 0.10.0 cut a second breaking change landed: the
`combined_class` cascade was retired in TideCluster in favour of a
leaner two-stage flow built directly on `kitehor rescore`. **kitehor
0.12.0 is the pinned conda dep** in `conda-deps.txt`; v0.11.0 added
the `hor_with_ssr` / `unresolved_with_ssr` classes and switched
`pure_ssr` to read `ssr_raw_total_coverage_pct`, but TideCluster
never shipped with it.

**Pipeline change.** TideCluster.py runs
`kitehor kite-periodicity --periodogram --out-peaks` first, then in
parallel `kitehor rescore --max-period 10000 --top-n 20`
and `kitehor ssr-scan`. No `analyze`, no `rule-classify`, no
`tandem-validate`, no `summary-merge`. Two new CLI flags expose the
rescore caps: `--kite-rescore-max-period` and `--kite-rescore-top-n`.

**Founder reassignment.** rescore's `founder_period` column is the
peak with the highest `identity_med` (gated at ≥ 0.7). On a real HOR
record that's typically the *HOR tile* (e.g. 1512 bp on CEN6's
classic 503/1512/×3 case), not the base monomer. TideCluster treats
rescore's `founder_period` as the **strongest** period and reassigns
**founder** to the smallest peak P such that
`strongest ≈ k·P` (k integer, 2 ≤ k ≤ 30, ±0.05 tolerance) with
`identity_med(P) ≥ 0.7`. When no qualifying divisor exists,
founder = strongest, multiplicity = 1. When rescore returns NA
`founder_period` (no peak passes 0.7, or every peak above
`--max-period`), TideCluster falls back to the kite rank-1 peak and
flags the array as a fallback.

**Subrepeat tiering.** Per-peak classifier ported from kitehor's
`docs/rule_proto.md` cheat-sheet:

| tier | condition (period ≤ founder/4 unless noted) |
|---|---|
| `HIGH` | `scan_occupancy_frac ≥ 0.15` AND (`subrepeat=true` OR `phaseC ≥ 0.10` OR `autoF ≥ 0.4`) |
| `LIKELY` | `scan_occupancy_frac ≥ 0.20` AND `scan_n_intervals ≥ 10` |
| `KMER_SUPPORT` | `phaseC ≥ 0.10` OR `autoF ≥ 0.4` (no per-base scan support) |
| `AMBIGUOUS` | `0.25 < ratio ≤ 0.33` with mild support |
| `WEAK` | period ≤ founder/4 but no signal fires |
| `OBSERVATIONAL` | founder NA but `scan_occupancy_frac ≥ 0.05` |
| `REJECT_*` | phantom / `period > founder/3` / `scan_occ = 0` & founder known / founder-itself / weak ambiguous |

Validated on the curated `test_data/IPIP200579_2026-04-14` corpus:
40 / 3024 records carry HIGH+LIKELY (~1.3 %), including the
cheat-sheet canonical TRC_104 (founder 180, sub 36, occ ~0.5) and
TRC_666 (founder 250, sub 36, occ 0.46). Subrepeats are
intentionally rare: a real nested subrepeat is a *short motif tiling
inside the founder monomer but occupying only part of it* — full
occupancy with integer multiplicity is HOR, not subrepeat.

**Report (v2) — per-TRA table** is now:
`▸ # · seqid · start · end · length · Founder · Δid · Strongest · ×k
 · Other periods · Subrepeat · SSR`. Coloured tier pills mark
subrepeat strength; the `▸` control opens a Details child row
(DataTables `row().child()` + a small handler in `tidecluster.js`)
showing the founder/strongest summary, the SSR scan, and a full
table of every rescored peak with its tier and diagnostics. The
`combined_class` driven roll-ups (per-TRC Classification card,
index "Array classes", KITE overview class bar / class mix column,
ideogram colours) are replaced by per-array counters: HOR (×k≥2),
Subrepeat, SSR, Fallback founder. Genome-distribution ideograms
colour by the dominant signal per array.

**`monomer_size_top3_estimats.csv` schema (v0.12-derived).** New
columns: `founder_period`, `strongest_period`, `multiplicity`,
`delta_id_pp`, `founder_id_med`, `strongest_id_med`,
`founder_fallback`, `subrepeat_{1,2}_period`/`occ`/`tier`. Empty
`hor_status` / `hor_confidence` cells are kept for the per-TRA
consensus R script's `HOR_*` alias shim. The legacy 0.10.x (`tv_*`,
`combined_class`) and 0.9.x (`subrepeat_*`, `tr_with_nested_tr`)
columns are still consumed by `tc_rerender_report.py` when present
so old output dirs still rerender.

**TideCluster code touched:** `TideCluster.py` (replace `analyze`
call, expose two CLI flags, run rescore + ssr-scan as parallel
subprocesses), `tc_utils.build_monomer_size_csv` (read rescored
peaks + ssr.tsv, reassign founder, tier subrepeats, emit the new
schema), `tc_rerender_report.py` (load the new CSV + the rescored
peaks TSV; rewrite `_arrays_table` + add `_details_html` /
`_other_periods_cell` / `_subrepeat_cell` / `_ssr_cell`; drop the
class-driven roll-ups), `tarean/assets/tidecluster.js`
(`initDetailsExpand` for the DataTables child rows; generalise the
floating tooltip), `tarean/assets/tidecluster.css` (tier pills,
details row, fallback marker).

## Addendum (cont.): TRC-aware two-pass founder reassignment

After the initial v0.12 cut shipped, we extended the per-array
founder reassignment into a **three-pass** algorithm that exploits
cluster-level context. The motivating case is CEN6's TRC_1:
`84848657..84876251` — strongest period 2386 bp (id_med 0.998), a
clear ~503 bp peak at id_med 0.736 (passes the 0.7 gate), but
2386/503 = 4.744 is well outside the ±0.05 strict integer-multiple
window. The cluster as a whole prefers a ~500 bp founder (median 502
across 31/33 arrays once we look TRC-wide).

**Algorithm (`tc_utils.build_monomer_size_csv`):**

1. **Pass 1 — per-TRA strict**: smallest divisor with `id_med ≥ 0.7`,
   integer k ∈ [2, 30], ratio tolerance ±0.05 (unchanged).
2. **Pass 2 — TRC consensus rescue**: per TRC, the consensus founder
   is the median of the largest ±10 % cluster of Pass-1 founders;
   gated by ≥ 3 arrays AND ≥ 25 % of the TRC's reassigned arrays.
   Arrays where Pass 1 left founder = strongest look for a rescored
   peak within `max(±10 %, ±20 bp)` of the consensus with
   `id_med ≥ 0.5` (relaxed because the cluster context provides a
   strong prior). The rescue is only adopted when the resulting
   multiplicity rounds to ≥ 2 (rescues where strongest ≈ founder add
   no HOR signal). The array gets `irregular_multiplicity = true`
   and the raw fractional k is stored in `multiplicity_raw`.
3. **Pass 3 — solo / no-consensus relaxed individual**: per-array
   tolerance relaxed to ±0.20 of integer, still `id_med ≥ 0.7`,
   k ∈ [2, 30]. Flagged irregular when it fires. Designed for TRCs
   with no consensus (e.g. solo-TRA TRCs).

**Validated on CEN6**:
- TRC_1 strict-only: 7/33 arrays reassigned.
- After Pass 2: 18/33 arrays carry an HOR multiplicity (11 of them
  rescued via the consensus, 7 strict), with the distribution
  `×2 (12), ×3 (2), ×4 (1), ×5 (3)`.
- Prevalent founder for TRC_1: 502 bp (in 31/33 arrays).
- Across the whole assembly: HOR call count rose from 7 (strict-only)
  to 19 (7 strict + 12 rescued).

**Report:**

- Per-TRA table `×k` cell appends a `~` pill (amber) when
  `irregular_multiplicity = true`; hovering shows the raw fractional
  k (e.g. `~` → "irregular multiplicity (k = 4.74)").
- Details child row header gains an `irreg` tag plus `(k = 4.744)`.
- Per-TRC Classification card adds **Prevalent founder**
  ("502 bp (in 31/33 arrays)"), **HOR multiplicities** distribution
  (`×2 (12), ×3 (2), ×4 (1), ×5 (3)`), and an **Arrays with
  irregular multiplicity** counter when non-zero.

**CSV schema additions** (`monomer_size_top3_estimats.csv`):
`multiplicity_raw`, `irregular_multiplicity`.
