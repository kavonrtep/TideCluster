# TideCluster @ ~90 Gbp — performance & scaling analysis

Original analysis date: 2026-07-22. The body below is the original
bottleneck analysis; the fixes it proposes are tracked here.

## Implementation status (shipped in 1.18.0)

Output-identical (byte-for-byte on the deterministic outputs; verified per
fix and end-to-end on A. thaliana `run_all --long` vs the prior release):

- **Fix 1 — DONE.** `IndexedFasta` random-access replaces the whole-genome
  `fasta_to_dict` load in TAREAN / rDNA extraction.
- **Fix 2 — DONE.** O(1) `token_row` remap + per-header grouping in
  `split_fasta_to_chunks`.
- **Fix 3 — DONE.** `filter_gff_remove_duplicates` single linear pass.
- **Fix 5 — DONE, with a refinement not in the original text below:** the
  single-threaded worker pool runs **only when parts > cores**; when
  parts ≤ cores the original serial "-t {cpu} per part" loop is kept (all
  cores per part). Forcing `-t 1` unconditionally *regressed* small
  genomes (A. thaliana `--long`: 23.5 min vs 8.7). The pool is
  memory-gated (measured first-part peak RSS vs budget).
- **Fix 6 — DONE.** `tidehunter_long` streams per-round GFF3 to disk.
- **Fix 4 — CLOSED, not needed.** Current assemblies/tools cap sequence
  length at 2 Gbp, so no single scaffold exceeds 2.147 Gbp and the
  `rtracklayer` per-scaffold GRanges overflow cannot trigger. The only
  residual (cumulative genome-length `cumsum`/`sum`) is in the *untracked*
  `tidecluster_viz.R` / `plot_karyotype.R`, not the shipped pipeline; the
  single-genome report (`tarean_report.R`, `read.table` auto-promotes) has
  no overflow.
- **Fix 7 — OPEN (future work).** `tarean_report.R` per-TRC full-GFF
  rescans + `rbind` growth; a report-speed issue with many TRCs, not a
  blocker.

`tools/tidehunter_mem_probe.py` (per-round TideHunter peak-RSS probe) ships
alongside this note.

---

This note records where the pipeline breaks or degrades on an extremely
large genome (~90 Gbp — e.g. a giant plant/amphibian assembly), with a
focus on the *FASTA → chunk → reassemble → complete-analysis* path, and a
prioritized set of suggested fixes (§ "Suggested fixes 1–7").

The reference targets are `TideCluster.py`, `tc_utils.py`, `tc_reannotate.py`
and the R helpers in `tarean/` + the top-level `tc_*.R` scripts.

---

## 0. Verdict

As written, `run_all` on a 90 Gbp genome **will not finish**. It hits
several hard walls (OOM, quadratic hangs, R 32-bit integer overflow)
*before* mere slowdowns matter.

The **chunked-RepeatMasker reannotate path**
(`run_repeatmasker_genome_chunked`, `tc_utils.py:1147`) is already
engineered correctly for this scale (streaming, bounded file-descriptor
LRU, O(1) coordinate remap, parallel `Pool`). It is the **template** the
rest of the pipeline should be brought up to. The problems are
concentrated in three areas:

1. the **TideHunter chunk↔reassemble path**,
2. the **TAREAN whole-genome load**, and
3. the **R report / viz / comparative layer** (integer overflow + O(n²)).

### Reference scale at current defaults

Defaults: `chunk_size=500000`, `overlap=50000` (`TideCluster.py:792,940`);
`file_size_limit=50 MB` (`tc_utils.py:1327`).

| Quantity | Value at 90 Gbp |
|---|---|
| `matching_table` rows from `split_fasta_to_chunks` | **~180,000** |
| `run_tidehunter` TideHunter invocations (serial) | **~1,800** |
| TR arrays detected (satellite-rich genome) | **millions – tens of millions** |
| R 32-bit integer ceiling | **2,147,483,647** (~2.1 Gbp) — exceeded cumulatively with certainty, possibly per-scaffold |

---

## 1. The FASTA → chunk → reassemble path (core focus)

Path: `tidehunter()` (`TideCluster.py:781`) and the 3-round
`tidehunter_long()` (`TideCluster.py:919`).
Data flow: genome → one concatenated chunk temp-FASTA → TideHunter (split
again into 50-MB parts) → `.out` → coordinate-remap back to genome
coordinates → `<prefix>_tidehunter.gff3`.

### 1.1 🔴 BLOCKER — O(features × chunks) coordinate remap

`TideHunterFeature.recalculate_coordinates` (`tc_utils.py:1451`) calls
`get_original_header_and_coordinates` (`tc_utils.py:929`), whose body is a
**linear scan of the entire matching_table, per feature**:

```python
# tc_utils.py:941
matching_table_part = [x for x in matching_table if x[4] == new_header]
```

Called once per feature at `TideCluster.py:822` (`tidehunter`) and
`:852` (`parse_tidehunter_results_to_gff3`). At ~180k chunk rows ×
millions of features → **~10¹²⁺ operations** = multi-day hang.

**This exact bottleneck is already solved in the RepeatMasker path**
(`tc_utils.py:1232`), with a comment noting "there can be millions":

```python
# tc_utils.py:1232 — O(1) token -> row lookup
token_row = {row[4]: row for row in matching_table}
```

The fix was simply never applied to the TideHunter path. See **Fix 2**.

### 1.2 🟡 MODERATE — `run_tidehunter` runs ~1,800 parts serially

`run_tidehunter` (`tc_utils.py:1310`) splits the chunk file into ~50-MB
parts and iterates them in a plain Python loop with **no `Pool`**:

```python
# tc_utils.py:1339
for f in fasta_file_parts:
    ...
    subprocess.check_call(tidehunter_cmd, shell=True)   # serial, but -t cpu inside
```

**This is NOT a naïve single-threaded idle loop.** Each part *does* get
`-t {cpu}` (`TideCluster.py:799`, passed through at `tc_utils.py:1341`),
and TideHunter parallelizes across the ~100 sequences in a 50-MB part, so
cores are busy *within* a part. The throughput leaks are second-order,
each multiplied by ~1,800:

- **serial barrier** — a part finishes only when its slowest (most
  repeat-dense) 500 kb chunk finishes; other threads idle at the barrier
  while the next part waits;
- **tail underutilization** — the last wave of each part has fewer
  sequences than threads. At ~100 seqs/part this is minor at 32 cores
  (~3 waves) but significant at 96–128 cores (~1 wave), which is the
  regime a 90 Gbp HPC run actually lives in — so the damage **scales with
  core count**;
- **per-process overhead** — spawn / thread-pool setup / IO / teardown,
  ×1,800;
- **sublinear `-t`** — the code already distrusts single-process internal
  threading for the RepeatMasker engine (`run_repeatmasker_genome_chunked`
  docstring cites Dfam #274, "`-pa` does not parallelise effectively", and
  runs N single-threaded processes in a pool instead); the TideHunter
  path bets the opposite way.

Net effect: a **moderate** (~10–40 %, core-count-weighted) throughput
cost — **not** a blocker; a run that would otherwise finish just takes
longer. The serial design is also a **deliberate memory choice** (peak
RAM = one part; the docstring at `tc_utils.py:1315` says the whole split
exists because "tidehunter is consuming excessive amount of memory"), so
any parallelization must be memory-gated. `tidehunter_long` re-runs this
whole thing **3×** over the (masked) genome (`TideCluster.py:956`).
See **Fix 5**.

### 1.3 🟠 SEVERE — disk amplification

The path materializes the genome to temp disk repeatedly:

- gz-decompress to a full temp copy — `prepare_fasta_input`
  (`tc_utils.py:2706`, `shutil.copyfileobj`);
- ~99 GB chunk temp-FASTA — `split_fasta_to_chunks` (`tc_utils.py:881`);
- ~90 GB of TideHunter parts — `split_fasta_to_parts` (`tc_utils.py:1387`);
- a **full masked genome copy per round** — `mask_fasta_with_gff3`
  (`tc_utils.py:2659`, `bedtools maskfasta`).

Peak transient temp disk for `run_all --long` on gz input is
**several × genome size** and can exhaust `$TMPDIR`. See **Fix 5**.

### 1.4 🟡 MODERATE — repeated full-file scans + O(seqs × chunks) planning

`read_fasta_sequence_size` (`tc_utils.py:781`) reads **every base of the
90 GB file** just to sum sequence lengths, and is called independently by
`split_fasta_to_chunks`, `split_fasta_to_parts`, and `tarean`
(`TideCluster.py:64`) — several full 90 GB I/O passes, no `.fai` reuse.
Inside `split_fasta_to_chunks` the per-sequence planning scan
(`tc_utils.py:884`) is O(seqs × chunks). Peak RAM here scales with the
**largest single chromosome** (one whole sequence string resident at a
time, `tc_utils.py:883`): a multi-Gbp scaffold is a multi-GB Python `str`.

### 1.5 🟠 SEVERE — `tidehunter_long` holds every feature in RAM

`all_features_for_masking` and `gff3_lists` (`TideCluster.py:952-977`)
accumulate **every `TideHunterFeature` object across all 3 rounds
simultaneously** — millions of objects, each carrying a consensus string.
See **Fix 6**.

---

## 2. Hard blockers elsewhere (crash / OOM / hang)

### 2.1 🔴 Whole 90 Gbp genome loaded into a Python dict (TAREAN + rDNA)

`tarean()` calls `extract_sequences_from_gff3` (`TideCluster.py:60`) →
`gff3_to_fasta(..., load_sequence=True)` → `fasta_to_dict`
(`tc_utils.py:1488`), which reads the **entire genome into a dict of
strings** (~90–180 GB with Python overhead) → guaranteed OOM.

`clustering()` was *already* fixed for exactly this by passing
`load_sequence=False` (`TideCluster.py:626`, with an explicit OOM
comment). But TAREAN, and the rDNA genomic fallback
`_write_trc_region_subject` (`tc_utils.py:691`), still do the full load.
See **Fix 1**.

### 2.2 🔴 O(N²) `filter_gff_remove_duplicates`

Runs first in `clustering()` (`TideCluster.py:605`). It loads every
feature into `gff_data`, then inside the dedup loop rebuilds the entire
items list on every iteration:

```python
# tc_utils.py:2445 — inside `for i, (k1, v1) in enumerate(gff_data.items())`
k2, v2 = list(gff_data.items())[i + 1]   # O(N) rebuild per iteration → O(N²)
```

Over millions of arrays → hang. See **Fix 3**.

### 2.3 🔴 R 32-bit integer overflow (two distinct failures)

R integers are 32-bit (max 2,147,483,647).

- **Silent → NA.** `tidecluster_viz.R:75-76`:
  ```r
  contigs$offset_bp <- cumsum(c(0, contigs$length[-nrow(contigs)]))  # line 75
  genome_length    <- sum(contigs$length)                           # line 76
  ```
  `read.table` keeps `length` as *integer*; once the cumulative offset
  passes ~2.1 Gbp (certain at 90 Gbp), `cumsum`/`sum` overflow to `NA`
  (with an "integer overflow" warning). `genome_length` then flows into
  the report payload (`tidecluster_viz.R:209`) as `NA`.
- **Hard error.** `rtracklayer::import(.gff3)` builds GRanges with strict
  32-bit start/end (`tc_comparative_analysis.R:1098,1155,1590`;
  `tidecluster_viz.R:81`). **Any single scaffold > 2.147 Gbp aborts the
  import.**

See **Fix 4**.

---

## 3. Severe slowdowns / large memory (may finish on a very large box)

- **C1 — whole-genome feature structures held in RAM at once:**
  `merge_overlapping_gff3_intervals` (`tc_utils.py:2261`),
  `get_cluster_size2` (`tc_utils.py:2148`), and
  `add_cluster_info_to_gff3` (`tc_utils.py:2203`, a dict-of-dicts holding
  **every array's consensus + its ×2 dimer**). `resolve_trc_overlaps`
  (`tc_utils.py:2320`) additionally has an **O(features² per seqid)**
  covering test (`tc_utils.py:2367`). `run_repeatmasker_genome_chunked`
  accumulates **all RM hits genome-wide** in `records` before writing
  (`tc_utils.py:1234`).
- **C2 — clustering scales with array count (millions):** mmseqs
  `easy-cluster` + all-vs-all `blastn` on representatives
  (`-max_target_seqs 1000000`, `tc_utils.py:1805`) + a NetworkX graph of
  all edges (`get_connected_component_clusters`, `tc_utils.py:1655`).
- **C3 — R report O(n²):** `tarean_report.R` reads the whole-genome GFF
  and rebuilds it via **per-row `rbind`** (`list_to_dataframe`, lines
  28-44), then **re-scans the full GFF per TRC across ~8 columns**
  (lines 187-222). `plot_karyotype.R:52-92` and
  `tc_summarize_comparative_analysis.R:305-354` grow data frames by
  one-row `rbind` over all features + per-contig full scans (O(n²)).
- **C4 — comparative all-vs-all:** `tc_comparative_analysis.R:52-65`
  (MMseqs self-search, `--max-seqs = max(10000, N)`) and
  `compare_trc_by_blast.R:444-446` + dotplot self-BLASTs
  (`-num_alignments 10000000`), each followed by a full-pairwise
  `read.table` + igraph in memory. Bounded by TRC count (thousands)
  rather than genome bp, so less acute than §1/§2 unless the genome
  yields huge TRC counts. See **Fix 7** for the report half.

---

## 4. What already scales well (the template)

- **`run_repeatmasker_genome_chunked` / `split_fasta_to_chunk_files`**
  (`tc_utils.py:980-1271`, reannotate path): streams the genome once,
  **bounded 256-handle LRU** to avoid FD exhaustion (comment explicitly
  cites ~1,800 files at 90 Gbp, `tc_utils.py:1042`), **O(1) `token_row`**
  remap (`tc_utils.py:1232`), parallel `Pool` with serial library warm-up.
- `clustering()` consensus extraction with `load_sequence=False`
  (`TideCluster.py:626`) — the memory fix TAREAN/rDNA still need.
- Streaming stages: `.out`→GFF3, `filter_gff_by_length`,
  `add_attribute_to_gff`, TideHunter part concatenation.

---

## Suggested fixes 1–7

Ordered by (severity, then effort). Severity: 🔴 blocker (crashes /
hangs / OOMs), 🟠 severe (very slow or huge memory).

Each fix lists the location, the current behavior, why it breaks at
90 Gbp, a concrete approach, and correctness/risk notes.

---

### Fix 1 — 🔴 Stop loading the whole genome into RAM in TAREAN / rDNA

**Locations**
- `TideCluster.py:60` → `extract_sequences_from_gff3` (`tc_utils.py:1683`)
  → `gff3_to_fasta(..., load_sequence=True)` → `fasta_to_dict`
  (`tc_utils.py:1488`).
- rDNA genomic fallback `_write_trc_region_subject` (`tc_utils.py:676`,
  the `gff3_to_fasta(tmp_gff, fasta, "Name")` at line 691).

**Current behavior.** `fasta_to_dict` reads the *entire* FASTA into a
`{seqid: sequence_string}` dict. On 90 Gbp that is ~90–180 GB of RAM →
OOM. The sequence is needed here (unlike clustering) because TAREAN
extracts the actual genomic array sequences per TRC.

**Approach — random access via a FASTA index instead of a full dict.**
1. Build/reuse a `.fai` (e.g. `samtools faidx` or `pysam.FastaFile`,
   which builds the index on first open). This gives O(1) byte-offset
   seek to any `seqid:start-end` without holding the genome in RAM.
2. Rewrite `extract_sequences_from_gff3` to iterate the GFF3 once
   (already streamed) and, per feature, fetch just that interval:
   `fh.fetch(seqid, start-1, end)`. Peak RAM becomes O(largest single
   array), not O(genome).
3. For `_write_trc_region_subject`, same treatment — fetch each TRC array
   region on demand.

**Dependency note.** `pysam` is not currently in `conda-deps.txt`; a
zero-dependency alternative is to build the `.fai` (12-col fixed format:
name, length, offset, linebases, linewidth) once via a single streaming
pass and do the byte-offset seeks by hand (the offset math is the same
`.fai` spec bedtools/samtools use). Either way it also removes the
repeated full-file `read_fasta_sequence_size` scans (§1.4) if the `.fai`
is reused for lengths.

**Correctness/risk.** Output is byte-identical (same subsequences, same
order). Main risk is soft-masked/lowercase handling and line-wrapping in
the manual `.fai` route — `pysam` avoids that. Verify on the bundled
`test_data/CEN6_ver_220406.fasta` that extracted arrays match the current
`fasta_to_dict` path exactly (diff the per-TRC FASTAs).

**Effort:** medium.

---

### Fix 2 — 🔴 O(1) coordinate remap in the TideHunter path

**Locations**
- `get_original_header_and_coordinates` (`tc_utils.py:929`) and its caller
  `TideHunterFeature.recalculate_coordinates` (`tc_utils.py:1451`).
- Feature-loop call sites: `TideCluster.py:822`, `TideCluster.py:852`.

**Current behavior.** `[x for x in matching_table if x[4] == new_header]`
(`tc_utils.py:941`) is an O(len(matching_table)) scan **per feature**.

**Approach — reuse the RepeatMasker path's dict.** The RM path already
does this (`tc_utils.py:1232`):
```python
token_row = {row[4]: row for row in matching_table}   # new_header -> row
```
Build that dict once per round (in `tidehunter()` before the feature
loop, and in `parse_tidehunter_results_to_gff3`), then have
`recalculate_coordinates(matching_table, token_row)` do a single dict
lookup instead of the scan. Two clean shapes:
- pass `token_row` alongside `matching_table` into
  `recalculate_coordinates`, or
- precompute it inside `get_original_header_and_coordinates` **only if**
  the caller passes a prebuilt index (keep the old signature working for
  other callers).

The lookup arithmetic is unchanged:
```python
row = token_row[new_header]         # [orig_header, i, start, end, new_header]
real_chunk_size = row[3] - row[2]
ori_header      = row[0]
ori_start       = new_start + row[2]
ori_end         = new_end   + row[2]
```

**Correctness/risk.** Pure O(N·M) → O(N) speedup; results identical
(same row is selected — `new_header`/token is unique per chunk). Also
worth applying the same dict to the per-sequence scan at
`tc_utils.py:884` in `split_fasta_to_chunks` (group `matching_table` by
`row[0]` into a `dict[str, list]`).

**Effort:** low.

---

### Fix 3 — 🔴 De-quadratic `filter_gff_remove_duplicates`

**Location.** `tc_utils.py:2423`, called from `clustering()`
(`TideCluster.py:605`).

**Current behavior.** After building an ordered `gff_data`, the dedup loop
does:
```python
# tc_utils.py:2442-2447
for i, (k1, v1) in enumerate(gff_data.items()):
    if i == len(gff_data) - 1:
        break
    k2, v2 = list(gff_data.items())[i + 1]   # rebuilds the full list each iter
    if v1.seqid == v2.seqid and v1.start == v2.start and v1.end == v2.end:
        duplicated_ids.add(k1)
```
`list(gff_data.items())` is rebuilt on **every** iteration → O(N²) time
(and heavy allocation) over millions of arrays.

**Approach — single linear pass over the already-sorted items.** Materialize
the sorted items **once** and compare neighbors:
```python
items = list(gff_data.items())            # once
for (k1, v1), (k2, v2) in zip(items, items[1:]):
    if v1.seqid == v2.seqid and v1.start == v2.start and v1.end == v2.end:
        duplicated_ids.add(k1)
```
Semantics are identical (adjacent-in-sorted-order duplicate detection,
keeping the second). This is O(N).

**Memory follow-up (optional).** `gff_data` still holds every
`Gff3Feature` (with its consensus string) in RAM. For extreme scale the
dedup could instead stream the GFF3 twice keyed on `(seqid,start,end)`
without retaining feature objects — but the O(N²)→O(N) change alone
removes the hang; the memory reduction is a secondary concern.

**Correctness/risk.** Byte-identical output. Confirm with
`tools/founder_diff.py`-style before/after on a real run that the set of
dropped IDs is unchanged.

**Effort:** low.

---

### Fix 4 — 🔴 64-bit-safe genomic coordinates in the R layer

**Locations**
- Silent overflow: `tidecluster_viz.R:67-78` (`cumsum`/`sum` on integer
  contig lengths).
- Hard-fail import: `rtracklayer::import`/`import.gff3` at
  `tc_comparative_analysis.R:1098,1155,1590`; `tidecluster_viz.R:81`;
  plus `width(gff)` derivations (`tc_comparative_analysis.R:1156`).
- `read.table`-based GFF/fai reads that then do integer arithmetic:
  `tarean_report.R:20`, `plot_karyotype.R:44,66-67,85`,
  `tc_summarize_comparative_analysis.R:297,320-321,347`.

**Current behavior.** R's 32-bit integers overflow at ~2.147 Gbp.
`cumsum`/`sum` overflow *silently to NA*; GRanges construction from
coordinates > 2^31 *errors out*.

**Approach.**
1. **Force numeric (double) coordinate columns on read.** After each
   `read.table`, coerce length/offset/start/end columns with
   `as.numeric()` (doubles are exact for integers up to 2^53 ≈ 9.0e15,
   far beyond 90 Gbp). For `tidecluster_viz.R:75-76`:
   ```r
   contigs$length    <- as.numeric(contigs$length)
   contigs$offset_bp <- cumsum(c(0, contigs$length[-nrow(contigs)]))
   genome_length     <- sum(contigs$length)
   ```
   This removes the silent-NA path. `.fai` column 2 and GFF start/end
   should all be read as numeric.
2. **Avoid GRanges for whole-genome coordinates where a single scaffold
   can exceed 2.147 Gbp.** Options, cheapest first:
   - Do the interval math on plain numeric data frames (the code already
     mostly needs `width = end - start` and per-`Name` sums, which don't
     require GRanges).
   - If GRanges is genuinely needed, note it cannot hold >2^31 coords at
     all; the realistic mitigation is to keep coordinates in doubles and
     restrict GRanges use to per-array widths (small), never absolute
     positions on giant scaffolds.
3. **Audit `integer(0)` accumulators** (`plot_karyotype.R:66-67`,
   `tc_summarize:320-321`) and switch to `numeric(0)`.

**Correctness/risk.** For genomes under 2.1 Gbp the numeric coercion is a
no-op on values; downstream formatting may need `format(..., scientific=FALSE)`
to avoid `1e+09`-style labels. The GRanges avoidance is the larger change
and only strictly required if a *single sequence* exceeds 2.147 Gbp
(cumulative overflow is fixed by step 1 alone). Recommend gating: detect
`max(length) > 2^31-1` and warn/branch.

**Effort:** medium (step 1 is low; step 2 is the medium part).

---

### Fix 5 — 🟡 (throughput) / 🟠 (temp disk) — memory-gated parallel `run_tidehunter` + bound temp disk

Two separable concerns are bundled here: **5a** parallelism (moderate,
memory-governed) and **5b** temp-disk amplification (severe — can abort
the run by exhausting `$TMPDIR`).

**Locations.** `run_tidehunter` (`tc_utils.py:1310`, serial loop at 1339);
disk amplification across `prepare_fasta_input` (`tc_utils.py:2706`),
`split_fasta_to_chunks` (`tc_utils.py:881`), `split_fasta_to_parts`
(`tc_utils.py:1387`), `mask_fasta_with_gff3` (`tc_utils.py:2659`).

#### 5a — Parallelize parts, **memory-gated** (throughput)

**Current behavior.** ~1,800 parts run one-at-a-time; each uses
TideHunter's `-t {cpu}` internally. Peak RAM = one part's multi-threaded
TideHunter working set.

**Approach.** A `Pool` of module-level picklable workers, each running one
**single-threaded** TideHunter on one part (mirror
`run_repeatmasker_genome_chunked`, `tc_utils.py:1200`). Search is
embarrassingly parallel across parts, so pool_size workers × `-t 1` gives
cleaner near-linear scaling than one process × `-t cpu`.

> **⚠ Memory is the governing constraint — watch this.** The parallel
> design can use **more** RAM than today's serial one for the *same*
> concurrency:
>
> | | concurrent seq working sets | address spaces |
> |---|---|---|
> | current (serial, `-t cpu`) | ~cpu | **1 shared** |
> | parallel (pool=cpu, `-t 1`) | ~cpu | **cpu separate** |
>
> Same number of sequences in flight, but the parallel version pays
> per-process duplication (base heap, input buffer, output buffers, any
> loaded tables) ×pool_size instead of sharing one address space. So a
> naïve `pool_size = cpu` can raise peak RAM and OOM where the serial
> path fit.

Controls (in priority order):

1. **Gate `pool_size` on a memory budget, not on `cpu`:**
   `pool_size = max(1, min(cpu, floor(budget / part_peak)))`, with
   `budget` from `$AGENT_MEMORY` minus headroom. Never hardcode `cpu`.
2. **Measure `part_peak` empirically — don't guess.** Use
   `tools/tidehunter_mem_probe.py` (added with this doc; see §5c for
   results). The sandbox has no `/usr/bin/time`; the reliable peak is
   `os.wait4(pid, 0)` → `rusage.ru_maxrss` (kernel-exact, no sampling),
   which the probe confirms agrees with TideHunter's own `Peak RSS`
   line to < 1 % and is reproducible across repeats. **Prefer it over a
   `/proc/<pid>/status:VmHWM` poller** — an early Popen + 10 ms-VmHWM-poll
   harness read 296 MB where TideHunter self-reported 535 MB on the same
   (first, cold) run; `ru_maxrss` avoids that class of disagreement
   entirely. **Log the measured per-part peak and the chosen pool_size**
   so memory use is observable during the run (the "watch memory"
   requirement).
3. **Shrink part size when running in parallel — but it only helps
   rounds 1–2.** The 50-MB `file_size_limit` (`tc_utils.py:1327`) was
   sized for the serial path. Smaller parts cut each worker's footprint
   (so `pool_size × part_peak` fits the budget) **and** improve load
   balance (finer granularity → less tail waste, §1.2). **Caveat from
   §5c:** round 3 (long periods) is *base-dominated* (~1.9 GB fixed
   regardless of part size), so shrinking parts barely lowers its
   footprint — rounds 1–2 do scale with size and benefit. Trade-off:
   more parts → more per-process overhead; tune once `part_peak` is
   known per round.
4. **Dynamic admission (optional, heavier).** Launch a new worker only
   when free system memory is above a watermark — robust to a
   repeat-dense part that spikes, but more complex than a static
   budget-derived `pool_size`. Consider only if step 3 can't keep peaks
   bounded (some satellite-dense 500 kb chunks can spike well above the
   median part).

> Note: the memory watch applies to the **current serial path too** — at
> very high `cpu`, `-t cpu` on a repeat-dense part can already spike; the
> 50-MB split bounds the *input* per part, not the *working set*.

**Correctness/risk.** Output-preserving (each part is an independent set
of sequences; the chunk-overlap guarantees are unchanged; merged output is
regenerated downstream so order doesn't matter). Validate feature counts
are identical serial-vs-parallel on `test_data`, and record peak RSS at a
couple of `cpu`/part-size settings before choosing defaults.

#### 5b — Bound peak temp disk (severe; can abort the run)

`run_all --long` on gz input can hold several × genome size in `$TMPDIR`
simultaneously (§1.3). Bound it:

1. **Stream-delete each TideHunter part** right after it is consumed into
   the merged `.out`, instead of only removing the parts dir at the very
   end (`tc_utils.py:1357`, currently after *all* parts exist).
2. **For `--long`, delete each round's masked FASTA** as soon as the next
   round's chunk file is written (currently removed later at
   `TideCluster.py:981`).
3. **Skip the intermediate concatenated chunk temp-FASTA** by having
   `split_fasta_to_chunks` write directly into part files (merge it with
   `split_fasta_to_parts`), halving the ~2× genome temp footprint.
4. **Reuse the `.fai`** from Fix 1 so `split_fasta_to_parts` /
   `read_fasta_sequence_size` don't re-scan 90 GB from disk.

#### 5c — Measured per-part footprint (empirical, drives 5a)

Measured with `tools/tidehunter_mem_probe.py` on
`test_data/CEN6_ver_220406.fasta` (180 Mb centromeric satellite — dense,
worst-case-ish), **single-threaded** (`-t 1` = the per-worker footprint the
pool formula needs). Peak RSS from `os.wait4().ru_maxrss` (kernel-exact;
cross-checked against VmHWM and TideHunter's own `Peak RSS` line — all
three agree to <1 %). Each part is a pack of 500 kb chunks, exactly as
`split_fasta_to_chunks` produces.

| round (period range) | 5 Mb part | 15 Mb part | 30 Mb part | wall @30 Mb | rough fit (sublinear) |
|---|---|---|---|---|---|
| **1** · `-p 40 -P 3000`      | 296 MB  | 590 MB  | 833 MB  | 48 s  | ~180 MB base + ~22 MB/Mbp |
| **2** · `-p 3001 -P 10000`   | 686 MB  | 1363 MB | 1620 MB | 112 s | ~480 MB base + ~40 MB/Mbp |
| **3** · `-p 10001 -P 25000`  | 2208 MB | 2791 MB | 3023 MB | 255 s | **~1.9 GB base** + ~33 MB/Mbp |

(Slopes fall as part size grows — all three rounds are sublinear in input
size; the "base" is the y-intercept of the two-endpoint fit.)

(Single-round `tidehunter` uses round-1 args by default,
`TideCluster.py:1107`, so its footprint = the round-1 row.)

Two findings reshape the 5a plan:

1. **The three rounds differ ~7–10× in peak** at equal input (296 MB vs
   2208 MB for a 5 Mb part). The pool MUST be sized on the **worst round
   (round 3)** — sizing on round 1 / the default args would admit ~7× too
   many workers and OOM the instant round 3 starts.
2. **Round 3 is base-dominated (~1.9 GB fixed).** Its peak barely grows
   with part size (5→30 Mb: 2208→3023 MB), so 5a-control-#3
   ("shrink parts") does **not** help round 3 — even a 1 Mb part costs
   ~1.9 GB/worker. Rounds 1–2 *do* scale with size, so shrinking helps
   them.

**Design consequence — size the pool _per round_, not once.** Because
`tidehunter_long` runs the rounds sequentially (`TideCluster.py:956`), each
round can take its own pool size from that round's measured peak:

```python
budget = AGENT_MEMORY_MB * 0.8                 # leave headroom
pool[r] = max(1, min(cpu, budget // peak_mb[r]))
```

At a 64 GB budget that is **~200 round-1 workers but only ~28 round-3
workers** (`tools/tidehunter_mem_probe.py` prints exactly this hint).
Round 3 is also the **slowest** (255 s vs 48 s for a 30 Mb part) while
finding the fewest arrays on this data (long-period satellites are rare) —
so it is simultaneously the memory *and* throughput bottleneck of the long
path, and the narrow round-3 pool is the binding constraint on the whole
parallelization. A single genome-wide pool sized on round 1 is the trap to
avoid.

**Reproducibility & safety margin.** With the clean `os.wait4` measurement,
warm runs are **essentially deterministic**: 3 repeats gave
round-1/5 Mb = 295.4 / 295.5 / 295.6 MB and
round-3/5 Mb = 2208.4 / 2208.2 / 2208.1 MB (spread < 0.1 %). The one
outlier — the very first, cold invocation self-reported 535 MB for
round-1/5 Mb vs ~289 MB warm — is a *single confounded data point*
(first-run/cold-cache **and** the discarded poller harness), so it isn't
proof of run-to-run variance; but it is a reason to keep headroom rather
than size to the last megabyte. The margins that actually matter:

- **Data variability** — different parts of a real genome carry different
  satellite content; a part with longer / more numerous long-period arrays
  can peak higher than this centromeric fixture's round-3 base. Measure on
  a *sample of parts from the actual target genome* and size on the worst.
- **Cold start / machine / allocator** — the first invocation, THP, glibc
  arena count, and total RAM differ across hosts; keep the `× 0.8` budget
  headroom (and consider `MALLOC_ARENA_MAX` to cap per-thread arena bloat
  if it shows up).

Net: size the per-round pool on the **worst measured part for that round on
the real genome**, with the 0.8 headroom — not on a single fixture
observation.

**Effort:** medium (5a is the memory-tuning work + per-round pool sizing;
5b is mostly delete-earlier plumbing).

---

### Fix 6 — 🟠 Stream `tidehunter_long` features to disk instead of RAM

**Location.** `tidehunter_long` (`TideCluster.py:919`); the RAM
accumulators `all_features_for_masking`, `gff3_lists`
(`TideCluster.py:952-977`).

**Current behavior.** Every `TideHunterFeature` from all 3 rounds is held
in Python lists simultaneously (millions of objects, each with a consensus
string) to (a) build the masking GFF3 for the next round and (b) write the
final merged GFF3.

**Approach.**
1. Have each round **write its features straight to a per-round GFF3 file**
   (already done for `keep_rounds`, `TideCluster.py:906`) and drop the
   in-memory list.
2. Build the next round's masking input by **concatenating the prior
   rounds' GFF3 files on disk** (a streamed `cat`) rather than
   `save_gff3_to_file(all_features_for_masking, ...)`
   (`TideCluster.py:965`).
3. Produce the final merged output by streaming-concatenating the 3
   per-round GFF3 files (`TideCluster.py:985`).

Peak RAM for the round loop becomes O(one round's parse buffer) instead of
O(all features across all rounds). Note `parse_tidehunter_results_to_gff3`
(`TideCluster.py:833`) itself still builds a per-round list — for full
streaming it could yield features and write as it goes.

**Correctness/risk.** Output GFF3 is identical (same features, merge is
just concatenation; downstream `merge_overlapping_gff3_intervals` already
handles ordering/dedup). Masking input is identical (same union of
regions). Low risk.

**Effort:** low.

---

### Fix 7 — 🟠 De-quadratic the R report / karyotype prep

**Locations.**
- `tarean_report.R`: `read_gff3` + `list_to_dataframe` (lines 18-44,
  per-row `rbind`); the per-TRC full-GFF rescan block (lines 187-222).
- `plot_karyotype.R:52-92` and
  `tc_summarize_comparative_analysis.R:305-354`: one-row `rbind` growth
  over all features + per-contig full scans.

**Current behavior.**
- `list_to_dataframe` does `do.call(rbind, lapply(input_list, as.data.frame))`
  over millions of 1-row frames → effectively O(n²) copies.
- The main TR summary computes ~8 columns each as
  `sapply(summary_df$TRC, function(x) ... gff[gff$attributes$Name == x, ] ...)`
  — a full logical scan of all arrays **per TRC** → O(n_TRC × n_arrays).
- `plot_karyotype.R` / `tc_summarize` grow a data frame by `rbind` one row
  at a time inside a `for` over all features → O(n²) reallocation, plus
  per-contig `max(end[seqname==contig])` full scans.

**Approach.**
1. **Parse GFF attributes vectorially.** Replace the per-row
   `process_attributes` + `rbind` with a vectorized `data.table::fread`
   (or `vroom`) read and a single `tstrsplit`/`regmatches` extraction of
   the needed attributes (`Name`, `ssr`, `annotation`, …) into columns.
   Avoid the nested-data.frame `gff$attributes` shape entirely.
2. **Group once, not per-TRC.** Replace the ~8 `sapply(TRC, ... subset ...)`
   columns with a single `split(gff, Name)` (the code already does this in
   `extract_summary_table`, `tarean_report.R:306`) or a `data.table`
   `by = Name` aggregation computing min/max/median array length, total
   size, type, SSR, annotation in one pass → O(n).
3. **Preallocate / vectorize the karyotype prep.** Replace the per-row
   `rbind` loops with a vectorized build: read the GFF as a data frame,
   then `aggregate`/`data.table` `by = seqname` for per-contig max end;
   build `contigs` with a single `data.frame(...)` call, not incremental
   `rbind`.

**Correctness/risk.** Numeric outputs unchanged; verify the rendered
report tables match on a real run (e.g. re-render an existing run dir and
diff the summary TSV/HTML tables). `data.table`/`vroom` are new R deps —
check whether the base-R `split()` + `vapply` route is enough before
adding a dependency (it removes the O(n²) `rbind` and the per-TRC rescan
without new packages). Combine with **Fix 4** (numeric coordinates) so the
same read path is 64-bit-safe.

**Effort:** medium.

---

## Priority summary

| # | Fix | Severity | Effort |
|---|-----|----------|--------|
| 1 | TAREAN/rDNA: `.fai` random access instead of `fasta_to_dict` full load | 🔴 blocker | med |
| 2 | TideHunter remap: O(1) `token_row` dict (as RM path already does) | 🔴 blocker | low |
| 3 | `filter_gff_remove_duplicates`: single linear pass, drop `list(...)[i+1]` | 🔴 blocker | low |
| 4 | R: numeric (double) coordinates; avoid GRanges >2.1 Gbp | 🔴 blocker | med |
| 5 | Memory-gated parallel `run_tidehunter` parts (🟡) + cap peak temp disk (🟠) | 🟡 / 🟠 | med |
| 6 | Stream `tidehunter_long` features to disk instead of RAM | 🟠 severe | low |
| 7 | R report/karyotype: vectorize GFF parse + group-once (no per-TRC rescan) | 🟠 severe | med |

Suggested landing order: **2 → 3 → 6** (small, high-value Python wins) →
**1** (the main OOM) → **4** (R overflow) → **5** (throughput) → **7**
(report). Each is independent; 1 and 5 share the `.fai` reuse.

### Cross-cutting: reuse the reannotate template

`run_repeatmasker_genome_chunked` (`tc_utils.py:1147`) already embodies the
right patterns for 90 Gbp — streaming single pass, bounded 256-handle LRU
(`tc_utils.py:1048`), O(1) `token_row` (`tc_utils.py:1232`), parallel
`Pool` with serial warm-up. Fixes 1, 2, 5 are largely "make the TideHunter
and TAREAN paths look like that one."
