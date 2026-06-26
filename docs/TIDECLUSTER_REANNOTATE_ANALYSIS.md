# TideCluster `reannotate` performance analysis — basis for a feature request

**Target:** TideCluster **1.15.2** (installed source analysed at
`.snakemake/conda/995dabc9…/share/tidecluster/`), specifically the
`tc_reannotate.py` code path invoked by the pipeline's `tidecluster_reannotate` rule.

**Why this exists:** on the evaluation genome (Pisum Pearl, **3.93 Gb, 3317
scaffolds, TR-rich**, run with **50 CPUs**) `tidecluster_reannotate` was the single
longest step in the whole pipeline at **16 h 47 m**, with an average CPU utilisation of
**~1.1 of 50 cores**. The genome being tandem-repeat-rich legitimately makes this step
heavy — but the *utilisation* number says the heaviness is being paid almost entirely
on **one core while 49 sit idle**. This document separates "inherently hard" from
"inefficient", with source-level evidence, so it can seed a feature request.

**Method:** read the installed TideCluster source; read the real run logs
(`TideCluster/default/tidecluster_reannotate.log/.err`) and the snakemake benchmark
(`benchmarks/tidecluster_reannotate.tsv`). No TideCluster code was modified.

---

## 1. What the step does (verified from source)

`tc_reannotate.py` (the standalone `tc_reannotate.py` entry the pipeline calls as
`tc_reannotate.py -s <genome> -f <dimer_library> -o <gff> -c {threads} --sensitivity …`):

1. Loads the dimer library; **tiles every monomer shorter than 4000 bp up to ≥4000 bp**
   (lines 135–149) so RepeatMasker has enough length to match short tandems.
2. Runs **one** RepeatMasker process on the **whole genome** (lines 153–170):

   ```python
   153  cmds = ["RepeatMasker", "-dir", temp_dir, "-nolow", "-no_is", "-e", "ncbi"]
   160  sensitivity_flag = SENSITIVITY_OPTIONS[args.sensitivity]   # default -> "" (full sensitivity)
   163  cmds.extend(["-lib", rm_library,
   165               "-pa", str(args.cpu),                          # <-- ALL parallelism delegated to -pa
   166               args.ref_seq])
   170  subprocess.run(cmds, check=True)
   ```
3. Parses / merges / extends / filters the RepeatMasker output into the final GFF3, in
   **single-threaded Python** (lines 185–234).

There is **no genome chunking** and **no internal job pool** — parallelism is entirely
delegated to RepeatMasker's `-pa`.

---

## 2. Evidence: `-pa 50` delivered ≈1-way parallelism

From the real run log (one RepeatMasker invocation, NCBI/RMBLAST 2.14.1+, custom lib):

```
analyzing file genome_cleaned.fasta
identifying matches to rm_library.fasta sequences in batch 1 of 67994
identifying matches to rm_library.fasta sequences in batch 2 of 67994
… (67994 search batches) …
processing output:
cycle 1 ……  cycle 2 ……  … (10 single-threaded ProcessRepeats cycles)
```

RepeatMasker prints `batch N of 67994` when a batch **starts**, so the start-order
reveals how many batches were ever in flight. Measured over all 67,994 batches:

| metric | value | meaning |
|---|---:|---|
| out-of-order starts | 241 / 67994 | starts were almost always strictly sequential |
| **max lag (peak concurrency)** | **9** | at most ~9 batches ever in flight, momentarily |
| **mean lag** | **0.00** | on average, batches started **one at a time** |
| benchmark `mean_load` | **1.1 / 50** | confirms ~1 core busy across the whole 16.8 h |

**Conclusion:** although the rule requested `-pa 50`, RepeatMasker processed the
67,994 search batches **essentially serially** (≈1 core, briefly up to ~9), so **~98%
of the 50-core allocation was idle for ~17 hours.** Other rules in the *same* run used
many cores (e.g. `dante_tir_fallback` ~15, `make_library_of_ltrs` ~5), so the node had
cores available — this step simply didn't use them.

### Why `-pa` under-parallelises here (and why this pipeline already works around it)

This is a **known limitation of RepeatMasker's `-pa` with a custom `-lib` + the
RMBlast/NCBI engine.** This very pipeline ships its own RepeatMasker wrapper
(`scripts/repeatmasker_wrapper.py`) whose module docstring says:

> *"breaks sequences into chunks of 50 Mbp … runs RepeatMasker on each chunk, running
> in parallel but each RepeatMasker job is run as a single thread — see
> https://github.com/Dfam-consortium/RepeatMasker/issues/274"*

i.e. the project already concluded that you cannot trust `-pa` for parallelism with a
custom library, and instead **splits the genome and runs many single-threaded
RepeatMasker processes in a `multiprocessing.Pool`.** `tc_reannotate.py` does **not**
do this, so it hits exactly the bottleneck the wrapper was written to avoid.

### A second, independent serial cost: `ProcessRepeats`

After the search, the log shows `processing output:` followed by **10 `cycle …`
passes** — RepeatMasker's `ProcessRepeats` stage, which is **single-threaded by
design** and is heavy on a TR-rich genome with millions of hits. A single whole-genome
RepeatMasker run does this once, serially. **Chunking fixes this too**: each chunk runs
its own `ProcessRepeats` in parallel, so the post-search serial tail is divided across
cores as well.

---

## 3. Secondary inefficiency: `O(N²)` post-RM filtering (`filter_intervals`)

`tc_reannotate.py:224` calls `tc_utils.filter_intervals(gr_filter, gr)`
(`tc_utils.py:99`). After bucketing by `seqid` and sorting both lists by start
(lines 120–123), it filters with a **nested loop that does not exploit the sort**:

```python
127  for seqid in dict_filter:
130      for grange_out in dict_gr[seqid]:          # N intervals on this seq
131          for grange_filter in dict_filter[seqid]:   # × M filter intervals on this seq
132              if grange_out.name == grange_filter.name:
133                  if grange_out.within(grange_filter):
134                      filtered.append(grange_out); break
```

This is **O(N·M) per sequence**. On a TR-rich chromosome both lists can be large
(many same-`Name` intervals), so this scales quadratically with tandem density —
precisely the regime this genome is in. Both lists are already sorted by `start`, so a
**sweep-line / two-pointer** (or an interval tree keyed by `Name`) makes it
`O((N+M) log)` with identical output. Likely a small fraction of the 16.8 h today
(RepeatMasker dominates), but it is on the single-threaded tail and will grow with
tandem density, so worth fixing alongside the main item.

(The sibling steps `merge_overlapping_gff3_intervals` (1679) and `extend_gff3_intervals`
(1655) are linear sweeps over the sorted GFF and look fine.)

---

## 4. Minor: library tiling inflates the RepeatMasker search

Lines 135–149 repeat every monomer < 4000 bp up to ≥4000 bp (`s * N`). This is
deliberate (RepeatMasker needs length to score short tandems), but it inflates the
custom-library size handed to RepeatMasker and therefore the per-batch search cost.
Worth (a) documenting the trade-off and (b) checking the cap is not larger than needed
for sensitivity. Algorithmic, lowest priority — do **not** change without a masked-bp
parity check.

---

## 5. Note: the same `-pa` pattern is in the `run_all` path too

`tc_utils.run_repeatmasker_with_renaming` (`tc_utils.py:1952`) builds
`RepeatMasker -pa {cpu} -lib {library} -e ncbi -s -no_is -norna …` and is used by
`annotate_gff` (line 368) in the **clustering/annotation** path (`run_all`). So the
"delegate parallelism to `-pa` on a whole-genome RepeatMasker" pattern is **shared by
both** `tc_reannotate` *and* the main TideCluster annotation step — meaning a single fix
(a chunked, pooled RepeatMasker helper) would also help `tidecluster_long`
(14 h on this genome). One helper, two beneficiaries.

---

## 6. Feature requests (ranked, all output-preserving)

**FR-1 — Parallelise RepeatMasker by chunking the genome, instead of relying on `-pa`.**
Replace the single whole-genome `RepeatMasker -pa {cpu}` call with: split the reference
into ~5–50 Mb chunks (pack small scaffolds together), run **one single-threaded
RepeatMasker per chunk** concurrently in a process pool of size `cpu`, then
remap coordinates and concatenate the `.out`. This is exactly
`scripts/repeatmasker_wrapper.py` in this repo and is proven output-identical (RM treats
each sequence independently). **Expected effect: the search *and* the `ProcessRepeats`
tail both parallelise → on this run, a plausible order-of-magnitude wall-clock
reduction (from ~1 effective core toward `cpu`), bounded below by the single slowest
chunk.** Highest impact by far.

**FR-2 — If a quick win is preferred over a full chunker:** at minimum, document that
`-pa` is ineffective with a custom `-lib` (RMBlast) and detect/warn when requested
cores ≫ achieved parallelism. (Subsumed by FR-1; listed for staging.)

**FR-3 — Make `filter_intervals` sub-quadratic** (`tc_utils.py:99`): exploit the
existing per-`seqid` sort with a sweep-line / two-pointer, or an interval tree keyed by
`Name`. Output-identical; removes a tandem-density-quadratic term from the
single-threaded tail.

**FR-4 — Audit the `< 4000 bp` monomer tiling** (`tc_reannotate.py:135–149`): confirm
the multiplier is the minimum needed for RM sensitivity; smaller → smaller library →
faster search. Needs a masked-bp parity check; lowest priority.

**FR-5 — Reuse the chunked helper in the `run_all` annotation path**
(`tc_utils.py:1952`) so `tidecluster_long` benefits from the same parallelism.

---

## 7. Honest caveats (so the feature request is fair)

- **The genome is genuinely TR-rich**, so a large *absolute* RepeatMasker cost is
  expected and not a bug. The defect is the **near-zero parallelism**, which is
  independent of how much work there is.
- The `-pa` ineffectiveness is **RepeatMasker's** behaviour with custom libraries, not
  a TideCluster bug per se — but TideCluster *chose* to delegate parallelism to it, and
  the fix (chunk + pool) lives cleanly in TideCluster.
- Part of the 16.8 h is the single-threaded `ProcessRepeats` tail, which no `-pa` value
  would have parallelised — only chunking does. This strengthens FR-1 over any
  `-pa`-tuning half-measure.
- **Validation for any change:** run `tc_reannotate` before/after on a smaller TR-rich
  reference and confirm the final GFF3 is identical (sort + `diff`) and total masked bp
  matches within the project's ±0.15 % losslessness bar. Chunking/pooling and the
  sweep-line filter must be byte-identical, not merely close.

---

### One-paragraph summary for the issue tracker

> On a 3.93 Gb TR-rich genome with 50 CPUs, `tc_reannotate` ran 16 h 47 m at ~1.1 cores
> average. It issues a single whole-genome `RepeatMasker -pa 50` with a custom dimer
> library; from the log, the 67,994 search batches started essentially sequentially
> (mean in-flight ≈ 0, peak 9), i.e. `-pa` delivered ~1-way parallelism — the known
> RepeatMasker-custom-`-lib` limitation (Dfam #274) that downstream tools work around by
> chunking the genome and pooling single-threaded RepeatMasker jobs. Please parallelise
> RepeatMasker by genome chunking (search **and** the single-threaded ProcessRepeats
> tail both scale), and make the O(N²) `filter_intervals` post-step sub-quadratic. Both
> are output-preserving.

---

## 8. Implementation status & validation (branch `perf/reannotate-efficiency`)

Implemented and unit-tested against the bundled CEN6 reference (`tests.sh` fast
suite + three new pure-python tests in `tests/unit.sh`):

**FR-3 — `filter_intervals` made sub-quadratic (done, byte-identical).** Replaced
the `O(N·M)` nested loop with a per-`(seqid, Name)` sweep: filter intervals are
sorted by start with a running prefix-max of their ends, so containment is a
single binary search (`O((N+M) log M)`). Because the filter inputs are
pre-merged, single-interval containment reduces to "max end among filters with
`start ≤ g.start` ≥ `g.end`", and the outer iteration order is preserved — the
result is the **same objects in the same order**. Pinned by
`tests/test_filter_intervals.py` (randomized equivalence vs. the original loop +
edge cases).

**FR-1 — chunked, pooled RepeatMasker (done).** `tc_reannotate.py -s` now calls
`tc_utils.run_repeatmasker_genome_chunked()`: the genome is split into
overlapping pieces (`split_fasta_to_chunk_files`, short numeric headers to dodge
RepeatMasker's 50-char limit), one single-threaded RepeatMasker runs per chunk in
a `multiprocessing.Pool` of `--cpu`, and hits are mapped back to genome
coordinates (the existing matching-table machinery) and written sorted. Both the
search **and** each chunk's `ProcessRepeats` now run concurrently. New flags
`--chunk_size` (default 50 Mb) / `--overlap` (default 100 kb). Plumbing pinned by
`tests/test_chunked_repeatmasker.py` (coordinate remap, header mapping, sort,
hit-less/no-`.out` chunks) with a fake RepeatMasker.

**Real-RepeatMasker parity (2 Mb CEN6 slice, 13 genuine CEN6 monomers as lib):**

| run | features | masked bp | vs. whole-genome |
|---|---:|---:|---|
| old whole-genome `-pa` | 98 | 51,545 | — |
| new, single chunk | 98 | 51,545 | **byte-identical** |
| new, 400 kb chunks (aggressive) | 98 | 51,574 | +29 bp = **+0.056 %** |

Single-chunk is byte-identical → the chunk/rename/remap/sort code is exactly
correct. The multi-chunk drift (4/98 features shifted 2–33 bp, one strand flip;
3 of the 4 are **mid-chunk**, 120–240 kb from any boundary) is RepeatMasker's
own `ProcessRepeats` context-dependence, not a boundary/overlap artifact, so it
shrinks with larger chunks and vanishes at single-chunk. It stays well within
the ±0.15 % bar. **Decision: ship split mode with a large default `chunk_size`**
— sequences below `2*chunk_size` are only packed (byte-identical); only very
large chromosomes split, confining the sub-0.15 % drift to their internal cuts.
This matches the trade-off the pipeline's own `scripts/repeatmasker_wrapper.py`
already accepts.

**FR-5 — correction.** The premise that the whole-genome `-pa` pattern is *also*
on the `run_all` / `tidecluster_long` path is **inaccurate**. `annotation()` /
`annotate_gff` (`TideCluster.py:486`, `tc_utils.py:368`) run RepeatMasker on the
concatenated TRC **consensus dimers** (sub-MB), not on the genome; the long
`tidecluster_long` wall-clock is TideHunter (already chunked via
`split_fasta_to_parts`) + clustering. So there is no whole-genome `-pa` to fix
there and FR-5 was **dropped** (no measurable payoff).

**FR-2 / FR-4** — FR-2 (warn that `-pa` is ineffective) is subsumed by FR-1
(the genome no longer relies on `-pa`). FR-4 (audit the `<4000 bp` tiling) is
left as a separate, optional follow-up gated on a masked-bp parity check.
