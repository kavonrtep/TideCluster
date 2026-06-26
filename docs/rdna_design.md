# rDNA identification — design & implementation notes

Implements the request in [`rdna_feature`](rdna_feature): flag which TRCs are
ribosomal DNA (45S vs 5S) by similarity to a reference library, surface it
cleanly in the clustering GFF3 and the HTML report.

## Decisions (agreed with maintainer)

- **No internal structure.** We label a TRC `45S` or `5S` only; we do not
  annotate 18S/ITS/5.8S/IGS/25S sub-structure inside the array.
- **Metric = best single-subunit reference coverage.** A TRC is rDNA if *one*
  library reference (e.g. an 18S or 25S gene) is matched near-full-length at
  high identity. This is robust to two facts confirmed on the test data:
  the library contains only the genes (so ITS/IGS spacer dilutes whole-monomer
  coverage), and matches are cross-species (~90 % identity), so the identity
  gate is 85 %, not 95 %.
- **Engine = blastn.** Reference coverage is what blastn yields directly (the
  reference is the query, so query-coverage == reference-coverage), using the
  same blast tooling the superfamily step already relies on.
- **Query tiers (consensus-first, genome fallback).** (1) TAREAN dimer library
  where available, else (2) the raw TideHunter clustering consensus (present for
  every TRC), and only (3) the TRC's genomic arrays for TRCs with no usable
  consensus.
- **Output = clean call + evidence.** `rDNA_type=45S|5S;rDNA_coverage=<frac>` in
  the clustering GFF3 (and the annotation GFF3 if present). Subunit detail is
  not written to the GFF3; the per-TRC table is in `<prefix>_rdna.tsv`.
- **Default-on** in `run_all`/`tarean` (disable with `--no_rdna`); bundled
  library at `data/rdna_library.fasta` (override with `--rdna_library`).

## Pipeline / code

- `tc_utils.blastn_rdna_reference_coverage()` — blastn library vs a subject,
  returns per-subject best reference coverage per top type.
- `tc_utils.assign_rdna_to_trcs()` — aggregate to per-TRC calls (max coverage
  per type across the TRC's sequences, threshold, pick the higher type).
- `tc_utils.identify_rdna()` — orchestration: consensus subject → blastn →
  genomic fallback for consensus-less TRCs → write attributes + `_rdna.tsv`.
- Wired into `TideCluster.tarean()` (via `_maybe_identify_rdna`) right before the
  v2 report is built, so `run_all` produces rDNA labels by default. Standalone
  `TideCluster.py rdna -pr <prefix> -f <fasta>` re-runs on an existing run (e.g.
  after extending the library) and re-renders the report.
- Report (`tc_rerender_report.py`): an "rRNA / rDNA tandem repeats" list and a
  summary row on the index page, an `rDNA` row on each TRC page, and
  `rDNA_type`/`rDNA_coverage` columns in `trc_table.tsv`.

Thresholds: `--rdna_min_coverage` (default 0.7), `--rdna_min_identity`
(default 85). On the calibration genome these separated perfectly (rDNA TRCs at
coverage 1.0 / ~90 % identity; all non-rDNA TRCs at 0.0).

## Calibration (test contig OZ408684.1, 27.8 Mb)

- The contig carries one ~5.4 Mb **45S** array (18S/5.8S/25S in canonical order)
  and only single 5S genes (no tandem 5S array).
- 45S monomers are ~9–12 kb, so a clean 45S TRC needs TideHunter `--long`
  (default caps period at 3000 bp).
- `run_all --long` fragmented the 5.4 Mb array into **four TRCs**
  (TRC_1/3/4/5); `identify_rdna` labelled exactly those four as 45S (coverage
  1.0) and left TRC_2/6/7 (other satellites) unlabelled — all via consensus, no
  fallback needed.

## Out of scope here — overlapping annotations (separate problem)

The four 45S TRCs above **overlap each other genomically** (related rDNA
variants over one array), which violates TideCluster's "annotate each region
once" aim. That is a general clustering-overlap problem, tracked separately from
rDNA labelling; this feature only *labels*. The calibration run above is the
concrete fixture for that workstream.

## Possible follow-ups

- Extend `data/rdna_library.fasta` to cover more Viridiplantae (current matches
  are ~90 % identity, implying room to broaden taxon sampling).
- 5S tandem arrays are not present on the test contig; validate 5S labelling on
  a genome that has a 5S array.
