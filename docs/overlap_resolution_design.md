# Cross-TRC overlap resolution — design notes

TideCluster's aim is to **annotate each genomic region once**, but the
clustering GFF3 could place several TRCs over the same bases: variant arrays of
a satellite (most visibly rDNA) get clustered into separate TRCs that interleave
and overlap at their boundaries. This is a general clustering-output issue,
tracked separately from rDNA labelling (which only *labels*; see
[`rdna_design.md`](rdna_design.md)).

## Decisions (agreed with maintainer)

- **Scope: general.** Resolve any cross-TRC overlap, not just rDNA.
- **Strategy: dominant TRC wins the contested span.** Split at overlap
  breakpoints; assign each elementary segment to the covering TRC with the
  largest **total array length** over the whole genome (tie-break: `Name`).
  Per-span TideHunter score is not present in the clustering GFF3, so total
  array length is the dominance weight (the maintainer's stated criterion).
- **Placement: rewrite the clustering GFF3, default-on.** Runs at the end of
  `clustering()` after the GFF3 is finalised; disabled with `--keep_overlaps`.

## Implementation

`tc_utils.resolve_trc_overlaps(gff_in, gff_out)` — a per-seqid breakpoint sweep.
No base of the union is lost (only reassigned); adjacent segments kept by the
same `(Name, strand)` are merged back into single features; all per-TRC
attributes (`repeat_type`, `ssr`, `rDNA_type`, …) are carried through. Output is
sorted and deterministic, and the operation is idempotent.

Wired into `TideCluster.clustering()` (default-on) for the `clustering` and
`run_all` subcommands; `--keep_overlaps` opts out.

**Caveat — absorbed TRCs.** A TRC whose arrays are *entirely* nested inside
larger dominant TRCs loses all its bases and is dropped from the GFF3; the run
logs `resolve_trc_overlaps: N TRC(s) fully absorbed …`. In practice this is rare
(fully-nested distinct TRCs would usually have clustered together); the observed
overlaps are boundary transitions, where no TRC vanishes.

## Validation (calibration genome OZ408684.1)

The 5.4 Mb 45S array fragments into four variant TRCs (TRC_1/3/4/5) that overlap
each other:

| | covered bp (union) | multiply-covered bp |
|---|---:|---:|
| before | 6,203,074 | 239,011 (3.9 %) |
| after  | 6,203,074 | 0 |

Union preserved exactly; all 239 kb of contested boundary span reassigned to the
dominant variant (TRC_1, the largest, keeps its full extent; TRC_3 −148 kb,
TRC_5 −68 kb, TRC_4 −22 kb; a single 84 bp TRC_2/TRC_7 touch resolved too). No
TRC was absorbed.

## Possible follow-up

- Upgrade the dominance weight from global total-array-length to a **local**
  signal (per-array TideHunter score, or the longer overlapping array at that
  locus) if boundary placement needs to be more locally faithful. Would require
  threading per-array scores from the TideHunter GFF3 into the clustering GFF3.
