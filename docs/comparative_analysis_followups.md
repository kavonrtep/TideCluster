# tc_comparative_analysis.R — deferred follow-ups

Captured 2026-06-01 alongside the silent-drop fix on branch
`fix/comparative-silent-drop`. The primary bug — TRCs with no surviving
MMseqs2 hit being silently dropped from every comparative output — is
fixed in that branch by (a) adding all input sequences as graph vertices
in `cluster_trc_sequences()` and (b) regrouping pure-SSR TRCs by
`ssrs_groups` pattern in `apply_ssr_grouping()`. A post-condition
warning (`check_trc_coverage()`) flags any TRC from `tc_clustering.gff3`
not represented in either `trc_satellite_families.tsv` or
`ssrs_groups.tsv`.

The two secondary issues called out by the original bug report were
**deferred** out of that branch to keep its blast radius small.

## 1. Threshold mismatch with the manuscript

`process_trc_analysis()` (`tc_comparative_analysis.R`) currently calls
`cluster_trc_sequences()` with:

```r
min_identity = 80,
min_coverage = 0.8       # applied as max(qcov, tcov) >= 0.8
```

The manuscript Implementation text states **"70 % identity and 20 %
bidirectional coverage"**. These disagree; one of them is wrong.
(Sensitivity 7.5 *does* match.) Action: pick one source of truth and
align the other. If the code wins, the manuscript text needs an
erratum; if the manuscript wins, the defaults change and clustering
output will shift for everyone — that change should land in its own
release with a CHANGELOG note.

## 2. Short-monomer satellites fall between classifiers

Several TRCs with 2–7 bp monomers (e.g. TRC_18, TRC_29, TRC_50, TRC_59
in the manuscript's *S. lycopersicum* example) are **not** in
`ssrs_groups.tsv` and were also dropped from `trc_satellite_families.tsv`
under the bug. After this branch lands they will appear as singletons,
but the underlying handoff issue remains: they are neither classified as
SSRs (under the 10 % `min_percentage` floor in
`cluster_ssrs_sequences()`) nor as proper satellites by the
similarity-based pipeline. Action: review the SSR/satellite hand-off
thresholds — `min_percentage`, and how the SSRS_summary.csv is built
upstream — so these TRCs land on one side or the other deliberately.

---

Both are tracked here rather than in a memory because they need
discussion before any behavioral change.
