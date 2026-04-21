## 1.9.0 (2026-04-21)
- KITE now classifies each tandem repeat array into one of three HOR
  (higher-order repeat) categories — `No HOR detected`, `HOR-visible`,
  `HOR-dominant` — using integer-multiple relationships between the
  top three monomer-size estimates. Per-TRC counts in each category
  are reported in the KITE HTML summary and TSV. Thresholds exposed
  as constants at the top of `tarean/kite.R`.
- New `tc_rerender_report.py` CLI builds an opt-in report v2 from an
  existing TideCluster output directory. Emits `<prefix>_report_v2/`
  containing modern DataTables-powered HTML pages (summary, all TRCs,
  TAREAN, KITE, superfamilies) plus one per-TRC dashboard each,
  offline-first, no new conda dependencies. The original reports are
  left untouched.
- Pipeline now writes `<prefix>_pipeline_stats.json` alongside the
  index HTML so downstream tools have the summary numbers without
  having to scrape HTML.
- Test suite gains a `rerender` tier.

## 1.8.3 (2026-04-18)
- Conda releases are now built from this repository via GitHub Actions on
  tag push; the external `kavonrtep/recipes/tidecluster` recipe is retired.
- Added `tests.sh` dispatcher with smoke/short/long tiers gating CI and
  release builds.

## 1.7.1 (2025-10-09)
- Gzip FASTA input support added
- Interactive visualization and comparative analysis HTML reports implemented
- Comparative analysis of multiple assemblies was added

## 1.6.2 (2025-07-10)
- RepeatMasker wrapper function added to handle long sequence names
- Documentation updated with additional output files
## 1.6.1 (2025-05-03)
- Bugfix in tarean run
- Extract consensus functionality added
## 1.6 (2024-11-05)
- TRC superfamily clustering added (consensus-vs-consensus BLAST analysis)
- HTML report improvements and bugfixes
- Documentation updated
## 1.5 (2024-03-07)
- K-mer interval tandem repeat estimation (KITE) added, this enables identification of higher order tandem repeats
- Analysis include estimation of monomer size for individual tandem repeat arrays
- Output of KITE analysis added to the HTML report
- Bug fix in sorting of output CSV tables

## 1.4 (2023-11-23)
- Better memory efficiency
- Limiting MMseq2 memory usage
## 1.3 (2023-10-23)
- Improved input FASTA parsing to decrease memory usage
## 1.2 (2023-10-06)
- Bug fix in SSRs detection
- Better SSRs reporting TAREAN output
## 1.1 (2023-09-07)
- Bugfix - better handling of SSR
## 1.0 (2023-08-17)
- Bugfix in HTML report
- Tidehunter run in batches to decrease memory usage
- Documentation updated
## 0.9 (2023-08-11)
- Documentation updated
- Info related to bug #1 added to documentation
- Exporting TR library added
## 0.0.8 (2023-05-02)
- Documentation update, HTML bugfix
- Workflow scheme added
- Workflow scheme formatted
- Update_gff script added
## 0.0.7 (2023-04-26)
- Output files added
## 0.0.6 (2023-04-25)
- Bug-fix - handle lower case character in DNA sequences
- --min_total_size option added, some bugfixes
## 0.0.5 (2023-04-17)
- Bug fix in logo plotting
- Clustering of SSRs added
- Format of output GFF corrected
## 0.0.4 (2023-04-05)
- Bugfixes, better defaults, documentation updated
## 0.0.3 (2023-03-23)
- Run all option added, bugfixes
## 0.0.2 (2023-03-23)
- Bug fixes
