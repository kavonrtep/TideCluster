<!-- Appears on: <prefix>_index.html (top). Edit freely. -->
## About this report

TideCluster identifies tandem repeats in genome assemblies using
TideHunter for detection and mmseqs2 / BLAST for clustering. Detected
repeats are grouped into **Tandem Repeat Clusters (TRCs)** — **families**
of related tandem repeats that share sequence similarity. Each TRC
contains one or more **Tandem Repeat Arrays (TRAs)** — the individual
genomic instances of the family. One family can therefore span many
chromosomes and contigs.

Report sections (top navigation):

- **TAREAN** — per-TRC (per-family) consensus sequences, array
  statistics, and HOR summary.
- **KITE** — per-TRA monomer-size estimates and HOR classification.
- **Superfamilies** — higher-level groupings of TRCs whose consensus
  sequences share BLAST similarity.
- **Per-TRC dashboards** — one page per family, reached from any TRC
  link in the three main tables.

The original v1 report is preserved under
`<prefix>_report_legacy/`; reach it via the **Legacy report**
link at the top-right of the navigation bar.

### Glossary

- **TRC (Tandem Repeat Cluster)** — a family of related tandem
  repeats. Numbered `TRC_1`, `TRC_2`, … in decreasing order of total
  array length.
- **TRA (Tandem Repeat Array)** — a single genomic instance
  belonging to a family. Each TRA has a seqid, start, end, and
  length.
- **HOR (Higher-Order Repeat)** — a tandem repeat whose repeating
  unit is itself a block of several basic monomers. Classified per
  TRA on a continuous confidence score; see
  `docs/hor_classification.md` for the formula.
- **Superfamily** — a group of TRCs whose TAREAN consensus
  sequences share significant BLAST similarity. A second level of
  the repeat hierarchy above the family.
- **Monomer** — the basic repeating unit of a tandem repeat. Its
  length in bp is the monomer size.
