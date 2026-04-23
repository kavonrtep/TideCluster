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

- **Family (TRC, Tandem Repeat Cluster)** — a group of tandem repeat
  arrays whose monomer sequences are nearly identical. Arrays are
  assigned to the same family when their monomers align at
  ≥ 75 % sequence identity over ≥ 80 % of the shorter sequence;
  similarity is propagated transitively, so if array A is similar to
  B and B to C at these thresholds, A, B and C are one family. TRCs
  are numbered `TRC_1`, `TRC_2`, … in decreasing order of total array
  length.
- **TRA (Tandem Repeat Array)** — a single genomic instance
  belonging to a family. Each TRA has a seqid, start, end, and
  length.
- **Superfamily** — a group of families whose consensus sequences
  share any BLAST-detectable similarity. Because this criterion is
  looser than family assignment, superfamilies can join families
  that share only a local region of their consensus, or families
  whose monomers differ in length — for example when one family's
  monomer is contained inside the longer monomer of another family.
  Higher-order repeat relationships are included but are not the
  only source of superfamily membership.
- **HOR (Higher-Order Repeat)** — a tandem repeat whose repeating
  unit is itself a block of several basic monomers. Classified per
  TRA on a continuous confidence score; see
  `docs/hor_classification.md` for the formula.
- **Monomer** — the basic repeating unit of a tandem repeat. Its
  length in bp is the monomer size.

> In satellite-rich genomes tandem repeats often form a continuum of
> related sequences rather than discrete groups; the family and
> superfamily boundaries reported here are operational summaries at
> fixed thresholds, not clear-cut biological clades.
