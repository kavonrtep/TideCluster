<!-- Appears on: <prefix>_report/tarean.html (top). Edit freely. -->
## TAREAN method

TAREAN (Tandem Repeat Analyzer) reconstructs a consensus sequence for
each TRC using a *k*-mer-based de Bruijn graph approach. For every
TRC above the minimum total length threshold, TAREAN extracts all
*k*-mers across several sizes (*k* = 11, 19, 23, 27), builds a de
Bruijn graph, searches for cycles representing candidate consensus
sequences, scores each candidate by *k*-mer coverage, and returns
the highest-scoring variant.

The table shows one row per TRC (family). Key columns:

- **Type** — `TR` for tandem repeat clusters, `SSR` for clusters whose
  monomer is a simple sequence repeat.
- **Monomer (TAREAN)** — primary TAREAN consensus length in bp.
- **Score (TAREAN)** — *k*-mer coverage score of the chosen
  consensus. Higher values indicate a tighter consensus.
- **Array sizes (min / med / max)** — length range of the TRAs in
  this family, in bp.
- **Graph / Logo** — thumbnails; click for the full image. Both are
  produced by TAREAN for the primary consensus variant.

For SSR-type TRCs and TRCs below the length threshold the TAREAN
columns (Monomer, Score, Graph, Logo, Consensus) are empty;
cluster statistics still show.

Per-array HOR / subrepeat / SSR signals (formerly summarised here)
now live on each TRC's dashboard, reachable by clicking the TRC id.
