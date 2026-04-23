<!-- Appears on: <prefix>_report/kite.html (top). Edit freely. -->
## KITE method

KITE (K-mer Interval Tandem Repeat Estimation) analyses each TRA
independently. For a given array KITE extracts all *k*-mers (default
*k* = 6), computes the distribution of distances between consecutive
occurrences of every *k*-mer, aggregates those into a combined
*k*-mer-interval profile, and identifies significant peaks. Each peak
position corresponds to a candidate monomer period. The top three
peaks by score are the primary, secondary, and tertiary monomer-size
estimates *m*₁, *m*₂, *m*₃, with weighted scores *s*₁, *s*₂, *s*₃.

KITE is complementary to TAREAN: TAREAN returns a single consensus
per TRC (family), while KITE returns a per-TRA (family member)
monomer estimate. The per-array view exposes variation in monomer
length across members of the family and feeds the HOR classifier
described in the next section. For the full algorithm, thresholds,
and worked examples see `docs/hor_classification.md`.
