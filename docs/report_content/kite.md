<!-- Appears on: <prefix>_report/kite.html (top). Edit freely. -->
## KITE method

KITE (K-mer Interval Tandem-repeat Estimation) analyses each TRA
independently. It collects the distances between consecutive
occurrences of every *k*-mer, builds a combined *k*-mer-interval
profile, and identifies its strongest peaks — each candidate monomer
period is reported as one of the per-array monomer-size estimates
(top peaks, ranked by score).

Since TideCluster 1.10.0 the step is performed by
[kitehor](https://github.com/kavonrtep/kitehor) (Rust). For each
array kitehor also returns a Founder period, the higher-order
Strongest period and its multiplicity ×k, candidate subrepeats with
an occupancy tier, and a dominant SSR motif. These signals drive the
per-TRC dashboards.

KITE is complementary to TAREAN: TAREAN returns one consensus per
TRC (family); KITE returns a per-TRA (family-member) monomer estimate
and structural signals. The per-array view exposes variation in
monomer length and higher-order structure across members of one
family.
