<!-- Appears on: <prefix>_report/kite.html (top). Edit freely. -->
## KITE method

KITE (K-mer Interval Tandem Repeat Estimation) analyses each TRA
independently. For a given array KITE extracts all *k*-mers (default
*k* = 6), computes the distribution of distances between consecutive
occurrences of every *k*-mer, aggregates those into a combined
*k*-mer-interval profile, and identifies significant peaks. Each peak
position corresponds to a candidate monomer period; the top peaks by
score are reported per array as "Monomer size estimates (Score)" (up
to five candidates, ranked best-first).

Since TideCluster 1.10.0 the KITE step is performed by
[kitehor](https://github.com/kavonrtep/kitehor), which augments the
periodicity scan with a structural classifier. Each array gets a
single `combined_class` — `hor`, `hor_with_ssr`, `tr`, `tr_with_ssr`,
`tr_with_subrepeat`, `pure_ssr`, or `unresolved` — together with the
HOR founder/tile periods and multiplicity, the tandem-validate
host/subrepeat periods, and an SSR scan. These drive the class badges
and per-array "Structure" cells throughout the report.

KITE is complementary to TAREAN: TAREAN returns a single consensus
per TRC (family), while KITE returns a per-TRA (family member)
monomer estimate and structural class. The per-array view exposes
variation in monomer length and higher-order structure across members
of the family.
