#!/usr/bin/env Rscript
# tests/test_trc_similarity.R
#
# Guards trc_similarity_tarean.tsv (compute_trc_similarity_table in
# tc_comparative_analysis.R): the pairwise identity / overlap table of TRCs with
# a TAREAN monomer, within and across samples.
#
# Asserted on two synthetic runs:
#   * only the FIRST library record of a TRC (best TAREAN variant) is used
#   * SSR TRCs and TRCs absent from the clustering GFF3 are excluded
#   * coverage is rotation-independent: a TRC whose consensus starts at a
#     different phase, or on the other strand, still covers 100% both ways.
#     A linear monomer query used to cut such alignments in two.
#   * overlap is asymmetric where it should be: a chimeric 2x-length monomer
#     (A + unrelated) contains all of A, but only half of it is in A
#   * unrelated TRCs produce no row; within-sample pairs are reported
#   * a run with nothing to compare writes a header-only table

suppressWarnings(suppressMessages({
  ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
  if (length(ROOT) == 0 || is.na(ROOT)) ROOT <- normalizePath("..")
  src <- file.path(ROOT, "tc_comparative_analysis.R")
  if (!file.exists(src)) src <- "tc_comparative_analysis.R"
  source(src)
  library(Biostrings)
}))

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")

set.seed(42)
rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
mutate <- function(s, rate) {
  x <- strsplit(s, "")[[1]]
  i <- which(runif(length(x)) < rate)
  x[i] <- sapply(x[i], function(b) sample(setdiff(c("A", "C", "G", "T"), b), 1))
  paste(x, collapse = "")
}
rotate <- function(s, k) paste0(substring(s, k + 1), substring(s, 1, k))
revcomp <- function(s) as.character(reverseComplement(DNAString(s)))
dimer <- function(s) paste0(s, s)

A <- rnd(180)        # satellite family A
U <- rnd(200)        # unrelated
X <- rnd(180)        # unrelated half of the chimera

make_run <- function(dir, lib, gff_rows) {
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  writeLines(c("##gff-version 3", gff_rows), file.path(dir, "tc_clustering.gff3"))
  ds <- DNAStringSet(unname(lib))
  names(ds) <- paste0(names(lib), "#", names(lib))
  writeXStringSet(ds, file.path(dir, "tc_consensus_dimer_library.fasta"))
}
gff_row <- function(trc, type = "TR", start = 1)
  sprintf("chr1\tTideCluster\ttandem_repeat\t%d\t%d\t1\t.\t.\tName=%s;repeat_type=%s",
          start, start + 999, trc, type)

base <- file.path(tempdir(), "trc_similarity")
d1 <- file.path(base, "s1")
d2 <- file.path(base, "s2")

# sample 1: TRC_1 = A (first record) plus a decoy second variant that matches U;
#           TRC_2 = U; TRC_3 = SSR; TRC_9 absorbed (not in GFF3)
lib1 <- c(TRC_1 = dimer(A), TRC_1 = dimer(U), TRC_2 = dimer(U),
          TRC_3 = dimer(strrep("AT", 60)), TRC_9 = dimer(A))
make_run(d1, lib1, c(gff_row("TRC_1"), gff_row("TRC_2", start = 5000),
                     gff_row("TRC_3", "SSR", start = 9000)))
# sample 2: TRC_5 = A at another phase, 5% diverged; TRC_6 = A reverse-complemented
#           and rotated; TRC_7 = chimera A + X (360 bp)
A5 <- mutate(A, 0.05)
lib2 <- c(TRC_5 = dimer(rotate(A5, 70)),
          TRC_6 = dimer(revcomp(rotate(A, 123))),
          TRC_7 = dimer(paste0(A, X)))
make_run(d2, lib2, c(gff_row("TRC_5"), gff_row("TRC_6", start = 3000),
                     gff_row("TRC_7", start = 6000)))

out <- file.path(base, "out")
suppressMessages(compute_trc_similarity_table(c(d1, d2), c("S1", "S2"), c("tc", "tc"),
                                              out, ncpu = 2))
tab <- read.delim(file.path(out, "trc_similarity_tarean.tsv"), stringsAsFactors = FALSE)

want_cols <- c("spec1", "spec2", "trc_spec1", "trc_spec2", "identity",
               "overlap1_in_2", "overlap2_in_1", "monomer_length1", "monomer_length2")
if (!identical(names(tab), want_cols))
  fail(paste("unexpected columns:", paste(names(tab), collapse = ",")))
ok("column layout")

row_of <- function(s1, t1, s2, t2) {
  r <- tab[tab$spec1 == s1 & tab$trc_spec1 == t1 & tab$spec2 == s2 & tab$trc_spec2 == t2, ]
  if (nrow(r) != 1) fail(sprintf("expected exactly one row %s:%s vs %s:%s, got %d",
                                 s1, t1, s2, t2, nrow(r)))
  r
}
present <- unique(c(paste(tab$spec1, tab$trc_spec1), paste(tab$spec2, tab$trc_spec2)))

if (any(c("S1 TRC_3", "S1 TRC_9") %in% present))
  fail("SSR TRC_3 or absorbed TRC_9 made it into the table")
ok("SSR TRCs and TRCs absent from the clustering GFF3 are excluded")

if (any(tab$trc_spec1 == "TRC_1" & tab$trc_spec2 == "TRC_2" & tab$spec2 == "S1"))
  fail("TRC_1 matched TRC_2 through its second (non-best) library variant")
ok("only the first (best) TAREAN variant of a TRC is used")

r <- row_of("S1", "TRC_1", "S2", "TRC_5")
if (r$overlap1_in_2 < 0.99 || r$overlap2_in_1 < 0.99)
  fail(sprintf("rotated copy not fully covered: %.3f / %.3f", r$overlap1_in_2, r$overlap2_in_1))
if (r$identity < 90 || r$identity > 99)
  fail(sprintf("5%%-diverged copy should be ~95%% identical, got %.2f", r$identity))
ok(sprintf("rotation-independent coverage (identity %.2f)", r$identity))

r <- row_of("S1", "TRC_1", "S2", "TRC_6")
if (r$overlap1_in_2 < 0.99 || r$overlap2_in_1 < 0.99 || r$identity < 99.9)
  fail(sprintf("reverse-complement copy: %.2f %.3f %.3f",
               r$identity, r$overlap1_in_2, r$overlap2_in_1))
ok("reverse-complemented, rotated copy is a full match")

r <- row_of("S1", "TRC_1", "S2", "TRC_7")
if (r$overlap1_in_2 < 0.99 || abs(r$overlap2_in_1 - 0.5) > 0.05)
  fail(sprintf("chimera overlaps should be ~1 / ~0.5, got %.3f / %.3f",
               r$overlap1_in_2, r$overlap2_in_1))
if (r$monomer_length1 != 180 || r$monomer_length2 != 360)
  fail("monomer lengths wrong")
ok("asymmetric overlap for a chimeric monomer, monomer lengths reported")

invisible(row_of("S2", "TRC_5", "S2", "TRC_6"))
ok("within-sample pairs are reported")

if (any(c(tab$trc_spec1, tab$trc_spec2)[c(tab$spec1, tab$spec2) == "S1"] == "TRC_2"))
  fail("unrelated TRC_2 produced a row")
ok("unrelated TRCs produce no row")

# nothing to compare: header-only table
make_run(d1, c(TRC_2 = dimer(U)), gff_row("TRC_2"))
out2 <- file.path(base, "out_empty")
suppressMessages(compute_trc_similarity_table(d1, "S1", "tc", out2))
lines <- readLines(file.path(out2, "trc_similarity_tarean.tsv"))
if (length(lines) != 1 || lines[1] != paste(want_cols, collapse = "\t"))
  fail("a run with a single TRC should write a header-only table")
ok("header-only table when there is nothing to compare")

cat("\nALL PASS\n")
