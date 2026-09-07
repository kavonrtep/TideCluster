#!/usr/bin/env Rscript
# tests/test_absorbed_trc.R
#
# Regression guard for issue #8: tc_comparative_analysis.R aborted with
#   Error in if (is.na(x)[1]) { : missing value where TRUE/FALSE needed
# after the (expensive) search, on a group whose annotation fraction was a
# zero-length vector rather than NA.
#
# How a run gets into that state, with nothing corrupt about it:
# `resolve_trc_overlaps` (default-on since 1.16) DROPS a TRC whose every span is
# won by a dominant overlapping neighbour, but `clustering()` freezes the
# consensus set BEFORE that step (TideCluster.py: cons_cls at ~924, resolution
# at ~941, save_consensus_files at ~955). So an absorbed TRC keeps its
# <prefix>_consensus/ files, gets annotated into <prefix>_annotation.tsv, and
# is absent only from <prefix>_clustering.gff3 -- which is what builds
# total_trc_length. Its weighted annotation length then stays numeric(0).
#
# Two independent guards are asserted:
#   1. get_seq_files() excludes consensus sequences of TRCs the clustering GFF3
#      does not contain (the root-cause fix)
#   2. create_annotation_dataframes() tolerates a zero-length fraction (the
#      backstop, for any other route to an empty vector)

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

# ---------------------------------------------------------------------------
# A miniature run directory: TRC_1 and TRC_2 are real; TRC_9 was absorbed by
# overlap resolution, so it is in the consensus dir but NOT in the GFF3.
# ---------------------------------------------------------------------------
d <- file.path(tempdir(), "absorbed_run")
unlink(d, recursive = TRUE)
dir.create(file.path(d, "tc_consensus"), recursive = TRUE, showWarnings = FALSE)

seq_for <- function(k) paste(rep(c("ACGTTGCA", "TTGGCCAA", "GATTACAG")[[k %% 3 + 1]], 20),
                             collapse = "")
writeLines(c(
  "##gff-version 3",
  "chr1\tTideCluster\ttandem_repeat\t1\t5000\t1\t.\t.\tName=TRC_1;repeat_type=TR",
  "chr1\tTideCluster\ttandem_repeat\t8000\t12000\t1\t.\t.\tName=TRC_2;repeat_type=TR"
), file.path(d, "tc_clustering.gff3"))            # note: no TRC_9

pool <- DNAStringSet(c(seq_for(1), seq_for(2), seq_for(9)))
names(pool) <- c("TRC_1_rep0_chr1_0_rnd1", "TRC_2_rep0_chr1_1_rnd1",
                 "TRC_9_rep0_chr1_2_rnd1")        # TRC_9 survives in the consensus dir
writeXStringSet(pool, file.path(d, "tc_consensus", "consensus_sequences_all.fasta"))

lib <- DNAStringSet(c(seq_for(1), seq_for(9)))
names(lib) <- c("TRC_1#TRC_1", "TRC_9#TRC_9")
writeXStringSet(lib, file.path(d, "tc_consensus_dimer_library.fasta"))

# --- 1. the filter ---------------------------------------------------------
res <- suppressMessages(get_seq_files(d, "S1", "tc"))
kept <- sub("^S1:", "", c(names(res$tc), names(res$th)))
kept_trc <- sort(unique(sub("^(TRC_[0-9]+).*", "\\1", kept)))

if ("TRC_9" %in% kept_trc)
  fail(paste("TRC_9 is absent from clustering.gff3 but survived into the pool:",
             paste(kept_trc, collapse = ", ")))
ok("a TRC absent from clustering.gff3 is excluded from the search pool")

if (!all(c("TRC_1", "TRC_2") %in% kept_trc))
  fail(paste("real TRCs were dropped; kept:", paste(kept_trc, collapse = ", ")))
ok("TRCs present in clustering.gff3 are kept")

# a run whose GFF3 lists every TRC must be untouched
writeLines(c(
  "##gff-version 3",
  "chr1\tTideCluster\ttandem_repeat\t1\t5000\t1\t.\t.\tName=TRC_1;repeat_type=TR",
  "chr1\tTideCluster\ttandem_repeat\t8000\t12000\t1\t.\t.\tName=TRC_2;repeat_type=TR",
  "chr1\tTideCluster\ttandem_repeat\t20000\t21000\t1\t.\t.\tName=TRC_9;repeat_type=TR"
), file.path(d, "tc_clustering.gff3"))
res2 <- suppressMessages(get_seq_files(d, "S1", "tc"))
if (length(res2$tc) != 2 || length(res2$th) != 3)
  fail(sprintf("a consistent run was filtered: tc=%d (want 2), th=%d (want 3)",
               length(res2$tc), length(res2$th)))
ok("a run with no absorbed TRC is passed through unchanged")

# --- 2. the backstop guard -------------------------------------------------
grps <- list(df = data.frame(group_id = c(1L, 2L, 3L), stringsAsFactors = FALSE))
annotation_results <- list(
  annot_fractions = list(S1 = list(c(SatA = 0.8),   # normal
                                   NA,              # no annotation
                                   numeric(0))),    # the issue-#8 case
  prevalent_annot = list(S1 = list("SatA", NA, NA)))

out <- tryCatch(create_annotation_dataframes(grps, annotation_results, "S1"),
                error = function(e)
                  fail(paste("zero-length annotation fraction still aborts:",
                             conditionMessage(e))))
if (nrow(out$annot_df) != 3) fail("annotation data frame lost rows")
if (!identical(out$annot_df[["S1_annot"]][3], ""))
  fail(paste("empty fraction should render as \"\", got:",
             dQuote(out$annot_df[["S1_annot"]][3])))
ok("a zero-length annotation fraction renders as empty, without aborting")

cat("\nALL PASS\n")
