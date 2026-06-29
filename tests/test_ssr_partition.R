#!/usr/bin/env Rscript
# tests/test_ssr_partition.R
#
# Regression guard for the Part 1 "double-assigned TRC" bug in
# tc_comparative_analysis.R (see docs/tidecluster_comparative_dev_report.md):
# a TRC that spanned several communities kept multiple rows in trc_groups,
# and apply_ssr_grouping() only reassigned the FIRST row, so the TRC ended up
# in two satellite families (the SSR family + a stale similarity family).
#
# This test feeds apply_ssr_grouping() a trc_groups that still has DUPLICATE
# rows for one TRC (simulating the pre-fix-B state) and asserts that the
# result is a clean partition: every TRC maps to exactly one group_id. It
# guards Fix A (reassign every matching row, not just the first). Fix B
# (collapse to one row per TRC in cluster_trc_sequences) makes the duplicate
# state impossible upstream; this test keeps Fix A honest if it ever recurs.

suppressWarnings(suppressMessages({
  ok <- TRUE

  # Source the pipeline; the sys.nframe() guard must stop main() from running.
  ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
  if (length(ROOT) == 0 || is.na(ROOT)) ROOT <- normalizePath("..")
  src <- file.path(ROOT, "tc_comparative_analysis.R")
  if (!file.exists(src)) src <- "tc_comparative_analysis.R"
  source(src)
}))

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }

if (!exists("apply_ssr_grouping", mode = "function")) {
  fail("apply_ssr_grouping() not found after source() (source guard broken?)")
}

# --- synthetic input -------------------------------------------------------
# Two samples A and B. A:TRC_1 deliberately appears in TWO rows (groups 1 & 2)
# -- the pre-fix-B duplicate state. B:TRC_1 is its cross-genome partner.
# A:TRC_2 is an unrelated similarity-only TRC that must stay put.
trc_groups <- data.frame(
  trc_id   = c("A:TRC_1", "A:TRC_1", "B:TRC_1", "A:TRC_2"),
  group_id = c(1L,         2L,        1L,        3L),
  stringsAsFactors = FALSE
)

# One SSR pattern shared by A:TRC_1 and B:TRC_1 -> one SSR family.
ssrs_groups <- data.frame(
  A_trc_id = "TRC_1",
  B_trc_id = "TRC_1",
  stringsAsFactors = FALSE
)

out_dir <- tempfile("ssr_partition_")
dir.create(out_dir)  # no trc_graph.rds inside -> graph sync is skipped

res <- suppressMessages(
  apply_ssr_grouping(trc_groups, ssrs_groups, prefix = c("A", "B"),
                     output_directory = out_dir)
)

# --- assertions ------------------------------------------------------------
# 1) Partition: every TRC maps to exactly one group_id.
groups_per_trc <- tapply(res$group_id, res$trc_id, function(x) length(unique(x)))
double_assigned <- names(groups_per_trc)[groups_per_trc > 1]
if (length(double_assigned) > 0) {
  cat("Result table:\n"); print(res)
  fail(sprintf("double-assigned TRC(s): %s", paste(double_assigned, collapse = ", ")))
}

# 2) A:TRC_1 and B:TRC_1 ended up in the SAME (SSR) family.
g_a <- unique(res$group_id[res$trc_id == "A:TRC_1"])
g_b <- unique(res$group_id[res$trc_id == "B:TRC_1"])
if (length(g_a) != 1 || length(g_b) != 1 || g_a != g_b) {
  cat("Result table:\n"); print(res)
  fail("A:TRC_1 and B:TRC_1 were not grouped into the same SSR family")
}

# 3) The unrelated similarity TRC is still present and on its own.
if (!("A:TRC_2" %in% res$trc_id)) fail("A:TRC_2 vanished from trc_groups")

cat("PASS: apply_ssr_grouping produces a clean partition (no double-assignment)\n")
