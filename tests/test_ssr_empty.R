#!/usr/bin/env Rscript
# tests/test_ssr_empty.R
#
# Regression guard: cluster_ssrs_sequences() used to crash when NO sample
# contained an SSR TRC.  `all_patterns` is then empty, so `result_df` was built
# with zero rows and the very next line — result_df[[<col>]] <- NA — failed with
#
#   Error in `[[<-.data.frame`(...) : replacement has 1 row, data has 0
#
# taking the whole comparative run down.  It is reachable on any genome pair
# whose TideCluster runs found no SSR (SSRS_summary.csv is then a 0-byte file),
# which is exactly what the bundled CEN6 fixture produces.
#
# Asserts the empty case returns a correctly SHAPED empty frame, so downstream
# code sees "no SSR groups" rather than an error.

suppressWarnings(suppressMessages({
  ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
  if (length(ROOT) == 0 || is.na(ROOT)) ROOT <- normalizePath("..")
  src <- file.path(ROOT, "tc_comparative_analysis.R")
  if (!file.exists(src)) src <- "tc_comparative_analysis.R"
  source(src)
}))

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }

# Two samples, no SSR rows at all — the state an empty SSRS_summary.csv produces.
empty_tbl <- data.frame(trc_id = character(0), ssrs_type = character(0),
                        stringsAsFactors = FALSE)
data_list <- list(sample_a = empty_tbl, sample_b = empty_tbl)

res <- tryCatch(cluster_ssrs_sequences(data_list),
                error = function(e) fail(paste("crashed on all-empty SSR input:",
                                               conditionMessage(e))))

if (!is.data.frame(res)) fail("result is not a data.frame")
if (nrow(res) != 0)      fail(paste("expected 0 rows, got", nrow(res)))

expected <- c("cluster_index", "sample_a_trc_id", "sample_b_trc_id",
              "sample_a_ssrs_type", "sample_b_ssrs_type", "major_pattern")
missing <- setdiff(expected, colnames(res))
if (length(missing) > 0) {
  cat("Columns:", paste(colnames(res), collapse = ", "), "\n")
  fail(paste("empty result is missing column(s):", paste(missing, collapse = ", ")))
}

# One sample with an SSR, the other without, must still work (the mixed case).
mixed <- list(
  sample_a = data.frame(trc_id = "TRC_5", ssrs_type = "AAT (90.6%)",
                        stringsAsFactors = FALSE),
  sample_b = empty_tbl)
res2 <- tryCatch(cluster_ssrs_sequences(mixed),
                 error = function(e) fail(paste("crashed on one-sided SSR input:",
                                                conditionMessage(e))))
if (nrow(res2) != 1) fail(paste("mixed case: expected 1 row, got", nrow(res2)))
if (is.na(res2$sample_a_trc_id[1])) fail("mixed case: sample_a TRC not recorded")
if (!is.na(res2$sample_b_trc_id[1])) fail("mixed case: sample_b should be NA")

cat("PASS: cluster_ssrs_sequences handles samples with no SSRs at all\n")
