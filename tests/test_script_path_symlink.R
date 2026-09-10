#!/usr/bin/env Rscript
# tests/test_script_path_symlink.R
#
# Regression guard for issue #9: tc_summarize_comparative_analysis.R could not
# find its html/ templates when invoked through the conda `bin` symlink.
#
#   Error in create_javascript_files(output_dir) :
#     Template directory not found: /opt/conda/bin/html
#
# R's `--file=` holds the path AS INVOKED and does not resolve symlinks. Both
# deployments put each CLI in bin/ as a symlink to the real tree --
#   conda: $PREFIX/bin/x -> $PREFIX/share/tidecluster/x
#   SIF:   /opt/conda/bin/x -> /opt/tidecluster/x
# -- so dirname() of the raw value is bin/, where html/ does not exist. Only a
# source checkout, where script and assets are siblings, happened to work.
#
# Builds that exact layout and runs each affected script's OWN resolution code
# through the symlink, by appending a one-line reporter to a copy of it. No
# parsing of the script's contents beyond that, so it cannot go stale silently.

fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
if (length(ROOT) == 0 || is.na(ROOT)) ROOT <- normalizePath("..")

# script -> the asset directory it looks up beside itself (NA = none).
# tc_summarize_comparative_analysis.R is the shipped CLI (conda build.sh and
# TideCluster.def both symlink it into bin/) and is REQUIRED here. The other two
# carry the same idiom but are untracked dev helpers, so they are checked only
# when present rather than failing a clean checkout.
SCRIPTS <- list(
  "tc_summarize_comparative_analysis.R" = "html",
  "tidecluster_viz.R"                   = "html_template",
  "plot_karyotype.R"                    = NA
)
REQUIRED <- "tc_summarize_comparative_analysis.R"

base  <- file.path(tempdir(), "tc_symlink_layout")
unlink(base, recursive = TRUE)
share <- file.path(base, "share", "tidecluster")
bin   <- file.path(base, "bin")
dir.create(share, recursive = TRUE); dir.create(bin, recursive = TRUE)

# A truncated copy of each script: everything up to and including the
# script_path block, plus a line that prints what it resolved to. Truncating
# avoids running the script's real work (which needs args, input dirs, ...).
truncate_after_resolution <- function(src, dst) {
  lines <- readLines(src, warn = FALSE)
  # Anchor on the script_path assignment itself, not on commandArgs: some
  # scripts parse their options first, and cutting there would run that code.
  # Matches the fixed form (`sub(...)` then normalizePath) and the pre-fix one
  # (`dirname(sub(...))`) alike, so this test genuinely fails against old code
  # instead of erroring out before it tests anything.
  start <- grep('^script_path <- (dirname\\()?sub\\("--file=', lines)[1]
  if (is.na(start)) fail(paste("no script_path block in", basename(src)))
  close <- grep("^\\}\\s*$", lines)
  close <- close[close > start][1]
  if (is.na(close)) fail(paste("could not find the end of the block in", basename(src)))
  writeLines(c('args <- commandArgs(trailingOnly = FALSE)',
               lines[seq(start, close)],
               'cat(script_path, "\n")'), dst)
}

present <- names(SCRIPTS)[file.exists(file.path(ROOT, names(SCRIPTS)))]
missing_required <- setdiff(REQUIRED, present)
if (length(missing_required) > 0)
  fail(paste("missing from the repo:", paste(missing_required, collapse = ", ")))
for (nm in setdiff(names(SCRIPTS), present))
  cat("SKIP:", nm, "(not in this checkout)\n")

for (nm in present) {
  src <- file.path(ROOT, nm)
  truncate_after_resolution(src, file.path(share, nm))
  if (!file.symlink(file.path(share, nm), file.path(bin, nm)))
    fail(paste("could not create the bin symlink for", nm))
}
for (d in unique(stats::na.omit(unlist(SCRIPTS[present])))) {
  file.copy(file.path(ROOT, d), share, recursive = TRUE)
  if (!dir.exists(file.path(share, d)))
    fail(paste("asset dir did not copy into the fake tree:", d))
}

want <- normalizePath(share)
for (nm in present) {
  for (how in c("bin symlink", "real path")) {
    invoke <- file.path(if (how == "bin symlink") bin else share, nm)
    out <- suppressWarnings(system2("Rscript", invoke, stdout = TRUE, stderr = TRUE))
    got <- trimws(paste(out, collapse = " "))
    if (!nzchar(got) || grepl("Error", got, fixed = TRUE))
      fail(sprintf("%s via %s: script errored: %s", nm, how, got))
    got <- normalizePath(got, mustWork = FALSE)
    if (!identical(got, want))
      fail(sprintf("%s via %s: script_path = %s, want %s", nm, how, got, want))
  }
  ok(sprintf("%-38s resolves to the real tree (symlink and real path)", nm))

  asset <- SCRIPTS[[nm]]
  if (!is.na(asset) && !dir.exists(file.path(want, asset)))
    fail(sprintf("%s: %s/ not found under the resolved path", nm, asset))
  if (!is.na(asset))
    ok(sprintf("%-38s finds its %s/ assets there", nm, asset))
}

cat("\nALL PASS\n")
