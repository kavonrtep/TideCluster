#!/usr/bin/env Rscript
# Phase-aligned PFM (Fourier-style) consensus prototype for KITE.
#
# CLI driver. Two modes:
#   --mode single : one TRA FASTA + m, write consensus / PFM / diagnostics.
#   --mode batch  : walk an existing KITE TSV + TRC fasta dir, produce
#                   combined per-TRA consensus FASTA + diagnostics TSV.
#
# Prototype only. Not wired into TideCluster.py or tarean/kite.R; see
# docs/kite_consensus_implementation_plan.md for the graduation path.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(parallel)
})

# Locate this script and source its lib/
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .args[grep("--file=", .args)]
.script_dir <- if (length(.file_arg) > 0) {
  normalizePath(dirname(sub("--file=", "", .file_arg)))
} else {
  getwd()
}
source(file.path(.script_dir, "lib", "consensus_core.R"))

# -------- CLI ----------------------------------------------------------

option_list <- list(
  make_option("--mode",          type = "character", default = "batch",
              help = "single | batch"),
  make_option("--fasta",         type = "character", default = NULL,
              help = "single: TRA FASTA file with one DNA record"),
  make_option("--m",             type = "integer",   default = NULL,
              help = "single: monomer length (bp)"),
  make_option("--out",           type = "character", default = NULL,
              help = "single: output prefix"),
  make_option("--kite-tsv",      type = "character", default = NULL,
              help = "batch: path to monomer_size_top3_estimats.csv"),
  make_option("--tarean-dir",    type = "character", default = NULL,
              help = "batch: path to <prefix>_tarean/fasta/ with TRC_*.fasta"),
  make_option("--out-dir",       type = "character", default = NULL,
              help = "batch: output directory (created if absent)"),
  make_option("--cpu",           type = "integer",   default = 4L,
              help = "batch: parallel worker count (default 4)"),
  make_option("--window-size",   type = "integer",   default = 50000L,
              help = "sampling window size in bp (default 50000)"),
  make_option("--n-windows",     type = "integer",   default = 16L,
              help = "number of sampling windows per TRA (default 16)"),
  make_option("--pseudo",        type = "double",    default = 0.5,
              help = "pseudocount for PFM (default 0.5)"),
  make_option("--write-pfm",     type = "logical",   default = FALSE,
              help = "batch: also emit one PFM TSV per TRA (default FALSE)",
              action = "store_true")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- helpers shared by both modes ----

format_fasta <- function(header, seq, width = 60L) {
  wrapped <- gsub(sprintf("(.{%d})", width), "\\1\n", seq, perl = TRUE)
  paste0(">", header, "\n", wrapped, "\n")
}

write_pfm_tsv <- function(pfm, path) {
  df <- data.frame(position = seq_len(nrow(pfm)),
                   A = pfm[, "A"], C = pfm[, "C"],
                   G = pfm[, "G"], T = pfm[, "T"])
  write.table(df, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE)
}

# Format a diagnostics row as a named list (TSV-friendly).
diag_row <- function(res, extras) {
  harm <- res$harmonic
  row <- list(
    length_bp     = res$length_bp,
    bases_counted = res$bases_counted,
    coverage      = sprintf("%.4f", res$coverage),
    mean_entropy  = sprintf("%.4f", res$mean_entropy),
    consensus_len = nchar(res$consensus)
  )
  for (k in names(harm)) row[[paste0("harm_", k)]] <- sprintf("%.4f", harm[k])
  modifyList(row, extras)
}

# -------- single-array mode -------------------------------------------

run_single <- function(opt) {
  stopifnot(!is.null(opt$fasta), !is.null(opt$m), !is.null(opt$out))
  if (!file.exists(opt$fasta)) stop("fasta not found: ", opt$fasta)
  ds <- readDNAStringSet(opt$fasta)
  if (length(ds) == 0) stop("empty FASTA")
  S <- ds[[1]]
  rec_name <- names(ds)[1]
  res <- compute_tra_consensus(S, opt$m,
                               window_size = opt$`window-size`,
                               n_windows   = opt$`n-windows`,
                               pseudo      = opt$pseudo)
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  writeLines(format_fasta(paste0(rec_name, "  m=", opt$m), res$consensus),
             con = paste0(opt$out, "_consensus.fasta"))
  write_pfm_tsv(res$pfm, paste0(opt$out, "_pfm.tsv"))
  diag <- diag_row(res, list(record = rec_name, m = opt$m))
  write.table(as.data.frame(diag, stringsAsFactors = FALSE),
              file = paste0(opt$out, "_diagnostics.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("wrote:", paste0(opt$out, "_consensus.fasta"), "\n",
      "      ", paste0(opt$out, "_pfm.tsv"), "\n",
      "      ", paste0(opt$out, "_diagnostics.tsv"), "\n")
}

# -------- batch mode --------------------------------------------------

read_kite_tsv <- function(path) {
  # The KITE TSV is tab-separated with a header row. We need TRC_ID,
  # seqid, start, end, monomer_size (as m1), plus any pre-computed HOR
  # fields just to pass through to the output.
  read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

# For a given TRC, map each fasta record to its (seqid, start, end).
index_trc_fasta <- function(fa_path) {
  ds <- readDNAStringSet(fa_path)
  nms <- names(ds)
  # Seqs in the TRC fasta are dimers; KITE takes the first half.
  halves <- lapply(ds, function(x) subseq(x, start = 1L, end = nchar(x) %/% 2L))
  parsed <- lapply(nms, parse_record_name)
  list(halves = halves, parsed = parsed, names = nms)
}

process_one <- function(row, trc_index_cache, window_size, n_windows, pseudo, write_pfm, pfm_dir) {
  trc <- row$TRC_ID
  if (!trc %in% names(trc_index_cache)) {
    return(list(ok = FALSE, reason = paste0("no fasta for ", trc), row = row))
  }
  idx <- trc_index_cache[[trc]]
  # Find the record matching seqid/start/end.
  match_i <- NA_integer_
  for (i in seq_along(idx$parsed)) {
    p <- idx$parsed[[i]]
    if (!is.na(p$seqid) && p$seqid == row$seqid &&
        p$start == row$start && p$end == row$end) {
      match_i <- i; break
    }
  }
  if (is.na(match_i)) {
    return(list(ok = FALSE,
                reason = sprintf("no record for %s:%d-%d in %s",
                                 row$seqid, row$start, row$end, trc),
                row = row))
  }
  S <- idx$halves[[match_i]]
  m <- as.integer(row$monomer_size)
  if (is.na(m) || m < 2L || m > length(S)) {
    return(list(ok = FALSE,
                reason = sprintf("bad m=%s (L=%d)", as.character(m), length(S)),
                row = row))
  }
  res <- compute_tra_consensus(S, m,
                               window_size = window_size,
                               n_windows   = n_windows,
                               pseudo      = pseudo)
  # per-TRA identifier used in FASTA header + PFM filename
  id <- sprintf("%s__%s__%d-%d", trc, row$seqid, row$start, row$end)
  if (write_pfm) {
    write_pfm_tsv(res$pfm, file.path(pfm_dir, paste0(id, ".tsv")))
  }
  diag <- diag_row(res, list(TRC_ID = trc, seqid = row$seqid,
                             start = row$start, end = row$end,
                             m = m, id = id))
  list(ok = TRUE, consensus = res$consensus, id = id, diag = diag)
}

run_batch <- function(opt) {
  for (f in c("kite-tsv", "tarean-dir", "out-dir")) {
    if (is.null(opt[[f]])) stop("batch requires --", f)
  }
  dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)
  pfm_dir <- file.path(opt$`out-dir`, "per_tra_pfm")
  if (isTRUE(opt$`write-pfm`)) dir.create(pfm_dir, showWarnings = FALSE)

  t0 <- Sys.time()
  tsv <- read_kite_tsv(opt$`kite-tsv`)
  message(sprintf("TSV rows: %d ; unique TRCs: %d",
                  nrow(tsv), length(unique(tsv$TRC_ID))))

  # Build per-TRC fasta index lazily.
  unique_trcs <- unique(tsv$TRC_ID)
  trc_index_cache <- setNames(vector("list", length(unique_trcs)), unique_trcs)
  for (trc in unique_trcs) {
    fa <- file.path(opt$`tarean-dir`, paste0(trc, ".fasta"))
    if (!file.exists(fa)) next
    trc_index_cache[[trc]] <- index_trc_fasta(fa)
  }
  present <- sum(!vapply(trc_index_cache, is.null, logical(1)))
  message(sprintf("TRC fastas found: %d / %d", present, length(unique_trcs)))

  rows <- split(tsv, seq_len(nrow(tsv)))
  rows <- lapply(rows, as.list)

  results <- mclapply(
    rows,
    function(r) process_one(r,
                            trc_index_cache = trc_index_cache,
                            window_size     = opt$`window-size`,
                            n_windows       = opt$`n-windows`,
                            pseudo          = opt$pseudo,
                            write_pfm       = isTRUE(opt$`write-pfm`),
                            pfm_dir         = pfm_dir),
    mc.cores = opt$cpu
  )

  ok_mask <- vapply(results, function(x) isTRUE(x$ok), logical(1))
  ok      <- results[ok_mask]
  failed  <- results[!ok_mask]

  fa_path <- file.path(opt$`out-dir`, "per_tra_consensus.fasta")
  con <- file(fa_path, "w")
  for (r in ok) cat(format_fasta(r$id, r$consensus), file = con)
  close(con)

  # Collect diagnostic rows into a single data frame
  if (length(ok) > 0) {
    common_cols <- c("TRC_ID", "seqid", "start", "end", "id", "m",
                     "length_bp", "bases_counted", "coverage",
                     "mean_entropy", "consensus_len")
    harm_cols <- setdiff(names(ok[[1]]$diag),
                         c(common_cols, "record"))
    all_cols <- c(common_cols, sort(harm_cols))
    diag_df <- do.call(rbind, lapply(ok, function(r) {
      as.data.frame(r$diag[all_cols], stringsAsFactors = FALSE)
    }))
    write.table(diag_df,
                file = file.path(opt$`out-dir`, "per_tra_diagnostics.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  fail_log <- file.path(opt$`out-dir`, "failed.tsv")
  if (length(failed) > 0) {
    df <- do.call(rbind, lapply(failed, function(x) {
      data.frame(TRC_ID = x$row$TRC_ID,
                 seqid  = x$row$seqid,
                 start  = x$row$start,
                 end    = x$row$end,
                 m      = x$row$monomer_size,
                 reason = x$reason,
                 stringsAsFactors = FALSE)
    }))
    write.table(df, file = fail_log, sep = "\t", quote = FALSE,
                row.names = FALSE)
  }

  dt <- difftime(Sys.time(), t0, units = "secs")
  summary_path <- file.path(opt$`out-dir`, "summary.log")
  writeLines(c(
    sprintf("total rows:       %d", nrow(tsv)),
    sprintf("processed ok:     %d", length(ok)),
    sprintf("failed:           %d", length(failed)),
    sprintf("wall time (s):    %.1f", as.numeric(dt)),
    sprintf("output FASTA:     %s", fa_path)
  ), summary_path)
  message("summary: ", summary_path)
  cat(readLines(summary_path), sep = "\n"); cat("\n")
}

# -------- main --------------------------------------------------------

if (opt$mode == "single") {
  run_single(opt)
} else if (opt$mode == "batch") {
  run_batch(opt)
} else {
  stop("--mode must be 'single' or 'batch'")
}
