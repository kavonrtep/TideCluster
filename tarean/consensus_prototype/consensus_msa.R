#!/usr/bin/env Rscript
# MSA-based per-TRA consensus prototype (M_B).
#
# CLI driver; parallels consensus_pfm.R so the BLAST validator can
# consume either. Two modes:
#   --mode single  one TRA FASTA + m, write consensus + diagnostics.
#   --mode batch   walk a KITE top-3 TSV + TRC fasta dir, emit one
#                  combined per_tra_consensus.fasta and per_tra_diagnostics.tsv.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(parallel)
})

.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .args[grep("--file=", .args)]
.script_dir <- if (length(.file_arg) > 0) {
  normalizePath(dirname(sub("--file=", "", .file_arg)))
} else getwd()
source(file.path(.script_dir, "lib", "consensus_core.R"))      # parse_record_name
source(file.path(.script_dir, "lib", "consensus_msa_core.R"))

option_list <- list(
  make_option("--mode",       type = "character", default = "batch"),
  make_option("--fasta",      type = "character", default = NULL),
  make_option("--m",          type = "integer",   default = NULL),
  make_option("--out",        type = "character", default = NULL),
  make_option("--kite-tsv",   type = "character", default = NULL),
  make_option("--tarean-dir", type = "character", default = NULL),
  make_option("--out-dir",    type = "character", default = NULL),
  make_option("--cpu",        type = "integer",   default = 4L),
  make_option("--k",                type = "integer", default = 8L),
  make_option("--max-monomers",     type = "integer", default = 150L),
  make_option("--min-anchor-count", type = "integer", default = 20L),
  make_option("--gap-cutoff",       type = "double",  default = 0.5),
  make_option("--min-base-freq",    type = "double",  default = 0.4),
  make_option("--mafft-args",       type = "character",
              default = "--retree 1 --maxiterate 0 --nuc --thread 1"),
  make_option("--mafft-timeout",    type = "integer", default = 300L))
opt <- parse_args(OptionParser(option_list = option_list))

format_fasta <- function(header, seq, width = 60L) {
  if (nchar(seq) == 0L) return(paste0(">", header, "\n\n"))
  wrapped <- gsub(sprintf("(.{%d})", width), "\\1\n", seq, perl = TRUE)
  paste0(">", header, "\n", wrapped, "\n")
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

diag_row <- function(res, extras) {
  fnum <- function(x, digits = 4)
    if (is.null(x) || (length(x) == 1 && is.na(x))) ""
    else sprintf(paste0("%.", digits, "f"), x)
  row <- list(
    method             = res$method,
    reason             = as.character(res$reason %||% ""),
    length_bp          = res$length_bp,
    m                  = res$m,
    anchor_kmer        = res$anchor_kmer %||% "",
    anchor_score       = fnum(res$anchor_score, 3),
    anchor_count       = res$anchor_count %||% 0L,
    n_monomers_extracted = res$n_monomers_extracted %||% 0L,
    n_monomers_aligned   = res$n_monomers_aligned %||% 0L,
    aln_length         = res$aln_length %||% 0L,
    consensus_length   = res$consensus_length %||% 0L,
    cons_mean_entropy  = fnum(res$cons_mean_entropy, 3),
    n_cols_kept        = res$n_cols_kept %||% 0L,
    n_cols_dropped     = res$n_cols_dropped %||% 0L,
    n_N                = res$n_N %||% 0L,
    wall_sec           = fnum(res$wall_sec, 2)
  )
  modifyList(row, extras)
}

# ---- single mode ----------------------------------------------------

run_single <- function(opt) {
  stopifnot(!is.null(opt$fasta), !is.null(opt$m), !is.null(opt$out))
  ds <- readDNAStringSet(opt$fasta)
  if (length(ds) == 0L) stop("empty FASTA")
  S <- ds[[1]]; rec_name <- names(ds)[1]
  res <- compute_tra_consensus_msa(
    S, opt$m, k = opt$k,
    max_monomers = opt$`max-monomers`,
    min_anchor_occurrences = opt$`min-anchor-count`,
    gap_cutoff = opt$`gap-cutoff`,
    min_base_freq = opt$`min-base-freq`,
    mafft_args = opt$`mafft-args`,
    mafft_timeout = opt$`mafft-timeout`)
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  writeLines(format_fasta(paste0(rec_name, "  m=", opt$m,
                                 "  [", res$method, "]"),
                          res$consensus),
             con = paste0(opt$out, "_consensus.fasta"))
  write.table(
    as.data.frame(diag_row(res, list(record = rec_name, m_cli = opt$m)),
                  stringsAsFactors = FALSE),
    file = paste0(opt$out, "_diagnostics.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("method=%s  reason=%s  consensus_length=%d  wall=%.1fs\n",
              res$method, as.character(res$reason),
              res$consensus_length, res$wall_sec))
}

# ---- batch mode -----------------------------------------------------

read_kite_tsv <- function(path) {
  read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

index_trc_fasta <- function(fa_path) {
  ds <- readDNAStringSet(fa_path)
  halves <- lapply(ds, function(x) subseq(x, start = 1L, end = nchar(x) %/% 2L))
  parsed <- lapply(names(ds), parse_record_name)
  list(halves = halves, parsed = parsed, names = names(ds))
}

process_one <- function(row, trc_index_cache, opt) {
  trc <- row$TRC_ID
  if (!trc %in% names(trc_index_cache)) {
    return(list(ok = FALSE, reason = paste0("no fasta for ", trc), row = row))
  }
  idx <- trc_index_cache[[trc]]
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
  res <- compute_tra_consensus_msa(
    S, m, k = opt$k,
    max_monomers = opt$`max-monomers`,
    min_anchor_occurrences = opt$`min-anchor-count`,
    gap_cutoff = opt$`gap-cutoff`,
    min_base_freq = opt$`min-base-freq`,
    mafft_args = opt$`mafft-args`,
    mafft_timeout = opt$`mafft-timeout`)
  id <- sprintf("%s__%s__%d-%d", trc, row$seqid, row$start, row$end)
  diag <- diag_row(res, list(TRC_ID = trc, seqid = row$seqid,
                             start = row$start, end = row$end, id = id))
  list(ok = TRUE, consensus = res$consensus, id = id, diag = diag,
       method = res$method)
}

run_batch <- function(opt) {
  for (f in c("kite-tsv", "tarean-dir", "out-dir"))
    if (is.null(opt[[f]])) stop("batch requires --", f)
  dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

  t0 <- Sys.time()
  tsv <- read_kite_tsv(opt$`kite-tsv`)
  message(sprintf("TSV rows: %d ; unique TRCs: %d",
                  nrow(tsv), length(unique(tsv$TRC_ID))))

  unique_trcs <- unique(tsv$TRC_ID)
  trc_index_cache <- setNames(vector("list", length(unique_trcs)),
                              unique_trcs)
  for (trc in unique_trcs) {
    fa <- file.path(opt$`tarean-dir`, paste0(trc, ".fasta"))
    if (!file.exists(fa)) next
    trc_index_cache[[trc]] <- index_trc_fasta(fa)
  }
  present <- sum(!vapply(trc_index_cache, is.null, logical(1)))
  message(sprintf("TRC fastas found: %d / %d", present, length(unique_trcs)))

  rows <- lapply(split(tsv, seq_len(nrow(tsv))), as.list)
  results <- mclapply(
    rows,
    function(r) process_one(r, trc_index_cache = trc_index_cache, opt = opt),
    mc.cores = opt$cpu)

  ok_mask <- vapply(results, function(x) isTRUE(x$ok), logical(1))
  ok     <- results[ok_mask]
  failed <- results[!ok_mask]

  fa_path <- file.path(opt$`out-dir`, "per_tra_consensus.fasta")
  con <- file(fa_path, "w")
  for (r in ok) if (nchar(r$consensus) > 0L)
    cat(format_fasta(r$id, r$consensus), file = con)
  close(con)

  if (length(ok) > 0L) {
    diag_df <- do.call(rbind, lapply(ok, function(r)
      as.data.frame(r$diag, stringsAsFactors = FALSE)))
    write.table(diag_df,
                file = file.path(opt$`out-dir`, "per_tra_diagnostics.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (length(failed) > 0L) {
    df <- do.call(rbind, lapply(failed, function(x)
      data.frame(TRC_ID = x$row$TRC_ID, seqid = x$row$seqid,
                 start = x$row$start, end = x$row$end,
                 m = x$row$monomer_size, reason = x$reason,
                 stringsAsFactors = FALSE)))
    write.table(df, file = file.path(opt$`out-dir`, "failed.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  method_tab <- table(vapply(ok, function(x) x$method, character(1)))
  dt <- difftime(Sys.time(), t0, units = "secs")
  summary_path <- file.path(opt$`out-dir`, "summary.log")
  writeLines(c(
    sprintf("total rows:            %d", nrow(tsv)),
    sprintf("processed ok:          %d", length(ok)),
    sprintf("processed failed:      %d", length(failed)),
    sprintf("msa_ok:                %d", method_tab["msa_ok"] %||% 0L),
    sprintf("msa_failed:            %d",
            sum(method_tab[names(method_tab) != "msa_ok"])),
    sprintf("wall time (s):         %.1f", as.numeric(dt)),
    sprintf("output FASTA:          %s", fa_path)
  ), summary_path)
  message("summary: ", summary_path)
  cat(readLines(summary_path), sep = "\n"); cat("\n")
}

if (opt$mode == "single") run_single(opt) else
if (opt$mode == "batch")  run_batch(opt)  else
  stop("--mode must be 'single' or 'batch'")
