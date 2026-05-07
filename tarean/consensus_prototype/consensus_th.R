#!/usr/bin/env Rscript
# H3 — TideHunter-fragment consensus CLI driver.
#
# --mode batch  walk a TideHunter GFF3 + clustering GFF3, build per-TRA
#               consensus from fragments, write per_tra_th_consensus.fasta
#               + per_tra_th_diagnostics.tsv + summary.log.
# --mode single one TRA id, same inputs; emit consensus + diagnostics.
#
# Mirrors the CLI shape of consensus_msa.R so the same BLAST validator
# can consume the output FASTA.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(parallel)
  library(GenomicRanges)
})

.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- .args[grep("--file=", .args)]
.script_dir <- if (length(.file_arg) > 0)
  normalizePath(dirname(sub("--file=", "", .file_arg))) else getwd()
source(file.path(.script_dir, "lib", "consensus_msa_core.R"))   # run_mafft
source(file.path(.script_dir, "lib", "consensus_th_core.R"))    # compute_tra_consensus_th

option_list <- list(
  make_option("--mode",            type = "character", default = "batch"),
  make_option("--tidehunter-gff3", type = "character", default = NULL),
  make_option("--clustering-gff3", type = "character", default = NULL),
  make_option("--tra-id",          type = "character", default = NULL,
              help = "single mode: TRC_X__seqid__start-end"),
  make_option("--out",             type = "character", default = NULL,
              help = "single mode: output prefix"),
  make_option("--out-dir",         type = "character", default = NULL,
              help = "batch mode: output directory"),
  make_option("--cpu",             type = "integer", default = 4L),
  make_option("--min-cn",                type = "integer", default = 2L),
  make_option("--min-relative-length",   type = "double",  default = 0.5),
  make_option("--length-tol",            type = "double",  default = 0.15),
  make_option("--gap-cutoff",            type = "double",  default = 0.5),
  make_option("--min-base-freq",         type = "double",  default = 0.4),
  make_option("--mafft-args",            type = "character",
              default = "--retree 1 --maxiterate 0 --nuc --thread 1"),
  make_option("--mafft-timeout",         type = "integer", default = 300L))
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
    method               = res$method,
    reason               = as.character(res$reason %||% ""),
    n_fragments          = res$n_fragments %||% 0L,
    n_fragments_used     = res$n_fragments_used %||% 0L,
    bucket_median_length = res$bucket_median_length %||% NA_integer_,
    consensus_length     = res$consensus_length %||% 0L,
    aln_length           = res$aln_length %||% NA_integer_,
    n_cols_kept          = res$n_cols_kept %||% 0L,
    n_cols_dropped       = res$n_cols_dropped %||% 0L,
    n_N                  = res$n_N %||% 0L,
    cons_mean_entropy    = fnum(res$cons_mean_entropy, 3),
    wall_sec             = fnum(res$wall_sec, 2))
  modifyList(row, extras)
}

# ---- TideHunter / clustering parsing -------------------------------

parse_attr <- function(a, k) {
  m <- regmatches(a, regexpr(paste0("(?<=", k, "=)[^;]+"), a, perl = TRUE))
  if (length(m) == 0L) NA_character_ else m
}

read_gff3 <- function(p) {
  d <- read.table(p, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                  comment.char = "#", quote = "")
  names(d) <- c("seqid", "src", "ftype", "start", "end", "score",
                "strand", "phase", "attrs")
  d
}

# Returns list:
#   tra:  data.frame with seqid, start, end, TRC, tra_id (1 row per TRA)
#   frag: data.frame with tra_id, seq, length, copy_number
build_tra_fragment_index <- function(th_path, cl_path) {
  th <- read_gff3(th_path)
  cl <- read_gff3(cl_path)
  th$cons_seq <- vapply(th$attrs, parse_attr, character(1), k = "consensus_sequence")
  th$cons_len <- as.integer(vapply(th$attrs, parse_attr, character(1),
                                   k = "consensus_length"))
  th$cn       <- as.numeric(vapply(th$attrs, parse_attr, character(1),
                                   k = "copy_number"))
  cl$TRC      <- vapply(cl$attrs, parse_attr, character(1), k = "Name")
  cl$tra_id   <- sprintf("%s__%s__%d-%d", cl$TRC, cl$seqid, cl$start, cl$end)

  th_gr <- GRanges(th$seqid, IRanges(th$start, th$end))
  cl_gr <- GRanges(cl$seqid, IRanges(cl$start, cl$end))
  ov <- findOverlaps(th_gr, cl_gr, minoverlap = 1L, type = "any")
  ovw <- pmin(th$end[queryHits(ov)],   cl$end[subjectHits(ov)]) -
         pmax(th$start[queryHits(ov)], cl$start[subjectHits(ov)]) + 1L
  df <- data.frame(th_i = queryHits(ov), cl_i = subjectHits(ov), w = ovw)
  df <- df[order(df$th_i, -df$w), ]
  df <- df[!duplicated(df$th_i), ]

  frag <- data.frame(
    tra_id      = cl$tra_id[df$cl_i],
    seq         = th$cons_seq[df$th_i],
    length      = th$cons_len[df$th_i],
    copy_number = th$cn[df$th_i],
    stringsAsFactors = FALSE)
  frag <- frag[!is.na(frag$seq) & nchar(frag$seq) > 0, ]
  list(tra  = cl[, c("seqid", "start", "end", "TRC", "tra_id")],
       frag = frag)
}

# ---- single mode ---------------------------------------------------

run_single <- function(opt) {
  stopifnot(!is.null(opt$`tidehunter-gff3`),
            !is.null(opt$`clustering-gff3`),
            !is.null(opt$`tra-id`),
            !is.null(opt$out))
  idx <- build_tra_fragment_index(opt$`tidehunter-gff3`,
                                  opt$`clustering-gff3`)
  fr <- idx$frag[idx$frag$tra_id == opt$`tra-id`, , drop = FALSE]
  if (nrow(fr) == 0L)
    stop("no TideHunter fragments mapped to ", opt$`tra-id`)
  res <- compute_tra_consensus_th(
    fr,
    min_cn = opt$`min-cn`,
    min_relative_length = opt$`min-relative-length`,
    length_tol = opt$`length-tol`,
    gap_cutoff = opt$`gap-cutoff`,
    min_base_freq = opt$`min-base-freq`,
    mafft_args = opt$`mafft-args`,
    mafft_timeout = opt$`mafft-timeout`)
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  writeLines(format_fasta(paste0(opt$`tra-id`, "  [", res$method, "]"),
                          res$consensus),
             con = paste0(opt$out, "_consensus.fasta"))
  write.table(
    as.data.frame(diag_row(res, list(id = opt$`tra-id`)),
                  stringsAsFactors = FALSE),
    file = paste0(opt$out, "_diagnostics.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("method=%s  reason=%s  consensus_length=%d  wall=%.1fs\n",
              res$method, as.character(res$reason),
              res$consensus_length, res$wall_sec))
}

# ---- batch mode ----------------------------------------------------

process_one <- function(tra_row, frag_by_tra, opt) {
  fr <- frag_by_tra[[tra_row$tra_id]]
  if (is.null(fr)) {
    return(list(ok = FALSE,
                reason = "no_fragments_mapped",
                tra_row = tra_row))
  }
  res <- compute_tra_consensus_th(
    fr,
    min_cn = opt$`min-cn`,
    min_relative_length = opt$`min-relative-length`,
    length_tol = opt$`length-tol`,
    gap_cutoff = opt$`gap-cutoff`,
    min_base_freq = opt$`min-base-freq`,
    mafft_args = opt$`mafft-args`,
    mafft_timeout = opt$`mafft-timeout`)
  diag <- diag_row(res, list(
    TRC_ID = tra_row$TRC,
    seqid  = tra_row$seqid,
    start  = tra_row$start,
    end    = tra_row$end,
    id     = tra_row$tra_id))
  list(ok = TRUE,
       consensus = res$consensus, id = tra_row$tra_id,
       diag = diag, method = res$method, reason = res$reason)
}

run_batch <- function(opt) {
  for (f in c("tidehunter-gff3", "clustering-gff3", "out-dir"))
    if (is.null(opt[[f]])) stop("batch requires --", f)
  dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

  t0 <- Sys.time()
  message("indexing TideHunter + clustering")
  idx <- build_tra_fragment_index(opt$`tidehunter-gff3`,
                                  opt$`clustering-gff3`)
  message(sprintf("  TRAs: %d ; fragments: %d", nrow(idx$tra), nrow(idx$frag)))
  frag_by_tra <- split(idx$frag, idx$frag$tra_id)
  tra_rows <- lapply(split(idx$tra, seq_len(nrow(idx$tra))), as.list)

  results <- mclapply(
    tra_rows,
    function(r) process_one(r, frag_by_tra = frag_by_tra, opt = opt),
    mc.cores = opt$cpu)

  ok_mask <- vapply(results, function(x) isTRUE(x$ok), logical(1))
  ok     <- results[ok_mask]
  failed <- results[!ok_mask]

  fa_path <- file.path(opt$`out-dir`, "per_tra_th_consensus.fasta")
  con <- file(fa_path, "w")
  for (r in ok) if (nchar(r$consensus) > 0L)
    cat(format_fasta(r$id, r$consensus), file = con)
  close(con)

  if (length(ok) > 0L) {
    diag_df <- do.call(rbind, lapply(ok, function(r)
      as.data.frame(r$diag, stringsAsFactors = FALSE)))
    write.table(diag_df,
                file = file.path(opt$`out-dir`, "per_tra_th_diagnostics.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (length(failed) > 0L) {
    df <- do.call(rbind, lapply(failed, function(x)
      data.frame(id     = x$tra_row$tra_id,
                 TRC_ID = x$tra_row$TRC,
                 seqid  = x$tra_row$seqid,
                 start  = x$tra_row$start,
                 end    = x$tra_row$end,
                 reason = x$reason,
                 stringsAsFactors = FALSE)))
    write.table(df, file = file.path(opt$`out-dir`, "failed.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  method_tab <- table(vapply(ok, function(x) x$method, character(1)))
  dt <- difftime(Sys.time(), t0, units = "secs")
  summary_path <- file.path(opt$`out-dir`, "summary.log")
  writeLines(c(
    sprintf("total TRAs:           %d", nrow(idx$tra)),
    sprintf("processed ok:         %d", length(ok)),
    sprintf("processed failed:     %d", length(failed)),
    sprintf("th_msa_ok:            %d", method_tab["th_msa_ok"]    %||% 0L),
    sprintf("th_single_ok:        %d", method_tab["th_single_ok"] %||% 0L),
    sprintf("th_failed:            %d", method_tab["th_failed"]    %||% 0L),
    sprintf("wall time (s):        %.1f", as.numeric(dt)),
    sprintf("output FASTA:         %s", fa_path)
  ), summary_path)
  message("summary: ", summary_path)
  cat(readLines(summary_path), sep = "\n"); cat("\n")
}

if (opt$mode == "single") run_single(opt) else
if (opt$mode == "batch")  run_batch(opt)  else
  stop("--mode must be 'single' or 'batch'")
