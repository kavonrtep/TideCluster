#!/usr/bin/env Rscript
# Core-coverage diagnostic.
#
# The "core repeat region" of a TRA is the genomic span between the
# leftmost and rightmost self-BLAST hit (= the portion actually
# containing the repeat). Coverage measured against the core span is
# tolerant to TideHunter boundary overextension.
#
# For every TRA with at least one self-hit from validate_drapa_blast,
# compute:
#   cov_full      = covered_bp / slen                       (primary metric)
#   core_span     = max(send) - min(sstart) + 1             (hit footprint)
#   core_coverage = covered_bp / core_span                  (coverage inside core)
#   core_fraction = core_span / slen                        (fraction of TRA that is core)
#   flank_bp      = slen - core_span                        (unhit flanks)
#
# Rationale: if TideHunter over-extended the TRA boundary, hits cluster
# in the middle and cov_full is depressed even when the consensus is
# fine. In that case core_fraction < 1 and core_coverage ~ 1.
#
# No production changes; writes core_metrics.tsv next to blast.tsv.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--blast-tsv",  type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/blast_validation/blast.tsv"),
  make_option("--failing-tsv", type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/failing_tras_annotated.tsv"),
  make_option("--out",        type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/core_metrics.tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- parsing helpers copied from validate_drapa_blast.R -------------
parse_q <- function(qseqid) {
  parts <- strsplit(qseqid, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(list(seqid = NA, start = NA, end = NA))
  r <- strsplit(parts[3], "-", fixed = TRUE)[[1]]
  list(seqid = parts[2],
       start = as.integer(r[1]), end = as.integer(r[2]))
}
parse_s <- function(sseqid) {
  parts <- strsplit(sseqid, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  if (n < 3) return(list(seqid = NA, start = NA, end = NA))
  list(seqid = paste(parts[seq_len(n - 2)], collapse = "_"),
       start = suppressWarnings(as.integer(parts[n - 1])),
       end   = suppressWarnings(as.integer(parts[n])))
}
merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) return(0L)
  lo <- pmin(starts, ends); hi <- pmax(starts, ends)
  o  <- order(lo); lo <- lo[o]; hi <- hi[o]
  cur_lo <- lo[1]; cur_hi <- hi[1]; covered <- 0L
  for (i in seq_along(lo)[-1]) {
    if (lo[i] <= cur_hi + 1L) cur_hi <- max(cur_hi, hi[i])
    else { covered <- covered + (cur_hi - cur_lo + 1L)
           cur_lo <- lo[i]; cur_hi <- hi[i] }
  }
  covered + (cur_hi - cur_lo + 1L)
}

# ---- load + self-hit filter -----------------------------------------
message("reading ", opt$`blast-tsv`)
bl <- read.table(opt$`blast-tsv`, sep = "\t", header = FALSE,
                 stringsAsFactors = FALSE,
                 col.names = c("qseqid", "sseqid", "pident", "length",
                               "mismatch", "gapopen", "qstart", "qend",
                               "sstart", "send", "evalue", "bitscore",
                               "qlen", "slen"))
message("  raw hits: ", nrow(bl))
uq <- unique(bl$qseqid); us <- unique(bl$sseqid)
qtab <- do.call(rbind, lapply(uq, function(x) as.data.frame(parse_q(x), stringsAsFactors = FALSE)))
rownames(qtab) <- uq
stab <- do.call(rbind, lapply(us, function(x) as.data.frame(parse_s(x), stringsAsFactors = FALSE)))
rownames(stab) <- us
bl$q_seqid <- qtab[bl$qseqid, "seqid"]; bl$q_start <- qtab[bl$qseqid, "start"]; bl$q_end <- qtab[bl$qseqid, "end"]
bl$s_seqid <- stab[bl$sseqid, "seqid"]; bl$s_start <- stab[bl$sseqid, "start"]; bl$s_end <- stab[bl$sseqid, "end"]
self <- bl$q_seqid == bl$s_seqid & bl$q_start == bl$s_start & bl$q_end == bl$s_end
bl <- bl[self, , drop = FALSE]
message("  self-hits: ", nrow(bl))

# ---- per-TRA core-coverage metrics ----------------------------------
by_q <- split(bl, bl$qseqid)
rows <- lapply(names(by_q), function(qid) {
  sub <- by_q[[qid]]
  slen <- sub$slen[1]
  lo   <- pmin(sub$sstart, sub$send)
  hi   <- pmax(sub$sstart, sub$send)
  core_start <- min(lo); core_end <- max(hi)
  core_span  <- core_end - core_start + 1L
  covered    <- merge_intervals(sub$sstart, sub$send)
  data.frame(
    id              = qid,
    slen            = slen,
    covered_bp      = covered,
    cov_full        = covered / slen,
    core_start      = core_start,
    core_end        = core_end,
    core_span       = core_span,
    core_fraction   = core_span / slen,
    core_coverage   = covered / core_span,
    flank_bp        = slen - core_span,
    flank_left_bp   = core_start - 1L,
    flank_right_bp  = slen - core_end,
    hit_count       = nrow(sub),
    stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)

write.table(res, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", opt$out, " (", nrow(res), " TRAs)")

# ---- sizing the effect on failing TRAs ------------------------------
if (file.exists(opt$`failing-tsv`)) {
  f <- read.table(opt$`failing-tsv`, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
  m <- merge(f, res, by = "id", all.x = TRUE)

  cat("\n=== failing (n=", nrow(f), ") — core-coverage metric ===\n", sep = "")
  cat("n with any self-hit: ", sum(!is.na(m$core_coverage)), "\n", sep = "")

  # How many are rescued at the 90% cutoff by switching to core coverage?
  full_pass  <- !is.na(m$cov_full)        & m$cov_full        >= 0.90
  core_pass  <- !is.na(m$core_coverage)   & m$core_coverage   >= 0.90
  rescued    <- !full_pass & core_pass
  lost       <- full_pass  & !core_pass
  cat(sprintf("  FAIL(cov_full) -> PASS(core_coverage): %d\n", sum(rescued)))
  cat(sprintf("  PASS(cov_full) -> FAIL(core_coverage): %d (sanity: 0 expected)\n",
              sum(lost)))

  cat("\n=== boundary-overextension signal (core_coverage >= 0.90 AND core_fraction < 0.95) ===\n")
  bo <- !is.na(m$core_coverage) & m$core_coverage >= 0.90 & m$core_fraction < 0.95
  cat(sprintf("  n: %d\n", sum(bo)))
  cat("  flank_bp stats (bp of array outside the core repeat region):\n")
  print(summary(m$flank_bp[bo]))
  cat("  core_fraction stats:\n")
  print(summary(m$core_fraction[bo]))

  cat("\n=== by category (n rescued / n total) ===\n")
  cats <- sort(unique(m$category))
  for (cn in cats) {
    idx    <- m$category == cn
    n_tot  <- sum(idx, na.rm = TRUE)
    n_resc <- sum(rescued[idx], na.rm = TRUE)
    cat(sprintf("  %-22s %4d / %4d rescued\n", cn, n_resc, n_tot))
  }

  cat("\n=== example rescued TRAs (top 10 by flank_bp) ===\n")
  r <- m[rescued, ]
  r <- r[order(-r$flank_bp), ]
  show_cols <- c("TRC_ID","id","m","length_bp","cov_full","core_coverage",
                 "core_fraction","flank_bp","flank_left_bp","flank_right_bp",
                 "HOR_status")
  print(head(r[, show_cols], 10), row.names = FALSE)

  cat("\n=== stuck TRAs (low core_coverage) — sub-repeat / internal gap candidates ===\n")
  stuck <- !is.na(m$core_coverage) & m$core_coverage < 0.90
  s <- m[stuck, ]
  s <- s[order(s$core_coverage), ]
  cat(sprintf("  n stuck: %d\n", nrow(s)))
  cat("  core_coverage dist:\n"); print(summary(s$core_coverage))
  cat("  top 10 most-gapped inside core:\n")
  print(head(s[, c("TRC_ID","id","m","length_bp","cov_full","core_coverage",
                   "core_fraction","category","HOR_status")], 10),
        row.names = FALSE)

  # Write a merged report
  out_merged <- sub("\\.tsv$", "_with_core.tsv", opt$`failing-tsv`)
  merged_cols <- c(colnames(f),
                   "core_start","core_end","core_span",
                   "core_fraction","core_coverage","flank_bp",
                   "flank_left_bp","flank_right_bp")
  merged_cols <- intersect(merged_cols, colnames(m))
  write.table(m[, merged_cols], out_merged, sep = "\t", quote = FALSE,
              row.names = FALSE)
  message("\nwrote ", out_merged)
}
