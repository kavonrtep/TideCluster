#!/usr/bin/env Rscript
# Consensus selector — emits the selected per-TRA consensus.
#
# Reads two pre-computed consensus FASTA files (array-MSA consensus from
# consensus_msa.R and TideHunter consensus from validate_th_single.R)
# plus their BLAST validator outputs, selects the better consensus per
# TRA by coverage_frac, and emits a unified FASTA + per-TRA metrics with
# extended variability diagnostics derived from the winning method's
# self-hits.
#
# No production wiring; prototype-only.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

option_list <- list(
  make_option("--mb-fasta",       type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/per_tra_consensus.fasta"),
  make_option("--mb-blast-tsv",   type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/blast_validation/blast.tsv"),
  make_option("--th-fasta",       type = "character",
              default = "tmp/consensus_prototype_drapa_th_single/per_tra_th_consensus.fasta"),
  make_option("--th-blast-tsv",   type = "character",
              default = "tmp/consensus_prototype_drapa_th_single/blast.tsv"),
  make_option("--kite-tsv",       type = "character",
              default = "tmp/drapa_run2/drapa_kite/monomer_size_top3_estimats.csv"),
  make_option("--out-dir",        type = "character",
              default = "tmp/consensus_prototype_drapa_ensemble"),
  # Quality grading thresholds (see grade assignment below).
  make_option("--cov-strict",        type = "double", default = 0.90,
              help = "minimum coverage_frac for grade A (clean pass)"),
  make_option("--cov-qualified",     type = "double", default = 0.80,
              help = "minimum coverage_frac for grade B (qualified pass)"),
  make_option("--cov-low",           type = "double", default = 0.50,
              help = "minimum coverage_frac for grade C (low coverage)"),
  make_option("--pident-sd-cutoff",  type = "double", default = 1.5,
              help = "pident_sd >= this raises flag_heterogeneous"),
  make_option("--gap-cutoff-bp",     type = "integer", default = 1000L,
              help = "max_internal_gap_bp >= this raises flag_internal_gap"),
  make_option("--core-cutoff",       type = "double", default = 0.95,
              help = "core_fraction < this AND core_coverage >= cov-strict raises flag_boundary_overext"),
  make_option("--pident-p05-cutoff", type = "double", default = 92.0,
              help = "pident_p05 < this raises flag_low_pident")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# ---- BLAST parsing helpers (shared with validators) -----------------

parse_q <- function(qseqid) {
  parts <- strsplit(qseqid, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(list(seqid = NA, start = NA, end = NA))
  r <- strsplit(parts[3], "-", fixed = TRUE)[[1]]
  list(seqid = parts[2], start = as.integer(r[1]), end = as.integer(r[2]))
}
parse_s <- function(sseqid) {
  parts <- strsplit(sseqid, "_", fixed = TRUE)[[1]]; n <- length(parts)
  if (n < 3) return(list(seqid = NA, start = NA, end = NA))
  list(seqid = paste(parts[seq_len(n - 2)], collapse = "_"),
       start = suppressWarnings(as.integer(parts[n - 1])),
       end   = suppressWarnings(as.integer(parts[n])))
}

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) return(list(covered = 0L, gaps = integer(0)))
  lo <- pmin(starts, ends); hi <- pmax(starts, ends)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  cur_lo <- lo[1]; cur_hi <- hi[1]; covered <- 0L; gaps <- integer(0)
  for (i in seq_along(lo)[-1]) {
    if (lo[i] <= cur_hi + 1L) {
      cur_hi <- max(cur_hi, hi[i])
    } else {
      covered <- covered + (cur_hi - cur_lo + 1L)
      gaps    <- c(gaps, lo[i] - cur_hi - 1L)
      cur_lo  <- lo[i]; cur_hi <- hi[i]
    }
  }
  covered <- covered + (cur_hi - cur_lo + 1L)
  list(covered = covered, gaps = gaps)
}

read_self_hits <- function(path) {
  message("reading ", path)
  bl <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("qseqid","sseqid","pident","length",
                                 "mismatch","gapopen","qstart","qend",
                                 "sstart","send","evalue","bitscore",
                                 "qlen","slen"))
  uq <- unique(bl$qseqid); us <- unique(bl$sseqid)
  qtab <- do.call(rbind, lapply(uq, function(x)
    as.data.frame(parse_q(x), stringsAsFactors = FALSE))); rownames(qtab) <- uq
  stab <- do.call(rbind, lapply(us, function(x)
    as.data.frame(parse_s(x), stringsAsFactors = FALSE))); rownames(stab) <- us
  bl$q_seqid <- qtab[bl$qseqid, "seqid"]; bl$q_start <- qtab[bl$qseqid, "start"]; bl$q_end <- qtab[bl$qseqid, "end"]
  bl$s_seqid <- stab[bl$sseqid, "seqid"]; bl$s_start <- stab[bl$sseqid, "start"]; bl$s_end <- stab[bl$sseqid, "end"]
  bl[bl$q_seqid == bl$s_seqid & bl$q_start == bl$s_start & bl$q_end == bl$s_end, , drop = FALSE]
}

# ---- per-TRA metrics with extended variability diagnostics ----------

compute_metrics <- function(self_hits, source_label) {
  by_q <- split(self_hits, self_hits$qseqid)
  do.call(rbind, lapply(names(by_q), function(qid) {
    sub  <- by_q[[qid]]
    iv   <- merge_intervals(sub$sstart, sub$send)
    slen <- sub$slen[1]; m <- sub$qlen[1]
    cov_full <- iv$covered / slen
    core_lo  <- min(pmin(sub$sstart, sub$send))
    core_hi  <- max(pmax(sub$sstart, sub$send))
    core_span <- core_hi - core_lo + 1L
    core_cov  <- iv$covered / core_span
    pident_w <- sum(sub$pident * sub$length) / sum(sub$length)
    pident_sd <- if (nrow(sub) > 1L) sd(sub$pident) else NA_real_
    pident_q  <- as.numeric(quantile(sub$pident, c(0.05, 0.5, 0.95), na.rm = TRUE))
    n_perfect <- sum(sub$pident >= 99.5)
    expected  <- slen / m
    hit_ratio <- nrow(sub) / expected
    data.frame(
      id              = qid,
      source          = source_label,
      array_length    = slen,
      cons_length     = m,
      hit_count       = nrow(sub),
      coverage_frac   = cov_full,
      core_coverage   = core_cov,
      core_fraction   = core_span / slen,
      flank_bp        = slen - core_span,
      mean_pident     = pident_w,
      pident_sd       = pident_sd,
      pident_p05      = pident_q[1],
      pident_p50      = pident_q[2],
      pident_p95      = pident_q[3],
      n_perfect_hits  = n_perfect,
      perfect_hit_frac = n_perfect / nrow(sub),
      max_internal_gap_bp = if (length(iv$gaps) > 0L) max(iv$gaps) else 0L,
      n_internal_gaps = length(iv$gaps),
      hit_ratio       = hit_ratio,
      stringsAsFactors = FALSE)
  }))
}

# ---- run -----------------------------------------------------------

mb_h <- read_self_hits(opt$`mb-blast-tsv`)
th_h <- read_self_hits(opt$`th-blast-tsv`)
mb_m <- compute_metrics(mb_h, "M_B")
th_m <- compute_metrics(th_h, "TH_single")
message(sprintf("M_B per-TRA rows:        %d", nrow(mb_m)))
message(sprintf("TH-single per-TRA rows:  %d", nrow(th_m)))

all_ids <- union(mb_m$id, th_m$id)
both <- merge(mb_m, th_m, by = "id", all = TRUE,
              suffixes = c("_mb", "_th"))

# Picker: choose the consensus with the higher coverage_frac.
# Ties broken by the higher mean_pident.
both$cov_mb_x <- ifelse(is.na(both$coverage_frac_mb), -Inf, both$coverage_frac_mb)
both$cov_th_x <- ifelse(is.na(both$coverage_frac_th), -Inf, both$coverage_frac_th)
both$pid_mb_x <- ifelse(is.na(both$mean_pident_mb),    0,    both$mean_pident_mb)
both$pid_th_x <- ifelse(is.na(both$mean_pident_th),    0,    both$mean_pident_th)
both$pick <- ifelse(
  abs(both$cov_mb_x - both$cov_th_x) < 1e-6,
  ifelse(both$pid_mb_x >= both$pid_th_x, "M_B", "TH_single"),
  ifelse(both$cov_mb_x >  both$cov_th_x, "M_B", "TH_single"))

# Pick the source's metrics
pick_col <- function(col) {
  ifelse(both$pick == "M_B",
         both[[paste0(col, "_mb")]],
         both[[paste0(col, "_th")]])
}

result <- data.frame(
  id              = both$id,
  source          = both$pick,
  array_length    = pick_col("array_length"),
  cons_length     = pick_col("cons_length"),
  hit_count       = pick_col("hit_count"),
  coverage_frac   = pick_col("coverage_frac"),
  core_coverage   = pick_col("core_coverage"),
  core_fraction   = pick_col("core_fraction"),
  flank_bp        = pick_col("flank_bp"),
  mean_pident     = pick_col("mean_pident"),
  pident_sd       = pick_col("pident_sd"),
  pident_p05      = pick_col("pident_p05"),
  pident_p50      = pick_col("pident_p50"),
  pident_p95      = pick_col("pident_p95"),
  n_perfect_hits  = pick_col("n_perfect_hits"),
  perfect_hit_frac= pick_col("perfect_hit_frac"),
  max_internal_gap_bp = pick_col("max_internal_gap_bp"),
  n_internal_gaps = pick_col("n_internal_gaps"),
  hit_ratio       = pick_col("hit_ratio"),
  cov_mb          = both$coverage_frac_mb,
  cov_th          = both$coverage_frac_th,
  stringsAsFactors = FALSE)

# ---- quality grade + diagnostic flags --------------------------------
#
# `quality_grade` summarises usability in a single label; `flags` is a
# comma-separated list of independent reasons. Flags can co-occur with
# any grade (a clean-coverage TRA can still be heterogeneous, etc.).

flag_boundary_overext <- !is.na(result$core_coverage) &
  !is.na(result$core_fraction) &
  result$core_coverage >= opt$`cov-strict` &
  result$core_fraction <  opt$`core-cutoff`
flag_internal_gap     <- !is.na(result$max_internal_gap_bp) &
  result$max_internal_gap_bp >= opt$`gap-cutoff-bp`
flag_heterogeneous    <- !is.na(result$pident_sd) &
  result$pident_sd >= opt$`pident-sd-cutoff`
flag_low_pident       <- !is.na(result$pident_p05) &
  result$pident_p05 <  opt$`pident-p05-cutoff`

result$flag_boundary_overext <- flag_boundary_overext
result$flag_internal_gap     <- flag_internal_gap
result$flag_heterogeneous    <- flag_heterogeneous
result$flag_low_pident       <- flag_low_pident

flag_names <- c("boundary_overext", "internal_gap",
                "heterogeneous", "low_pident")
flag_mat <- cbind(flag_boundary_overext, flag_internal_gap,
                  flag_heterogeneous, flag_low_pident)
result$flags <- vapply(seq_len(nrow(result)), function(i) {
  set <- flag_names[flag_mat[i, ]]
  if (length(set) == 0L) "" else paste(set, collapse = ",")
}, character(1))

# Grade is a coverage tier; flags are independent diagnostics. A grade-A
# TRA can still carry flags (e.g. a HOR array that tiles fully but is
# heterogeneous). A grade-B TRA almost always has at least one flag
# explaining why coverage is below the strict cutoff.
cov <- result$coverage_frac
result$quality_grade <- ifelse(is.na(cov),                       "D",
                        ifelse(cov >= opt$`cov-strict`,          "A",
                        ifelse(cov >= opt$`cov-qualified`,       "B",
                        ifelse(cov >= opt$`cov-low`,             "C",
                                                                  "D"))))

# Attach KITE annotation (HOR_status, monomer_size, TRC)
kite <- read.table(opt$`kite-tsv`, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE,
                   quote = "", comment.char = "")
kite$id <- sprintf("%s__%s__%d-%d", kite$TRC_ID, kite$seqid, kite$start, kite$end)
result <- merge(result, kite[, c("id","TRC_ID","monomer_size",
                                 "HOR_status","HOR_confidence")],
                by = "id", all.x = TRUE)

# Order by source then TRC then id
result <- result[order(result$source, result$TRC_ID, result$id), ]

write.table(result, file.path(opt$`out-dir`, "per_tra_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ---- write the picked-FASTA file ------------------------------------

mb_fa <- readDNAStringSet(opt$`mb-fasta`)
th_fa <- readDNAStringSet(opt$`th-fasta`)
out_fa <- file.path(opt$`out-dir`, "per_tra_consensus.fasta")
con <- file(out_fa, "w")
n_written <- 0L
n_missing <- 0L
for (i in seq_len(nrow(result))) {
  id <- result$id[i]
  src <- result$source[i]
  set <- if (src == "M_B") mb_fa else th_fa
  if (id %in% names(set)) {
    seq <- as.character(set[[id]])
    cat(">", id, " [", src, "]\n",
        gsub("(.{60})", "\\1\n", seq, perl = TRUE), "\n",
        sep = "", file = con)
    n_written <- n_written + 1L
  } else {
    n_missing <- n_missing + 1L
  }
}
close(con)
message(sprintf("wrote %d records to %s (%d missing)",
                n_written, out_fa, n_missing))

# ---- summary -------------------------------------------------------

n_total <- 2073
sumlines <- character()
sumlines <- c(sumlines,
  sprintf("total TRAs:            %d", n_total),
  sprintf("metrics rows:          %d", nrow(result)),
  sprintf("source = M_B:          %d", sum(result$source == "M_B")),
  sprintf("source = TH_single:    %d", sum(result$source == "TH_single")))

for (cut in c(0.99, 0.95, 0.90)) {
  n_pass <- sum(result$coverage_frac >= cut, na.rm = TRUE)
  sumlines <- c(sumlines,
    sprintf(">= %.0f%% coverage:      %4d (%.1f%%)",
            100*cut, n_pass, 100*n_pass/n_total))
}
sumlines <- c(sumlines,
  sprintf("median pident_sd:      %.2f", median(result$pident_sd, na.rm = TRUE)),
  sprintf("median perfect_frac:   %.2f", median(result$perfect_hit_frac, na.rm = TRUE)))

# Grade and flag counts (over all TRAs in the pipeline, including
# missing-from-metrics ones counted as grade D).
sumlines <- c(sumlines,
  "",
  "thresholds:",
  sprintf("  cov-strict       = %.2f", opt$`cov-strict`),
  sprintf("  cov-qualified    = %.2f", opt$`cov-qualified`),
  sprintf("  cov-low          = %.2f", opt$`cov-low`),
  sprintf("  pident-sd-cutoff = %.2f", opt$`pident-sd-cutoff`),
  sprintf("  gap-cutoff-bp    = %d",   opt$`gap-cutoff-bp`),
  sprintf("  core-cutoff      = %.2f", opt$`core-cutoff`),
  sprintf("  pident-p05-cutoff= %.1f", opt$`pident-p05-cutoff`),
  "",
  "quality_grade (over n_total):")
n_missing <- n_total - nrow(result)
grade_tab <- table(result$quality_grade)
for (g in c("A", "B", "C", "D")) {
  base <- if (g %in% names(grade_tab)) as.integer(grade_tab[g]) else 0L
  n <- base + (if (g == "D") n_missing else 0L)
  sumlines <- c(sumlines,
    sprintf("  %s: %4d (%.1f%%)", g, n, 100*n/n_total))
}
sumlines <- c(sumlines,
  "",
  "flags (over n_total; co-occurrence allowed):",
  sprintf("  boundary_overext: %4d", sum(result$flag_boundary_overext)),
  sprintf("  internal_gap:     %4d", sum(result$flag_internal_gap)),
  sprintf("  heterogeneous:    %4d", sum(result$flag_heterogeneous)),
  sprintf("  low_pident:       %4d", sum(result$flag_low_pident)))

writeLines(sumlines, file.path(opt$`out-dir`, "summary.log"))
cat(paste(sumlines, collapse = "\n"), "\n")

# Source x HOR status
cat("\nsource x HOR_status:\n")
print(table(result$source, result$HOR_status, useNA = "ifany"))

# Coverage bins by source
cat("\ncoverage_frac >=0.90 by source:\n")
for (s in c("M_B", "TH_single")) {
  sub <- result[result$source == s, ]
  cat(sprintf("  %-12s %d / %d (%.1f%%)\n", s,
              sum(sub$coverage_frac >= 0.90, na.rm = TRUE), nrow(sub),
              100 * mean(sub$coverage_frac >= 0.90, na.rm = TRUE)))
}

# Variability profile per HOR class
cat("\nvariability metrics by HOR_status (median across TRAs):\n")
for (h in unique(na.omit(result$HOR_status))) {
  sub <- result[!is.na(result$HOR_status) & result$HOR_status == h, ]
  cat(sprintf("  %-14s n=%4d  pident_sd=%.2f  perfect_frac=%.2f  pident_p05=%.1f  cov=%.3f\n",
              h, nrow(sub),
              median(sub$pident_sd, na.rm = TRUE),
              median(sub$perfect_hit_frac, na.rm = TRUE),
              median(sub$pident_p05, na.rm = TRUE),
              median(sub$coverage_frac, na.rm = TRUE)))
}
