#!/usr/bin/env Rscript
# BLAST the per-TRA consensus FASTA against each consensus's source
# tandem repeat array and report how well the consensus matches the
# array that produced it. Independent of TAREAN's consensus.
#
# Pipeline:
#   1. Read all drapa_tarean/fasta/TRC_*.fasta, take first half of each
#      record (KITE's analysis sequence), concatenate into
#      all_arrays.fasta.
#   2. makeblastdb on that flat file.
#   3. blastn per_tra_consensus.fasta vs the DB at strict settings.
#   4. Parse tabular output; for each consensus, keep only hits to its
#      own TRA (matched by (seqid, start, end)); compute:
#        hit_count          ~ array_length / m
#        total_align_bp     sum of hit alignment lengths
#        covered_bp         merged hit footprint on the array
#        coverage_fraction  covered_bp / array_length
#        mean_pident        length-weighted percent identity
#        strand_plus        fraction of hits on the plus strand
#        expected_hits      array_length / m
#        hit_ratio          hit_count / expected_hits
#   5. Aggregate distributions; write a REPORT.md and per-TRA TSV.
#
# Run from repo root:
#   Rscript tarean/consensus_prototype/tests/validate_drapa_blast.R

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
})

option_list <- list(
  make_option("--consensus-fasta", type = "character",
              default = "tmp/consensus_prototype_drapa/per_tra_consensus.fasta"),
  make_option("--diagnostics-tsv", type = "character",
              default = "tmp/consensus_prototype_drapa/per_tra_diagnostics.tsv"),
  make_option("--tarean-dir",      type = "character",
              default = "tmp/drapa_run2/drapa_tarean/fasta"),
  make_option("--kite-tsv",        type = "character",
              default = "tmp/drapa_run2/drapa_kite/monomer_size_top3_estimats.csv"),
  make_option("--out-dir",         type = "character",
              default = "tmp/consensus_prototype_drapa/blast_validation"),
  make_option("--pident",          type = "integer", default = 85L),
  make_option("--evalue",          type = "character", default = "1e-10"),
  make_option("--cpu",             type = "integer", default = 4L),
  make_option("--force",           type = "logical", default = FALSE,
              action = "store_true", help = "rebuild cached intermediate files")
)
opt <- parse_args(OptionParser(option_list = option_list))

repo_root <- normalizePath(getwd())
resolve_path <- function(p) {
  if (startsWith(p, "/")) p
  else normalizePath(file.path(repo_root, p), mustWork = FALSE)
}
opt$`consensus-fasta` <- resolve_path(opt$`consensus-fasta`)
opt$`diagnostics-tsv` <- resolve_path(opt$`diagnostics-tsv`)
opt$`tarean-dir`      <- resolve_path(opt$`tarean-dir`)
opt$`kite-tsv`        <- resolve_path(opt$`kite-tsv`)
opt$`out-dir`         <- resolve_path(opt$`out-dir`)

dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. Flat DB of arrays ----------------------------------

db_fasta <- file.path(opt$`out-dir`, "all_arrays.fasta")
if (opt$force || !file.exists(db_fasta)) {
  message("building flat array FASTA")
  trc_files <- list.files(opt$`tarean-dir`, pattern = "\\.fasta$", full.names = TRUE)
  con <- file(db_fasta, "w")
  total_records <- 0L
  total_bp <- 0
  for (f in trc_files) {
    ds <- readDNAStringSet(f)
    for (i in seq_along(ds)) {
      h <- subseq(ds[[i]], start = 1L, end = nchar(ds[[i]]) %/% 2L)
      cat(">", names(ds)[i], "\n", as.character(h), "\n",
          sep = "", file = con)
      total_bp <- total_bp + length(h)
      total_records <- total_records + 1L
    }
  }
  close(con)
  message(sprintf("  wrote %d array records (total %.1f Mb)",
                  total_records, total_bp / 1e6))
}

# ---------- 2. BLAST DB ---------------------------------------------

db_prefix <- file.path(opt$`out-dir`, "all_arrays")
if (opt$force || !file.exists(paste0(db_prefix, ".nhr"))) {
  message("building BLAST DB")
  cmd <- sprintf("makeblastdb -in %s -dbtype nucl -out %s -logfile %s",
                 shQuote(db_fasta), shQuote(db_prefix),
                 shQuote(file.path(opt$`out-dir`, "makeblastdb.log")))
  if (system(cmd) != 0) stop("makeblastdb failed")
}

# ---------- 3. blastn -----------------------------------------------

blast_tsv <- file.path(opt$`out-dir`, "blast.tsv")
if (opt$force || !file.exists(blast_tsv)) {
  message(sprintf("running blastn (pident >= %d, evalue %s)",
                  opt$pident, opt$evalue))
  outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
  cmd <- sprintf(
    "blastn -task blastn -db %s -query %s -outfmt %s -perc_identity %d -evalue %s -num_threads %d -max_target_seqs 10000 -out %s",
    shQuote(db_prefix), shQuote(opt$`consensus-fasta`),
    outfmt, opt$pident, opt$evalue, opt$cpu, shQuote(blast_tsv))
  t0 <- Sys.time()
  if (system(cmd) != 0) stop("blastn failed")
  message(sprintf("  blastn wall time: %.1f s",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---------- 4. Parse and compute per-TRA metrics --------------------

message("parsing BLAST output")
bl <- read.table(blast_tsv, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("qseqid", "sseqid", "pident", "length",
                               "mismatch", "gapopen", "qstart", "qend",
                               "sstart", "send", "evalue", "bitscore",
                               "qlen", "slen"))
message(sprintf("  %d raw hits", nrow(bl)))

parse_q <- function(qseqid) {
  parts <- strsplit(qseqid, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(list(trc = NA, seqid = NA, start = NA, end = NA))
  r <- strsplit(parts[3], "-", fixed = TRUE)[[1]]
  list(trc = parts[1], seqid = parts[2],
       start = as.integer(r[1]), end = as.integer(r[2]))
}
parse_s <- function(sseqid) {
  parts <- strsplit(sseqid, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  if (n < 3) return(list(seqid = NA, start = NA, end = NA))
  end_   <- suppressWarnings(as.integer(parts[n]))
  start_ <- suppressWarnings(as.integer(parts[n - 1]))
  seqid  <- paste(parts[seq_len(n - 2)], collapse = "_")
  list(seqid = seqid, start = start_, end = end_)
}
uq <- unique(bl$qseqid); us <- unique(bl$sseqid)
qtab <- do.call(rbind, lapply(uq, function(x) as.data.frame(parse_q(x), stringsAsFactors = FALSE)))
rownames(qtab) <- uq
stab <- do.call(rbind, lapply(us, function(x) as.data.frame(parse_s(x), stringsAsFactors = FALSE)))
rownames(stab) <- us
bl$q_seqid <- qtab[bl$qseqid, "seqid"]
bl$q_start <- qtab[bl$qseqid, "start"]
bl$q_end   <- qtab[bl$qseqid, "end"]
bl$s_seqid <- stab[bl$sseqid, "seqid"]
bl$s_start <- stab[bl$sseqid, "start"]
bl$s_end   <- stab[bl$sseqid, "end"]

self <- bl$q_seqid == bl$s_seqid & bl$q_start == bl$s_start & bl$q_end == bl$s_end
bl_self <- bl[self, , drop = FALSE]
message(sprintf("  %d self-hits", nrow(bl_self)))

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) return(list(covered = 0L, n = 0L))
  lo <- pmin(starts, ends); hi <- pmax(starts, ends)
  o <- order(lo)
  lo <- lo[o]; hi <- hi[o]
  cur_lo <- lo[1]; cur_hi <- hi[1]
  covered <- 0L
  for (i in seq_along(lo)[-1]) {
    if (lo[i] <= cur_hi + 1L) {
      cur_hi <- max(cur_hi, hi[i])
    } else {
      covered <- covered + (cur_hi - cur_lo + 1L)
      cur_lo <- lo[i]; cur_hi <- hi[i]
    }
  }
  covered <- covered + (cur_hi - cur_lo + 1L)
  list(covered = covered)
}

kite <- read.table(opt$`kite-tsv`, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE, quote = "",
                   comment.char = "")
# Accept both CSV schemas. kitehor (TideCluster >=1.10) emits lowercase
# `hor_*` column names; the legacy R kite.R emitted capitalised `HOR_*`.
# The rest of this script was written for the capitalised names, so we
# alias the new lowercase columns up when the old ones are absent.
.kite_aliases <- list(
  c("hor_status",       "HOR_status"),
  c("hor_confidence",   "HOR_confidence"),
  c("hor_founder",      "HOR_base_monomer"),
  c("hor_tile",         "HOR_hor_period"),
  c("hor_multiplicity", "HOR_n_harmonics"))
for (.pair in .kite_aliases) {
  if (.pair[1] %in% colnames(kite) && !(.pair[2] %in% colnames(kite))) {
    kite[[.pair[2]]] <- kite[[.pair[1]]]
  }
}
kite_key <- sprintf("%s__%s__%d-%d", kite$TRC_ID, kite$seqid,
                    kite$start, kite$end)
kite_idx <- setNames(seq_len(nrow(kite)), kite_key)

message("aggregating per-TRA metrics")
by_q <- split(bl_self, bl_self$qseqid)
res_list <- lapply(names(by_q), function(qid) {
  sub <- by_q[[qid]]
  iv <- merge_intervals(sub$sstart, sub$send)
  slen <- sub$slen[1]
  m <- sub$qlen[1]
  pident_wm <- sum(sub$pident * sub$length) / sum(sub$length)
  strand_plus <- mean(sub$sstart <= sub$send)
  ki <- kite_idx[qid]
  list(
    id              = qid,
    TRC_ID          = strsplit(qid, "__", fixed = TRUE)[[1]][1],
    m               = m,
    array_length    = slen,
    hit_count       = nrow(sub),
    total_align_bp  = sum(sub$length),
    covered_bp      = iv$covered,
    coverage_frac   = iv$covered / slen,
    mean_pident     = pident_wm,
    strand_plus     = strand_plus,
    expected_hits   = slen / m,
    hit_ratio       = nrow(sub) / (slen / m),
    HOR_status      = if (!is.na(ki)) kite$HOR_status[ki] else NA_character_,
    HOR_confidence  = if (!is.na(ki)) kite$HOR_confidence[ki] else NA_real_
  )
})
res_df <- do.call(rbind, lapply(res_list, as.data.frame, stringsAsFactors = FALSE))

all_ids <- names(readDNAStringSet(opt$`consensus-fasta`))
missing_ids <- setdiff(all_ids, names(by_q))
if (length(missing_ids) > 0) {
  message(sprintf("  %d TRAs produced zero self-hits at pident>=%d",
                  length(missing_ids), opt$pident))
  blanks <- data.frame(
    id              = missing_ids,
    TRC_ID          = vapply(missing_ids,
                             function(x) strsplit(x, "__", fixed = TRUE)[[1]][1],
                             character(1)),
    m               = NA_integer_,
    array_length    = NA_integer_,
    hit_count       = 0L,
    total_align_bp  = 0L,
    covered_bp      = 0L,
    coverage_frac   = 0,
    mean_pident     = NA_real_,
    strand_plus     = NA_real_,
    expected_hits   = NA_real_,
    hit_ratio       = 0,
    HOR_status      = NA_character_,
    HOR_confidence  = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(blanks))) {
    ki <- kite_idx[blanks$id[i]]
    if (!is.na(ki)) {
      blanks$m[i]             <- kite$monomer_size[ki]
      blanks$array_length[i]  <- kite$array_length[ki]
      blanks$expected_hits[i] <- kite$array_length[ki] / kite$monomer_size[ki]
      blanks$HOR_status[i]    <- kite$HOR_status[ki]
      blanks$HOR_confidence[i]<- kite$HOR_confidence[ki]
    }
  }
  res_df <- rbind(res_df, blanks)
}

out_tsv <- file.path(opt$`out-dir`, "per_tra_blast_metrics.tsv")
write.table(res_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", out_tsv, sprintf(" (%d TRAs)", nrow(res_df)))

# ---------- 5. Report ------------------------------------------------

message("\n=== headline distributions ===")
qprint <- function(v, label) {
  v <- v[!is.na(v) & is.finite(v)]
  q <- quantile(v, probs = c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95))
  message(sprintf("  %-22s n=%4d p05=%.3f p25=%.3f p50=%.3f p75=%.3f p90=%.3f p95=%.3f",
                  label, length(v), q[1], q[2], q[3], q[4], q[5], q[6]))
  q
}
qcov <- qprint(res_df$coverage_frac,     "coverage_fraction")
qpid <- qprint(res_df$mean_pident / 100, "mean_pident (frac)")
qhr  <- qprint(res_df$hit_ratio,         "hit_ratio (=1 perfect)")

bins <- cut(res_df$coverage_frac,
            breaks = c(-Inf, 0.50, 0.80, 0.90, 0.95, 0.99, Inf),
            labels = c("<50%", "50-80%", "80-90%", "90-95%", "95-99%", ">=99%"),
            right  = FALSE)
cov_tab <- table(bins)
message("\n=== coverage bins (TRA count) ===")
print(cov_tab)

by_hor <- split(res_df, res_df$HOR_status)
if (length(by_hor) > 0) {
  message("\n=== coverage / identity by HOR status ===")
  for (h in names(by_hor)) {
    sub <- by_hor[[h]]
    m1 <- median(sub$coverage_frac, na.rm = TRUE)
    m2 <- median(sub$mean_pident,   na.rm = TRUE)
    message(sprintf("  %-14s n=%4d median_coverage=%.3f median_pident=%.1f%%",
                    h, nrow(sub), m1, m2))
  }
}

ord <- order(res_df$coverage_frac)
worst <- res_df[head(ord, 10), c("id", "m", "array_length", "hit_count",
                                 "coverage_frac", "mean_pident", "HOR_status")]
best  <- res_df[tail(order(res_df$coverage_frac), 10),
                c("id", "m", "array_length", "hit_count",
                  "coverage_frac", "mean_pident", "HOR_status")]

report <- file.path(opt$`out-dir`, "REPORT.md")
f <- file(report, "w")
cat("# BLAST validation — per-TRA consensus vs source array\n\n", file = f)
cat(sprintf("Run parameters: pident >= %d, evalue <= %s, cpu %d\n\n",
            opt$pident, opt$evalue, opt$cpu), file = f)
cat(sprintf("TRAs with at least one self-hit: %d / %d (%.1f%%)\n\n",
            sum(res_df$hit_count > 0), nrow(res_df),
            100 * mean(res_df$hit_count > 0)),
    file = f)
cat("## Coverage distribution\n\n", file = f)
cat("| metric | p05 | p25 | p50 | p75 | p90 | p95 |\n",
    "|---|---|---|---|---|---|---|\n", file = f)
fmt_row <- function(v) {
  v <- v[!is.na(v) & is.finite(v)]
  q <- quantile(v, probs = c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95))
  sprintf("| %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
          q[1], q[2], q[3], q[4], q[5], q[6])
}
cat("| coverage_fraction    ", fmt_row(res_df$coverage_frac), "\n", file = f)
cat("| mean_pident / 100    ", fmt_row(res_df$mean_pident / 100), "\n", file = f)
cat("| hit_ratio (=1 perfect)", fmt_row(res_df$hit_ratio), "\n", file = f)
cat("\n## Coverage bins\n\n", file = f)
for (nm in names(cov_tab)) {
  cat(sprintf("- %s : %d TRAs\n", nm, cov_tab[[nm]]), file = f)
}
if (length(by_hor) > 0) {
  cat("\n## Breakdown by HOR status\n\n", file = f)
  cat("| HOR status | n | median coverage | median pident |\n",
      "|---|---|---|---|\n", file = f)
  for (h in names(by_hor)) {
    sub <- by_hor[[h]]
    cat(sprintf("| %s | %d | %.3f | %.2f%% |\n", h, nrow(sub),
                median(sub$coverage_frac, na.rm = TRUE),
                median(sub$mean_pident, na.rm = TRUE)), file = f)
  }
}
cat("\n## Worst 10 TRAs by coverage\n\n", file = f)
cat("| id | m | array_length | hits | coverage | pident | HOR |\n",
    "|---|---|---|---|---|---|---|\n", file = f)
for (i in seq_len(nrow(worst))) {
  r <- worst[i, ]
  cat(sprintf("| %s | %s | %s | %d | %.3f | %s | %s |\n",
              r$id, as.character(r$m), as.character(r$array_length),
              as.integer(r$hit_count), r$coverage_frac,
              if (is.na(r$mean_pident)) "—" else sprintf("%.1f", r$mean_pident),
              ifelse(is.na(r$HOR_status), "—", r$HOR_status)),
      file = f)
}
cat("\n## Best 10 TRAs by coverage\n\n", file = f)
cat("| id | m | array_length | hits | coverage | pident | HOR |\n",
    "|---|---|---|---|---|---|---|\n", file = f)
for (i in seq_len(nrow(best))) {
  r <- best[i, ]
  cat(sprintf("| %s | %s | %s | %d | %.3f | %s | %s |\n",
              r$id, as.character(r$m), as.character(r$array_length),
              as.integer(r$hit_count), r$coverage_frac,
              if (is.na(r$mean_pident)) "—" else sprintf("%.1f", r$mean_pident),
              ifelse(is.na(r$HOR_status), "—", r$HOR_status)),
      file = f)
}
cat("\n## Files\n\n", file = f)
cat("- `per_tra_blast_metrics.tsv` — per-TRA table of all metrics above.\n", file = f)
cat("- `blast.tsv` — raw blastn output.\n", file = f)
cat("- `all_arrays.fasta` — flat DB of TRA sequences.\n", file = f)
close(f)
message("\nwrote ", report)
