#!/usr/bin/env Rscript
# Single-fragment baseline (Option H pre-step).
#
# For each TRA, take the TideHunter fragment with the highest copy_number
# (= most-supported per-fragment consensus). BLAST the resulting per-TRA
# consensus FASTA against the same flat array DB used by
# validate_drapa_blast.R. Compute coverage_frac per TRA. Compare to the
# M_B (MSA) baseline.
#
# No combining of fragments. This is the "best single fragment" baseline
# that any combining approach has to beat.

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--tidehunter-gff3", type = "character",
              default = "tmp/drapa_run2/drapa_tidehunter.gff3"),
  make_option("--clustering-gff3", type = "character",
              default = "tmp/drapa_run2/drapa_clustering.gff3"),
  make_option("--blast-db",        type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/blast_validation/all_arrays"),
  make_option("--baseline-tsv",    type = "character",
              default = "tmp/consensus_prototype_drapa_msa_v2/blast_validation/per_tra_blast_metrics.tsv"),
  make_option("--out-dir",         type = "character",
              default = "tmp/consensus_prototype_drapa_th_single"),
  make_option("--pident",          type = "integer", default = 85L),
  make_option("--evalue",          type = "character", default = "1e-10"),
  make_option("--cpu",             type = "integer", default = 4L),
  make_option("--force",           type = "logical", default = FALSE,
              action = "store_true")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Map TideHunter fragments to TRAs ---------------------------

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

message("reading TideHunter and clustering GFF3")
th <- read_gff3(opt$`tidehunter-gff3`)
cl <- read_gff3(opt$`clustering-gff3`)
th$cons_seq <- vapply(th$attrs, parse_attr, character(1), k = "consensus_sequence")
th$cons_len <- as.integer(vapply(th$attrs, parse_attr, character(1),
                                 k = "consensus_length"))
th$cn       <- as.numeric(vapply(th$attrs, parse_attr, character(1),
                                 k = "copy_number"))
cl$TRC      <- vapply(cl$attrs, parse_attr, character(1), k = "Name")

th_gr <- GRanges(th$seqid, IRanges(th$start, th$end))
cl_gr <- GRanges(cl$seqid, IRanges(cl$start, cl$end))
ov <- findOverlaps(th_gr, cl_gr, minoverlap = 1L, type = "any")
ovw <- pmin(th$end[queryHits(ov)],   cl$end[subjectHits(ov)]) -
       pmax(th$start[queryHits(ov)], cl$start[subjectHits(ov)]) + 1L
df <- data.frame(th_i = queryHits(ov), cl_i = subjectHits(ov), w = ovw)
df <- df[order(df$th_i, -df$w), ]
df <- df[!duplicated(df$th_i), ]

th$TRC    <- NA_character_
th$tra_id <- NA_character_
th$TRC[df$th_i]    <- cl$TRC[df$cl_i]
th$tra_id[df$th_i] <- sprintf("%s__%s__%d-%d",
                              cl$TRC[df$cl_i], cl$seqid[df$cl_i],
                              cl$start[df$cl_i], cl$end[df$cl_i])
mapped <- th[!is.na(th$tra_id), ]
message(sprintf("  %d TideHunter rows mapped to %d TRAs",
                nrow(mapped), length(unique(mapped$tra_id))))

# ---- 2. Pick best fragment per TRA ---------------------------------
# Best = highest copy_number * consensus_length (= total bp accounted for).
# Tiebreak by consensus_length (longer monomer is usually the true one when
# TideHunter splits a HOR into fractional periods).

mapped$score <- mapped$cn * mapped$cons_len
mapped <- mapped[order(mapped$tra_id, -mapped$score, -mapped$cons_len), ]
best <- mapped[!duplicated(mapped$tra_id), ]
message(sprintf("  picked 1 fragment per TRA: %d fragments", nrow(best)))

cat("\n=== best-fragment stats ===\n")
cat("consensus_length: "); print(summary(best$cons_len))
cat("copy_number:      "); print(summary(best$cn))

# ---- 3. Write FASTA ------------------------------------------------

fa_path <- file.path(opt$`out-dir`, "per_tra_th_consensus.fasta")
con <- file(fa_path, "w")
for (i in seq_len(nrow(best))) {
  cat(">", best$tra_id[i], "\n", best$cons_seq[i], "\n", sep = "", file = con)
}
close(con)
message(sprintf("wrote %s (%d records)", fa_path, nrow(best)))

# Diagnostics TSV
diag <- data.frame(
  id           = best$tra_id,
  TRC_ID       = best$TRC,
  th_start     = best$start,
  th_end       = best$end,
  th_pident    = best$score,  # placeholder; overwritten below
  cons_length  = best$cons_len,
  copy_number  = best$cn,
  stringsAsFactors = FALSE)
diag$th_pident <- as.numeric(best$score)  # th-row score column
diag$th_score  <- best$cn * best$cons_len
write.table(diag, file.path(opt$`out-dir`, "per_tra_th_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Run blastn -------------------------------------------------

blast_tsv <- file.path(opt$`out-dir`, "blast.tsv")
if (opt$force || !file.exists(blast_tsv)) {
  message(sprintf("running blastn (pident >= %d, evalue %s)",
                  opt$pident, opt$evalue))
  outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
  cmd <- sprintf(
    "blastn -task blastn -db %s -query %s -outfmt %s -perc_identity %d -evalue %s -num_threads %d -max_target_seqs 10000 -out %s",
    shQuote(opt$`blast-db`), shQuote(fa_path),
    outfmt, opt$pident, opt$evalue, opt$cpu, shQuote(blast_tsv))
  t0 <- Sys.time()
  if (system(cmd) != 0) stop("blastn failed")
  message(sprintf("  blastn wall time: %.1f s",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---- 5. Parse + per-TRA metrics ------------------------------------

bl <- read.table(blast_tsv, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("qseqid", "sseqid", "pident", "length",
                               "mismatch", "gapopen", "qstart", "qend",
                               "sstart", "send", "evalue", "bitscore",
                               "qlen", "slen"))
message(sprintf("  raw hits: %d", nrow(bl)))

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
uq <- unique(bl$qseqid); us <- unique(bl$sseqid)
qtab <- do.call(rbind, lapply(uq, function(x) as.data.frame(parse_q(x), stringsAsFactors = FALSE)))
rownames(qtab) <- uq
stab <- do.call(rbind, lapply(us, function(x) as.data.frame(parse_s(x), stringsAsFactors = FALSE)))
rownames(stab) <- us
bl$q_seqid <- qtab[bl$qseqid, "seqid"]; bl$q_start <- qtab[bl$qseqid, "start"]; bl$q_end <- qtab[bl$qseqid, "end"]
bl$s_seqid <- stab[bl$sseqid, "seqid"]; bl$s_start <- stab[bl$sseqid, "start"]; bl$s_end <- stab[bl$sseqid, "end"]
bl <- bl[bl$q_seqid == bl$s_seqid & bl$q_start == bl$s_start & bl$q_end == bl$s_end, ]
message(sprintf("  self-hits: %d", nrow(bl)))

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) return(0L)
  lo <- pmin(starts, ends); hi <- pmax(starts, ends)
  o <- order(lo); lo <- lo[o]; hi <- hi[o]
  cur_lo <- lo[1]; cur_hi <- hi[1]; covered <- 0L
  for (i in seq_along(lo)[-1]) {
    if (lo[i] <= cur_hi + 1L) cur_hi <- max(cur_hi, hi[i])
    else { covered <- covered + (cur_hi - cur_lo + 1L)
           cur_lo <- lo[i]; cur_hi <- hi[i] }
  }
  covered + (cur_hi - cur_lo + 1L)
}
by_q <- split(bl, bl$qseqid)
res <- do.call(rbind, lapply(names(by_q), function(qid) {
  sub <- by_q[[qid]]
  data.frame(
    id            = qid,
    array_length  = sub$slen[1],
    cons_length   = sub$qlen[1],
    hit_count     = nrow(sub),
    covered_bp    = merge_intervals(sub$sstart, sub$send),
    coverage_frac = merge_intervals(sub$sstart, sub$send) / sub$slen[1],
    mean_pident   = sum(sub$pident * sub$length) / sum(sub$length),
    stringsAsFactors = FALSE)
}))
write.table(res, file.path(opt$`out-dir`, "per_tra_th_blast_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_total <- nrow(best)
ok99 <- sum(res$coverage_frac >= 0.99)
ok95 <- sum(res$coverage_frac >= 0.95)
ok90 <- sum(res$coverage_frac >= 0.90)
cat(sprintf("\n=== TH single-fragment baseline (n=%d TRAs queried) ===\n", n_total))
cat(sprintf("  >=99%% coverage: %d (%.1f%%)\n", ok99, 100*ok99/n_total))
cat(sprintf("  >=95%% coverage: %d (%.1f%%)\n", ok95, 100*ok95/n_total))
cat(sprintf("  >=90%% coverage: %d (%.1f%%)\n", ok90, 100*ok90/n_total))
cat(sprintf("  no self-hit:    %d\n", n_total - nrow(res)))

# ---- 6. Compare to M_B baseline ------------------------------------

if (file.exists(opt$`baseline-tsv`)) {
  base_df <- read.table(opt$`baseline-tsv`, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "", comment.char = "")
  m <- merge(base_df[, c("id", "coverage_frac", "mean_pident")],
             res[,    c("id", "coverage_frac", "mean_pident")],
             by = "id", suffixes = c("_mb", "_th"))
  cat(sprintf("\n=== head-to-head on %d TRAs in BOTH baselines ===\n", nrow(m)))
  for (cut in c(0.90, 0.95, 0.99)) {
    a <- sum(m$coverage_frac_mb >= cut)
    b <- sum(m$coverage_frac_th >= cut)
    both <- sum(m$coverage_frac_mb >= cut & m$coverage_frac_th >= cut)
    only_mb <- sum(m$coverage_frac_mb >= cut & m$coverage_frac_th <  cut)
    only_th <- sum(m$coverage_frac_mb <  cut & m$coverage_frac_th >= cut)
    cat(sprintf(" cut>=%.2f  M_B:%4d  TH:%4d  both:%4d  M_B-only:%3d  TH-only:%3d\n",
                cut, a, b, both, only_mb, only_th))
  }

  # M_B failures recovered by TH
  mb_fail <- m$coverage_frac_mb < 0.90
  th_pass <- m$coverage_frac_th >= 0.90
  cat(sprintf("\n  M_B fails (<0.90) but TH passes: %d\n",
              sum(mb_fail & th_pass)))
  cat(sprintf("  M_B passes (>=0.90) but TH fails: %d\n",
              sum(!mb_fail & !th_pass)))

  # Where TH wins — examples
  rescued <- m[mb_fail & th_pass, ]
  rescued$delta <- rescued$coverage_frac_th - rescued$coverage_frac_mb
  rescued <- rescued[order(-rescued$delta), ]
  cat("\n  top 10 TH wins (delta = TH - M_B):\n")
  print(head(rescued[, c("id", "coverage_frac_mb", "coverage_frac_th",
                         "mean_pident_mb", "mean_pident_th", "delta")], 10),
        row.names = FALSE)

  # Where TH loses
  losses <- m[!mb_fail & !th_pass, ]
  losses$delta <- losses$coverage_frac_th - losses$coverage_frac_mb
  losses <- losses[order(losses$delta), ]
  cat("\n  top 10 TH losses (delta = TH - M_B):\n")
  print(head(losses[, c("id", "coverage_frac_mb", "coverage_frac_th",
                        "mean_pident_mb", "mean_pident_th", "delta")], 10),
        row.names = FALSE)

  # Save merged for downstream eyeballing
  write.table(m, file.path(opt$`out-dir`, "th_vs_mb.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
