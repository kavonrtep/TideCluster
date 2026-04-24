#!/usr/bin/env Rscript
# diagnose_drift.R — identify the dominant drift mechanism preventing
# Option 2's PFM consensus from matching its source array.
#
# Stratified sample of TRAs from the Drapa BLAST validation output:
# succ = coverage >= 0.95, fail = coverage == 0 / no-hit. Per-TRA
# diagnostics are computed against the top-periodicity k-mer at the
# KITE m1 period; summaries compare succ vs fail.
#
# Six hypothesised drift causes (see design discussion):
#   1 start-offset  — modal_offset != 0, otherwise signal is clean
#   2 fractional m  — drift_slope != 0 (linear cumulative drift)
#   3 indel noise   — drift_rmse and gap_sd elevated
#   4 HOR blending  — mod_entropy_norm elevated, top_kmer count low
#   5 interruptions — gc_sd elevated, step changes in drift
#   6 rearrangements— visible in GC profile as steps
#
# Run from repo root:
#   Rscript tarean/consensus_prototype/tests/diagnose_drift.R

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
})

option_list <- list(
  make_option("--blast-metrics", type = "character",
              default = "tmp/consensus_prototype_drapa/blast_validation/per_tra_blast_metrics.tsv"),
  make_option("--tarean-dir",    type = "character",
              default = "tmp/drapa_run2/drapa_tarean/fasta"),
  make_option("--out-dir",       type = "character",
              default = "tmp/consensus_prototype_drapa/diagnose_drift"),
  make_option("--n-per-bin",     type = "integer", default = 5L,
              help = "TRAs per (HOR_status x stratum) cell; default 5"),
  make_option("--k",             type = "integer", default = 8L,
              help = "k-mer length for anchor search; default 8"),
  make_option("--seed",          type = "integer", default = 42L)
)
opt <- parse_args(OptionParser(option_list = option_list))

repo_root <- normalizePath(getwd())
resolve_path <- function(p) {
  if (startsWith(p, "/")) p
  else normalizePath(file.path(repo_root, p), mustWork = FALSE)
}
opt$`blast-metrics` <- resolve_path(opt$`blast-metrics`)
opt$`tarean-dir`    <- resolve_path(opt$`tarean-dir`)
opt$`out-dir`       <- resolve_path(opt$`out-dir`)

dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(opt$`out-dir`, "plots")
dir.create(plot_dir, showWarnings = FALSE)

set.seed(opt$seed)

# ---------- sample --------------------------------------------------

metrics <- read.table(opt$`blast-metrics`, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "")
metrics$stratum <- NA_character_
metrics$stratum[metrics$coverage_frac >= 0.95] <- "succ"
metrics$stratum[metrics$hit_count == 0]        <- "fail"
metrics <- metrics[!is.na(metrics$stratum) & !is.na(metrics$m), ]
metrics <- metrics[metrics$m >= 50, ]  # skip tiny SSRs to keep k-mer analysis meaningful

hors <- unique(metrics$HOR_status)
picks <- list()
for (h in hors) {
  for (s in c("succ", "fail")) {
    pool <- metrics[metrics$HOR_status == h & metrics$stratum == s, , drop = FALSE]
    if (nrow(pool) == 0) next
    n <- min(opt$`n-per-bin`, nrow(pool))
    idx <- sample.int(nrow(pool), n)
    picks[[paste0(h, "|", s)]] <- pool[idx, ]
  }
}
sampled <- do.call(rbind, picks)
message(sprintf("sampled %d TRAs across %d strata",
                nrow(sampled), length(picks)))

# ---------- load TRC FASTAs (cache on demand) -----------------------

trc_cache <- new.env()

load_trc <- function(trc_id) {
  if (exists(trc_id, envir = trc_cache)) return(get(trc_id, envir = trc_cache))
  fp <- file.path(opt$`tarean-dir`, paste0(trc_id, ".fasta"))
  if (!file.exists(fp)) return(NULL)
  ds <- readDNAStringSet(fp)
  halves <- lapply(ds, function(x) subseq(x, start = 1L, end = nchar(x) %/% 2L))
  names(halves) <- names(ds)
  assign(trc_id, halves, envir = trc_cache)
  halves
}

# ---------- per-TRA diagnostics -------------------------------------

diagnose_one <- function(S, m, k = 8L) {
  L <- length(S)
  if (L < max(k, m * 2L)) return(NULL)
  dna_str <- as.character(S)

  start_positions <- 1L:(L - k + 1L)
  kmers <- substring(dna_str, start_positions, start_positions + k - 1L)
  # Only k-mers with at least 20 occurrences — noise filter.
  kmer_positions <- split(start_positions, kmers)
  kmer_positions <- kmer_positions[lengths(kmer_positions) >= 20L]
  if (length(kmer_positions) == 0L) return(NULL)

  # Pick anchor: k-mer whose (position mod m) distribution is most
  # peaked. Score = max bin probability.
  score_anchor <- function(positions) {
    mods <- (positions - 1L) %% m
    tab <- tabulate(mods + 1L, nbins = m)
    max(tab) / length(positions)
  }
  scores <- vapply(kmer_positions, score_anchor, numeric(1))
  top_i <- which.max(scores)
  top_positions <- kmer_positions[[top_i]]
  top_kmer <- names(kmer_positions)[top_i]

  mods  <- (top_positions - 1L) %% m
  h_mod <- tabulate(mods + 1L, nbins = m)
  p_mod <- h_mod / sum(h_mod)
  mod_entropy <- -sum(p_mod[p_mod > 0] * log2(p_mod[p_mod > 0]))
  mod_entropy_norm <- mod_entropy / log2(m)     # 0 = peaked, 1 = uniform
  modal_offset <- (which.max(h_mod) - 1L)       # 0-based

  # Gap analysis — keep only gaps within 50% of m
  gaps <- diff(top_positions)
  close_gaps <- gaps[abs(gaps - m) < m * 0.5]
  gap_mean   <- if (length(close_gaps) > 0) mean(close_gaps)     else NA
  gap_sd     <- if (length(close_gaps) > 1) sd(close_gaps)        else NA

  # Cumulative drift vs best linear fit.
  drift <- cumsum(close_gaps - m)
  if (length(drift) >= 10) {
    fit <- lm(drift ~ seq_along(drift))
    drift_slope <- unname(coef(fit)[2])
    drift_rmse  <- sqrt(mean(resid(fit)^2))
  } else {
    drift_slope <- NA_real_
    drift_rmse  <- NA_real_
  }

  # Sliding-window GC as a crude interruption detector.
  win <- 5000L
  n_win <- max(1L, L %/% win)
  gc_profile <- numeric(n_win)
  for (w in seq_len(n_win)) {
    sub <- subseq(S, start = (w - 1L) * win + 1L,
                     end   = min(w * win, L))
    gc_profile[w] <- letterFrequency(sub, "GC", as.prob = TRUE)
  }
  gc_sd <- if (n_win > 1L) sd(gc_profile) else NA

  list(
    L = L, m = m,
    top_kmer = top_kmer,
    anchor_count = length(top_positions),
    anchor_score = scores[top_i],
    modal_offset = modal_offset,
    mod_entropy_norm = mod_entropy_norm,
    gap_mean = gap_mean, gap_sd = gap_sd,
    drift_slope = drift_slope, drift_rmse = drift_rmse,
    n_gaps = length(close_gaps),
    gc_sd = gc_sd,
    # raw data for plotting
    mods = mods,
    gaps = close_gaps,
    drift = drift,
    gc_profile = gc_profile
  )
}

plot_one <- function(png_path, r, d) {
  png(png_path, width = 1200, height = 800)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  # mod-m distribution
  barplot(tabulate(d$mods + 1L, nbins = d$m),
          main = sprintf("%s  %s  cov=%.3f  anchor=%s (n=%d)",
                         r$id, r$stratum, r$coverage_frac,
                         d$top_kmer, d$anchor_count),
          xlab = "position mod m", ylab = "count", border = NA)
  abline(v = d$modal_offset + 0.5, col = "red", lty = 2)
  # gap hist
  if (length(d$gaps) >= 2) {
    hist(d$gaps, breaks = 60, main = sprintf("gap distribution (m=%d)", d$m),
         xlab = "gap size (bp)", col = "gray80", border = NA)
    abline(v = d$m, col = "red")
    if (!is.na(d$gap_mean)) abline(v = d$gap_mean, col = "blue", lty = 2)
  } else {
    plot.new(); title(main = sprintf("gap distribution (m=%d) — no gaps", d$m))
  }
  # cumulative drift
  if (length(d$drift) >= 2) {
    plot(d$drift, type = "l", main = "cumulative drift",
         xlab = "anchor k-mer index", ylab = "drift (bp)")
    abline(h = 0, col = "gray")
    if (!is.na(d$drift_slope)) {
      abline(a = 0, b = d$drift_slope, col = "blue", lty = 2)
    }
  } else {
    plot.new(); title(main = "cumulative drift — no data")
  }
  # sliding GC
  plot(d$gc_profile, type = "l", main = "sliding GC (5 kb windows)",
       xlab = "window", ylab = "GC fraction")
  dev.off()
}

# ---------- run -----------------------------------------------------

all_rows <- list()
for (i in seq_len(nrow(sampled))) {
  r <- sampled[i, ]
  parts <- strsplit(r$id, "__", fixed = TRUE)[[1]]
  trc <- parts[1]; seqid <- parts[2]
  rng <- strsplit(parts[3], "-", fixed = TRUE)[[1]]
  start_ <- as.integer(rng[1]); end_ <- as.integer(rng[2])
  halves <- load_trc(trc)
  if (is.null(halves)) next
  target_key <- paste0(seqid, "_", start_, "_", end_)
  if (!target_key %in% names(halves)) next
  S <- halves[[target_key]]
  d <- diagnose_one(S, as.integer(r$m), k = opt$k)
  if (is.null(d)) next

  all_rows[[length(all_rows) + 1L]] <- list(
    id = r$id, stratum = r$stratum, HOR = r$HOR_status,
    coverage = r$coverage_frac, m = as.integer(r$m), L = d$L,
    anchor_kmer = d$top_kmer,
    anchor_count = d$anchor_count,
    anchor_score = d$anchor_score,
    modal_offset = d$modal_offset,
    mod_entropy_norm = d$mod_entropy_norm,
    gap_mean = d$gap_mean, gap_sd = d$gap_sd,
    drift_slope = d$drift_slope, drift_rmse = d$drift_rmse,
    gc_sd = d$gc_sd, n_gaps = d$n_gaps
  )
  png_path <- file.path(plot_dir,
                        paste0(gsub("[^A-Za-z0-9]", "_", r$id), ".png"))
  plot_one(png_path, r, d)
}

df <- do.call(rbind, lapply(all_rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
out_tsv <- file.path(opt$`out-dir`, "diagnostics.tsv")
write.table(df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", out_tsv, sprintf(" (%d rows)", nrow(df)))

# ---------- summarise -----------------------------------------------

fmt <- function(v) {
  v <- v[!is.na(v) & is.finite(v)]
  if (length(v) == 0) return("NA")
  sprintf("%.4f (n=%d)", median(v), length(v))
}

cat("\n=== succ vs fail medians (all HOR combined) ===\n")
for (col in c("mod_entropy_norm", "modal_offset",
              "gap_mean", "gap_sd",
              "drift_slope", "drift_rmse", "gc_sd",
              "anchor_score", "anchor_count")) {
  s <- df[df$stratum == "succ", col]
  f <- df[df$stratum == "fail", col]
  cat(sprintf("  %-20s  succ=%-22s  fail=%s\n", col, fmt(s), fmt(f)))
}

cat("\n=== per HOR x stratum (mod_entropy_norm, drift_slope, gap_sd) ===\n")
cat(sprintf("  %-14s %-5s  n  mod_ent  drift_slope  gap_sd  modal_offset\n",
            "HOR", "strat"))
for (h in unique(df$HOR)) {
  for (s in c("succ", "fail")) {
    sub <- df[df$HOR == h & df$stratum == s, , drop = FALSE]
    if (nrow(sub) == 0) next
    cat(sprintf("  %-14s %-5s  %d  %7.3f  %11.3f  %6.2f  %4.0f\n",
                h, s, nrow(sub),
                median(sub$mod_entropy_norm, na.rm = TRUE),
                median(sub$drift_slope, na.rm = TRUE),
                median(sub$gap_sd, na.rm = TRUE),
                median(sub$modal_offset, na.rm = TRUE)))
  }
}

# Per-row table (wide view)
cat("\n=== per-TRA table ===\n")
for (i in seq_len(nrow(df))) {
  cat(sprintf("  %-45s %-14s %-4s cov=%.2f  m=%4d  mod_ent=%.3f  drift_slope=%+7.3f  gap_sd=%5.2f  offset=%3d\n",
              substr(df$id[i], 1, 45), df$HOR[i], df$stratum[i],
              df$coverage[i], df$m[i],
              df$mod_entropy_norm[i], df$drift_slope[i],
              df$gap_sd[i], df$modal_offset[i]))
}

cat("\nplots: ", plot_dir, "\n")
