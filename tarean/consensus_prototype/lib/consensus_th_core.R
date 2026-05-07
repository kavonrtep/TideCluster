# H3 — TideHunter-fragment consensus core helpers.
#
# For each TRA with N >= 1 TideHunter fragments, combine the per-fragment
# consensuses into a single per-TRA monomer:
#   1. filter fragments (drop low copy_number, drop fractional-period
#      detections that are < half the median length);
#   2. length-bucket the survivors and pick the bucket with the highest
#      total copy_number * length (the dominant period);
#   3. dimerise every fragment in that bucket (F + F) so MAFFT can
#      discover rotation;
#   4. MAFFT-align dimers;
#   5. trim the alignment to the most-supported monomer-length window
#      and emit the column-wise majority consensus;
#   6. fall back to the best single fragment when MAFFT fails, only one
#      fragment survives, or the trimmed consensus is too short.
#
# Reuses run_mafft from consensus_msa_core.R; source that file before
# this one.

suppressPackageStartupMessages({
  library(Biostrings)
})

# ---- 1. Filter ------------------------------------------------------

# Drop fragments below a copy-number floor or shorter than `min_relative_length`
# of the median length (typical TideHunter fractional-period detections).
filter_fragments <- function(lengths, cns,
                             min_cn = 2,
                             min_relative_length = 0.5) {
  if (length(lengths) == 0L) return(logical(0))
  med <- median(lengths)
  cns >= min_cn & lengths >= min_relative_length * med
}

# ---- 2. Length cluster ----------------------------------------------

# Greedy single-link clustering by relative length. Walks fragments
# longest-first; assigns each to the first existing bucket whose median
# is within `tol` (relative), or starts a new bucket. Returns an integer
# vector of bucket ids of length(lengths).
cluster_by_length <- function(lengths, tol = 0.15) {
  n <- length(lengths)
  if (n == 0L) return(integer(0))
  ord <- order(-lengths)
  bucket <- integer(n)
  bucket[ord[1]] <- 1L
  meds <- lengths[ord[1]]
  for (i in seq_along(ord)[-1]) {
    j <- ord[i]; L <- lengths[j]
    diffs <- abs(L - meds) / meds
    pick <- which(diffs <= tol)
    if (length(pick) == 0L) {
      bucket[j] <- length(meds) + 1L
      meds <- c(meds, L)
    } else {
      b <- pick[1]
      bucket[j] <- b
      meds[b] <- median(lengths[bucket == b])
    }
  }
  bucket
}

# Bucket weight = sum(cn * length); pick the heaviest. Returns the
# indices of the fragments in the dominant bucket.
pick_dominant_bucket <- function(buckets, lengths, cns) {
  if (length(buckets) == 0L) return(integer(0))
  w <- tapply(lengths * cns, buckets, sum)
  dom_id <- as.integer(names(w)[which.max(w)])
  which(buckets == dom_id)
}

# ---- 3. Dimerise ----------------------------------------------------

dimerise <- function(seqs) {
  s <- as.character(seqs)
  DNAStringSet(paste0(s, s))
}

# ---- 5. Trim dimer alignment ---------------------------------------

# Slide a `target_length`-wide window across the alignment, pick the
# window with the lowest mean per-column gap fraction. Returns the
# sub-matrix and per-column gap fractions for that window. The window
# is monomer-length wide because the input is dimerised; one full copy
# of the monomer always lies inside any target_length-wide slice.
trim_dimer_alignment <- function(aln,
                                 target_length,
                                 gap_cutoff = 0.5,
                                 length_tol = 0.3) {
  if (length(aln) == 0L) return(NULL)
  mat <- do.call(rbind, lapply(seq_along(aln), function(i)
    strsplit(toupper(as.character(aln[[i]])), "", fixed = TRUE)[[1]]))
  n_rows <- nrow(mat); n_cols <- ncol(mat)
  gap_frac <- vapply(seq_len(n_cols), function(j) {
    col <- mat[, j]
    sum(col == "-" | col == "." | col == " " | col == "") / n_rows
  }, numeric(1))
  w <- as.integer(round(target_length))
  if (w < 1L) w <- 1L
  if (w >= n_cols) {
    return(list(start = 1L, end = n_cols,
                mat   = mat, gap_frac = gap_frac))
  }
  csum <- c(0, cumsum(gap_frac))
  win_gap_sum <- csum[(w + 1):length(csum)] - csum[1:(length(csum) - w)]
  s <- which.min(win_gap_sum)
  e <- s + w - 1L
  list(start    = s,
       end      = e,
       mat      = mat[, s:e, drop = FALSE],
       gap_frac = gap_frac[s:e])
}

# Column-wise majority over a kept block; emits N when the top base
# fraction falls below min_base_freq.
majority_consensus <- function(mat, min_base_freq = 0.4) {
  n_cols <- ncol(mat)
  out <- character(n_cols); n_N <- 0L; ent <- numeric(n_cols)
  for (j in seq_len(n_cols)) {
    col <- mat[, j]
    non_gap <- col[col %in% c("A", "C", "G", "T")]
    if (length(non_gap) == 0L) { out[j] <- "N"; n_N <- n_N + 1L; next }
    tab <- table(non_gap)
    top_b <- names(tab)[which.max(tab)]
    top_f <- max(tab) / length(non_gap)
    if (top_f < min_base_freq) { out[j] <- "N"; n_N <- n_N + 1L }
    else out[j] <- top_b
    p <- tab / sum(tab); ent[j] <- -sum(p * log2(p))
  }
  list(consensus    = paste(out, collapse = ""),
       n_N          = n_N,
       mean_entropy = if (length(ent) > 0) mean(ent) else NA_real_)
}

# ---- 6. End-to-end one TRA -----------------------------------------

# `fragments` is a data.frame with columns: seq (character), length
# (integer), copy_number (numeric). Returns a list with the consensus +
# diagnostic fields (flattened into one TSV row by the CLI layer).
compute_tra_consensus_th <- function(fragments,
                                     min_cn = 2,
                                     min_relative_length = 0.5,
                                     length_tol = 0.15,
                                     gap_cutoff = 0.5,
                                     min_base_freq = 0.4,
                                     mafft_args =
                                       "--retree 1 --maxiterate 0 --nuc --thread 1",
                                     mafft_timeout = 300L) {
  t0 <- Sys.time()
  base_out <- list(consensus = "", method = "th_failed",
                   n_fragments = nrow(fragments),
                   n_fragments_used = 0L,
                   bucket_median_length = NA_integer_,
                   consensus_length = 0L,
                   aln_length = NA_integer_,
                   n_cols_kept = 0L, n_cols_dropped = 0L, n_N = 0L,
                   cons_mean_entropy = NA_real_,
                   wall_sec = NA_real_, reason = NA_character_)
  finalize <- function(out) {
    out$wall_sec <- as.numeric(Sys.time() - t0, units = "secs"); out
  }

  if (is.null(fragments) || nrow(fragments) == 0L) {
    return(finalize(modifyList(base_out, list(reason = "no_fragments"))))
  }

  keep <- filter_fragments(fragments$length, fragments$copy_number,
                           min_cn = min_cn,
                           min_relative_length = min_relative_length)
  if (sum(keep) == 0L) {
    return(finalize(modifyList(base_out, list(reason = "all_filtered"))))
  }
  fr <- fragments[keep, , drop = FALSE]

  buckets <- cluster_by_length(fr$length, tol = length_tol)
  dom_idx <- pick_dominant_bucket(buckets, fr$length, fr$copy_number)
  fr <- fr[dom_idx, , drop = FALSE]
  bucket_med <- as.integer(round(median(fr$length)))

  fallback_single <- function(reason) {
    best <- which.max(fr$length * fr$copy_number)
    finalize(modifyList(base_out, list(
      consensus            = fr$seq[best],
      method               = "th_single_ok",
      reason               = reason,
      n_fragments_used     = nrow(fr),
      bucket_median_length = bucket_med,
      consensus_length     = nchar(fr$seq[best]))))
  }

  if (nrow(fr) == 1L) return(fallback_single("single_fragment_after_filter"))

  dimers <- dimerise(DNAStringSet(fr$seq))
  names(dimers) <- sprintf("f%d_len%d_cn%g",
                           seq_len(nrow(fr)), fr$length, fr$copy_number)
  aln <- run_mafft(dimers, mafft_args = mafft_args,
                   timeout_sec = mafft_timeout)
  if (is.null(aln)) return(fallback_single("mafft_failed"))

  trim <- trim_dimer_alignment(aln,
                               target_length = bucket_med,
                               gap_cutoff = gap_cutoff)
  if (is.null(trim) || ncol(trim$mat) < 0.7 * bucket_med) {
    out <- fallback_single("trim_too_short")
    out$aln_length <- length(aln[[1]])
    return(out)
  }

  cs <- majority_consensus(trim$mat, min_base_freq = min_base_freq)
  finalize(list(
    consensus            = cs$consensus,
    method               = "th_msa_ok",
    n_fragments          = nrow(fragments),
    n_fragments_used     = nrow(fr),
    bucket_median_length = bucket_med,
    consensus_length     = nchar(cs$consensus),
    aln_length           = length(aln[[1]]),
    n_cols_kept          = ncol(trim$mat),
    n_cols_dropped       = length(aln[[1]]) - ncol(trim$mat),
    n_N                  = cs$n_N,
    cons_mean_entropy    = cs$mean_entropy,
    reason               = "ok"))
}
