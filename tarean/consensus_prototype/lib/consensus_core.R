# Phase-aligned PFM consensus — core helpers.
# See docs/kite_consensus_implementation_plan.md for the algorithm.
# Source this file via:
#   source(file.path(dirname(sys.frame(1)$ofile), "lib", "consensus_core.R"))
# or from a driver script with explicit path.

suppressPackageStartupMessages({
  library(Biostrings)
})

# Select sampling windows inside a TRA of length L. Returns a list of
# integer 2-vectors c(start, end), 1-based inclusive. If the whole
# array fits within n_windows * window_size the function returns one
# whole-array window.
select_windows <- function(L, window_size = 50000L, n_windows = 16L) {
  if (L <= window_size || L <= window_size * n_windows) {
    return(list(c(1L, L)))
  }
  starts <- as.integer(seq(1L, L - window_size + 1L, length.out = n_windows))
  lapply(starts, function(s) c(s, s + window_size - 1L))
}

# Tally A/C/G/T counts into an m x 4 matrix for one sequence S over
# the given window list at period m. Ns and non-ACGT bases are
# ignored. Positions are 1-based absolute in the array.
tally_bases <- function(S, m, windows) {
  counts <- matrix(0L, nrow = m, ncol = 4L,
                   dimnames = list(NULL, c("A", "C", "G", "T")))
  for (w in windows) {
    w_start <- w[1]; w_end <- w[2]
    sub <- subseq(S, start = w_start, end = w_end)
    bases <- strsplit(as.character(sub), "", fixed = TRUE)[[1]]
    positions <- seq.int(w_start, w_end)
    row_idx <- ((positions - 1L) %% m) + 1L
    for (b in c("A", "C", "G", "T")) {
      hits <- which(bases == b)
      if (length(hits) == 0L) next
      counts[, b] <- counts[, b] + tabulate(row_idx[hits], nbins = m)
    }
  }
  counts
}

# Convert a count matrix into PFM (row-normalised) and consensus
# sequence (argmax base per row). Returns a list.
pfm_and_consensus <- function(counts, pseudo = 0.5) {
  denom <- rowSums(counts) + 4 * pseudo
  pfm <- sweep(counts + pseudo, 1, denom, "/")
  consensus_vec <- colnames(pfm)[max.col(pfm, ties.method = "first")]
  list(pfm = pfm, consensus = paste(consensus_vec, collapse = ""))
}

# Per-column Shannon entropy (base-2) of a PFM. Length m vector.
column_entropy <- function(pfm) {
  apply(pfm, 1, function(p) {
    p <- p[p > 0]
    if (length(p) == 0) return(0)
    -sum(p * log2(p))
  })
}

# Harmonic diagnostic: for each non-DC FFT bin k = 2..k_max (cycles
# per m positions), report the maximum relative amplitude across the
# four base proportion signals, normalised by the total non-DC
# amplitude of that base.  A high value at k = 2 indicates a sub-
# repeat of period m / 2 (for example an HOR with n = 2); at k = 3
# it indicates a sub-repeat of period m / 3.
harmonic_diagnostic <- function(pfm, k_max = 5L) {
  m <- nrow(pfm)
  ks <- 2:k_max
  result <- numeric(length(ks))
  names(result) <- paste0("k", ks)
  if (m < 2 * k_max) return(result)            # too short to measure
  for (k_i in seq_along(ks)) {
    k <- ks[k_i]
    strength <- 0
    for (b in c("A", "C", "G", "T")) {
      sig <- pfm[, b] - mean(pfm[, b])
      ft  <- Mod(fft(sig))
      # FFT index 1 = DC; index 2 = 1 cycle per m; index (k+1) = k cycles per m.
      total <- sum(ft[2:(floor(m / 2) + 1L)])
      if (total > 0) {
        strength <- max(strength, ft[k + 1L] / total)
      }
    }
    result[k_i] <- strength
  }
  result
}

# Full computation for one TRA. Returns all fields the prototype will
# write to disk as one TSV row + the consensus string and PFM matrix.
# `coverage_window_size` and `n_windows` drive windowing.
compute_tra_consensus <- function(dna_string, m,
                                  window_size = 50000L,
                                  n_windows = 16L,
                                  pseudo = 0.5,
                                  k_max = 5L) {
  L <- length(dna_string)
  windows <- select_windows(L, window_size, n_windows)
  counts  <- tally_bases(dna_string, m, windows)
  bases_counted <- sum(counts)
  denom_cov     <- sum(vapply(windows, function(w) w[2] - w[1] + 1L, integer(1)))
  coverage      <- bases_counted / max(1L, denom_cov)
  pc            <- pfm_and_consensus(counts, pseudo = pseudo)
  entropy_vec   <- column_entropy(pc$pfm)
  harm          <- harmonic_diagnostic(pc$pfm, k_max = k_max)
  list(
    consensus       = pc$consensus,
    pfm             = pc$pfm,
    length_bp       = L,
    bases_counted   = bases_counted,
    coverage        = coverage,
    mean_entropy    = mean(entropy_vec),
    entropy_per_col = entropy_vec,
    harmonic        = harm
  )
}

# -------- S4: windowed PFM with per-window rotation alignment ------
# The absolute-position fold in compute_tra_consensus smears the PFM
# when indel noise accumulates across the array (see
# tmp/consensus_prototype_drapa/diagnose_drift/). S4 splits the array
# into windows each short enough that accumulated drift is negligible,
# computes a local count matrix per window, aligns each window to a
# common frame via FFT cross-correlation, and sums the aligned counts.

# Rotate a count matrix by integer offset r (0-based). The rotated
# matrix row i corresponds to row ((i - 1 + r) mod m) + 1 of the
# input (i.e. cyclic left-shift by r).
rotate_counts <- function(counts, r) {
  m <- nrow(counts)
  if (m == 0L) return(counts)
  r <- r %% m
  if (r == 0L) return(counts)
  counts[c((r + 1L):m, seq_len(r)), , drop = FALSE]
}

# Find the integer rotation in [0, m) of `counts_i` that maximises
# column-wise cross-correlation with `counts_ref`. Uses FFT so the
# cost is O(m log m) per base column, well below the naive O(m^2)
# Hamming scan. Returns list(r, score).
align_counts_fft <- function(counts_i, counts_ref) {
  m <- nrow(counts_i)
  stopifnot(nrow(counts_ref) == m, ncol(counts_i) == ncol(counts_ref))
  total <- numeric(m)
  for (b in seq_len(ncol(counts_i))) {
    xi <- counts_i[, b]
    xr <- counts_ref[, b]
    # Mean-centring removes the DC component so the argmax is
    # dominated by true covariation, not raw base composition.
    xi <- xi - mean(xi); xr <- xr - mean(xr)
    ft_i <- fft(xi); ft_r <- fft(xr)
    cc <- Re(fft(Conj(ft_r) * ft_i, inverse = TRUE)) / m
    total <- total + cc
  }
  best_idx <- which.max(total)
  list(r = as.integer(best_idx - 1L), score = total[best_idx])
}

# Tally bases inside a [w_start, w_end] range at period m. Internal
# helper shared by the single- and windowed-PFM code paths.
.tally_range <- function(S, m, w_start, w_end) {
  counts <- matrix(0L, nrow = m, ncol = 4L,
                   dimnames = list(NULL, c("A", "C", "G", "T")))
  sub <- subseq(S, start = w_start, end = w_end)
  bases <- strsplit(as.character(sub), "", fixed = TRUE)[[1]]
  positions <- seq.int(w_start, w_end)
  row_idx <- ((positions - 1L) %% m) + 1L
  for (b in c("A", "C", "G", "T")) {
    hits <- which(bases == b)
    if (length(hits) == 0L) next
    counts[, b] <- counts[, b] + tabulate(row_idx[hits], nbins = m)
  }
  counts
}

# Windowed-PFM consensus. `target_copies_per_window` targets ~N
# monomer copies per window; window size is clamped to
# [min_window_bp, max_window_bp] so tiny m (SSR-like) still gets
# enough sequence per window and huge m does not produce only 1-2
# copies per window.
compute_tra_consensus_windowed <- function(dna_string, m,
                                           target_copies_per_window = 50L,
                                           min_window_bp = 10000L,
                                           max_window_bp = 50000L,
                                           pseudo = 0.5, k_max = 5L) {
  L <- length(dna_string)
  if (m < 2L || L < 2L * m) {
    # Fall back to Option 2 for very short arrays.
    return(compute_tra_consensus(dna_string, m,
                                 window_size = L, n_windows = 1L,
                                 pseudo = pseudo, k_max = k_max))
  }
  W <- max(min_window_bp, min(max_window_bp, target_copies_per_window * m))
  if (W > L) W <- L
  n_windows <- max(1L, L %/% W)
  if (n_windows == 1L) {
    return(compute_tra_consensus(dna_string, m,
                                 window_size = L, n_windows = 1L,
                                 pseudo = pseudo, k_max = k_max))
  }

  # Per-window counts, contiguous non-overlapping cover of the array.
  # Any tail shorter than W is appended to the last full window.
  window_counts <- vector("list", n_windows)
  window_bp     <- integer(n_windows)
  for (w in seq_len(n_windows)) {
    w_start <- (w - 1L) * W + 1L
    w_end <- if (w == n_windows) L else w * W
    window_counts[[w]] <- .tally_range(dna_string, m, w_start, w_end)
    window_bp[w] <- w_end - w_start + 1L
  }

  # Reference = window with the tallest per-row maximum (cleanest
  # periodicity signal). Ties broken by first-in-list. Using the
  # sharpest window as reference improves alignment quality of the
  # noisier ones; using window 1 by default would anchor on whatever
  # drift the array starts with.
  ref_sharpness <- vapply(window_counts, function(c) {
    r <- rowSums(c)
    if (sum(r) == 0L) return(0)
    max(r) / sum(r) * m
  }, numeric(1))
  ref_w <- which.max(ref_sharpness)

  offsets <- integer(n_windows)
  scores  <- numeric(n_windows)
  global  <- window_counts[[ref_w]]
  offsets[ref_w] <- 0L
  scores[ref_w]  <- NA_real_
  for (w in setdiff(seq_len(n_windows), ref_w)) {
    a <- align_counts_fft(window_counts[[w]], window_counts[[ref_w]])
    offsets[w] <- a$r
    scores[w]  <- a$score
    global <- global + rotate_counts(window_counts[[w]], a$r)
  }

  pc <- pfm_and_consensus(global, pseudo = pseudo)
  entropy_vec <- column_entropy(pc$pfm)
  harm <- harmonic_diagnostic(pc$pfm, k_max = k_max)

  bases_counted <- sum(global)
  list(
    consensus       = pc$consensus,
    pfm             = pc$pfm,
    length_bp       = L,
    bases_counted   = bases_counted,
    coverage        = bases_counted / sum(window_bp),
    mean_entropy    = mean(entropy_vec),
    entropy_per_col = entropy_vec,
    harmonic        = harm,
    n_windows       = n_windows,
    window_size_bp  = W,
    ref_window      = ref_w,
    window_offsets  = offsets,
    window_align_scores = scores
  )
}


# Parse a FASTA record name of the form "<seqid>_<start>_<end>" where
# <seqid> may itself contain underscores. Returns list with seqid,
# start, end. If parsing fails, returns all NAs.
parse_record_name <- function(name) {
  parts <- strsplit(name, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  if (n < 3) return(list(seqid = NA_character_, start = NA_integer_, end = NA_integer_))
  end_  <- suppressWarnings(as.integer(parts[n]))
  start_<- suppressWarnings(as.integer(parts[n - 1]))
  if (is.na(end_) || is.na(start_)) {
    return(list(seqid = NA_character_, start = NA_integer_, end = NA_integer_))
  }
  seqid <- paste(parts[seq_len(n - 2)], collapse = "_")
  list(seqid = seqid, start = start_, end = end_)
}
