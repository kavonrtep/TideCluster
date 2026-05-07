# MSA-based consensus core helpers (M_B).
#
# Pipeline for one TRA:
#   1. Find the highest-periodicity k-mer at period m (B2 anchor).
#   2. Extract m-bp monomer copies starting at each anchor occurrence.
#   3. Stratified sampling to a cap (default 150 copies) for MAFFT speed.
#   4. MAFFT alignment.
#   5. Consensus: argmax base per column, skip columns with
#      gap fraction > threshold (C2 "gap-collapse" consensus).
#
# MAFFT is already in the tidecluster conda env. No new dependencies.

suppressPackageStartupMessages({
  library(Biostrings)
})

# ---- anchor k-mer ---------------------------------------------------

# Pick the k-mer whose occurrence positions have the tightest
# distribution modulo m (i.e. the k-mer that fires at the most
# consistent phase within the monomer). Returns NULL when the array
# has no k-mer with at least `min_occurrences` occurrences.
find_anchor_kmer <- function(S, m, k = 8L, min_occurrences = 20L) {
  L <- length(S)
  if (L < k || L < 2L * m) return(NULL)
  dna_str <- as.character(S)
  start_positions <- 1L:(L - k + 1L)
  kmers <- substring(dna_str, start_positions, start_positions + k - 1L)
  kmer_positions <- split(start_positions, kmers)
  kmer_positions <- kmer_positions[lengths(kmer_positions) >= min_occurrences]
  if (length(kmer_positions) == 0L) return(NULL)
  score_anchor <- function(positions) {
    mods <- (positions - 1L) %% m
    tab <- tabulate(mods + 1L, nbins = m)
    max(tab) / length(positions)
  }
  scores <- vapply(kmer_positions, score_anchor, numeric(1))
  top_i <- which.max(scores)
  list(
    kmer      = names(kmer_positions)[top_i],
    positions = sort(kmer_positions[[top_i]]),
    score     = scores[top_i]
  )
}

# ---- monomer extraction + sampling ---------------------------------

# Extract m-bp copies starting at each anchor position. Copies that
# would run past the array end are dropped.
extract_monomers <- function(S, anchor_positions, m) {
  L <- length(S)
  keep <- anchor_positions + m - 1L <= L
  starts <- anchor_positions[keep]
  if (length(starts) == 0L) return(DNAStringSet())
  seqs <- vapply(starts, function(p)
    as.character(subseq(S, start = p, end = p + m - 1L)),
    character(1))
  out <- DNAStringSet(seqs)
  names(out) <- paste0("m", seq_along(out), "_at_", starts)
  out
}

# Stratified down-sample of a DNAStringSet: sort by position, pick
# `n` indices evenly spaced across the sorted list. Keeps the full
# positional span represented.
sample_stratified <- function(ds, n) {
  if (length(ds) <= n) return(ds)
  positions <- as.integer(sub(".*_at_", "", names(ds)))
  ord <- order(positions)
  ds_sorted <- ds[ord]
  idx <- round(seq(1, length(ds_sorted), length.out = n))
  ds_sorted[idx]
}

# ---- MAFFT ----------------------------------------------------------

# Invoke MAFFT via temp files. `mafft_args` is the full CLI tail;
# default is fast single-pass (no iterative refinement).
run_mafft <- function(seqs,
                      mafft_args = "--retree 1 --maxiterate 0 --nuc --thread 1",
                      timeout_sec = 300L) {
  if (length(seqs) < 2L) return(NULL)
  in_fa  <- tempfile(fileext = ".fa")
  out_fa <- tempfile(fileext = ".fa")
  on.exit(unlink(c(in_fa, out_fa)))
  writeXStringSet(seqs, in_fa)
  cmd <- sprintf("timeout %d mafft %s %s > %s 2>/dev/null",
                 as.integer(timeout_sec), mafft_args,
                 shQuote(in_fa), shQuote(out_fa))
  ret <- system(cmd)
  if (ret != 0L || !file.exists(out_fa) || file.size(out_fa) == 0L) {
    return(NULL)
  }
  readDNAStringSet(out_fa)
}

# ---- consensus from MSA --------------------------------------------

# Majority per column, skip columns where the gap fraction exceeds
# `gap_cutoff`. Columns where no single base reaches `min_base_freq`
# are emitted as `N` (still kept — a "consensus" position with low
# confidence). Returns consensus string + diagnostics.
consensus_from_msa <- function(aln, gap_cutoff = 0.5, min_base_freq = 0.4) {
  if (length(aln) == 0L) {
    return(list(consensus = "", n_kept = 0L, n_dropped_gap = 0L,
                mean_entropy = NA_real_, aln_depth = 0L, aln_length = 0L,
                n_N = 0L))
  }
  mat <- do.call(rbind, lapply(seq_along(aln), function(i) {
    strsplit(toupper(as.character(aln[[i]])), "", fixed = TRUE)[[1]]
  }))
  n_rows <- nrow(mat); n_cols <- ncol(mat)
  keep_col <- logical(n_cols)
  out_base <- character(n_cols)
  ent_vec  <- numeric(n_cols)
  n_N      <- 0L
  for (j in seq_len(n_cols)) {
    col <- mat[, j]
    n_gap <- sum(col == "-" | col == "." | col == " " | col == "")
    if (n_gap / n_rows > gap_cutoff) next
    non_gap <- col[col %in% c("A", "C", "G", "T")]
    if (length(non_gap) == 0L) next
    tab <- table(non_gap)
    top_b <- names(tab)[which.max(tab)]
    top_f <- max(tab) / length(non_gap)
    if (top_f < min_base_freq) {
      out_base[j] <- "N"; n_N <- n_N + 1L
    } else {
      out_base[j] <- top_b
    }
    p <- tab / sum(tab)
    ent_vec[j] <- -sum(p * log2(p))
    keep_col[j] <- TRUE
  }
  list(
    consensus     = paste(out_base[keep_col], collapse = ""),
    n_kept        = sum(keep_col),
    n_dropped_gap = n_cols - sum(keep_col),
    mean_entropy  = if (sum(keep_col) > 0) mean(ent_vec[keep_col]) else NA_real_,
    aln_depth     = n_rows,
    aln_length    = n_cols,
    n_N           = n_N
  )
}

# ---- end-to-end one TRA --------------------------------------------

# Returns a list with the consensus + a bag of diagnostic fields that
# the CLI layer flattens into a TSV row.
compute_tra_consensus_msa <- function(S, m,
                                      k = 8L,
                                      max_monomers = 150L,
                                      min_anchor_occurrences = 0L,
                                      gap_cutoff = 0.5,
                                      min_base_freq = 0.4,
                                      mafft_args =
                                        "--retree 1 --maxiterate 0 --nuc --thread 1",
                                      mafft_timeout = 300L) {
  L <- length(S); t0 <- Sys.time()
  # Dynamic anchor threshold: scale with expected copies so short arrays
  # don't fail automatically. `min_anchor_occurrences <= 0` means auto.
  expected_copies <- if (m > 0L) L %/% m else 0L
  if (is.null(min_anchor_occurrences) || is.na(min_anchor_occurrences) ||
      min_anchor_occurrences <= 0L) {
    anchor_min_occ <- max(4L, as.integer(floor(0.5 * expected_copies)))
  } else {
    anchor_min_occ <- as.integer(min_anchor_occurrences)
  }
  base_out <- list(consensus = "", method = "msa_failed",
                   length_bp = L, m = m,
                   anchor_kmer = NA_character_, anchor_score = NA_real_,
                   anchor_count = 0L, anchor_min_occ = anchor_min_occ,
                   expected_copies = expected_copies,
                   n_monomers_extracted = 0L, n_monomers_aligned = 0L,
                   aln_length = NA_integer_, consensus_length = 0L,
                   cons_mean_entropy = NA_real_,
                   n_cols_kept = 0L, n_cols_dropped = 0L, n_N = 0L,
                   wall_sec = NA_real_, reason = NA_character_)

  anchor <- find_anchor_kmer(S, m, k = k, min_occurrences = anchor_min_occ)
  if (is.null(anchor)) {
    out <- base_out; out$reason <- "no_anchor"
    out$wall_sec <- as.numeric(Sys.time() - t0, units = "secs"); return(out)
  }
  mons <- extract_monomers(S, anchor$positions, m)
  if (length(mons) < 4L) {
    out <- base_out
    out$anchor_kmer <- anchor$kmer; out$anchor_score <- anchor$score
    out$anchor_count <- length(anchor$positions)
    out$n_monomers_extracted <- length(mons)
    out$reason <- "too_few_monomers"
    out$wall_sec <- as.numeric(Sys.time() - t0, units = "secs"); return(out)
  }
  sampled <- sample_stratified(mons, max_monomers)
  aln <- run_mafft(sampled, mafft_args = mafft_args, timeout_sec = mafft_timeout)
  if (is.null(aln)) {
    out <- base_out
    out$anchor_kmer <- anchor$kmer; out$anchor_score <- anchor$score
    out$anchor_count <- length(anchor$positions)
    out$n_monomers_extracted <- length(mons)
    out$n_monomers_aligned <- length(sampled)
    out$reason <- "mafft_failed"
    out$wall_sec <- as.numeric(Sys.time() - t0, units = "secs"); return(out)
  }
  cons <- consensus_from_msa(aln,
                             gap_cutoff   = gap_cutoff,
                             min_base_freq = min_base_freq)
  list(
    consensus            = cons$consensus,
    method               = "msa_ok",
    length_bp            = L,
    m                    = m,
    anchor_kmer          = anchor$kmer,
    anchor_score         = anchor$score,
    anchor_count         = length(anchor$positions),
    anchor_min_occ       = anchor_min_occ,
    expected_copies      = expected_copies,
    n_monomers_extracted = length(mons),
    n_monomers_aligned   = length(sampled),
    aln_length           = cons$aln_length,
    consensus_length     = nchar(cons$consensus),
    cons_mean_entropy    = cons$mean_entropy,
    n_cols_kept          = cons$n_kept,
    n_cols_dropped       = cons$n_dropped_gap,
    n_N                  = cons$n_N,
    wall_sec             = as.numeric(Sys.time() - t0, units = "secs"),
    reason               = "ok"
  )
}
