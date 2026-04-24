#!/usr/bin/env Rscript
# Validation diagnostics for the prototype's Drapa batch output.
# Answers the four questions in docs/kite_consensus_implementation_plan.md §5.2.
#
# Run from repo root:
#   Rscript tarean/consensus_prototype/tests/validate_drapa.R

suppressPackageStartupMessages({
  library(Biostrings)
})

OUT   <- "/home/petr/PycharmProjects/TideCluster/tmp/consensus_prototype_drapa"
DRAPA <- "/home/petr/PycharmProjects/TideCluster/tmp/drapa_run2"

diag <- read.table(file.path(OUT, "per_tra_diagnostics.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
cons <- readDNAStringSet(file.path(OUT, "per_tra_consensus.fasta"))
# FASTA headers carry a trailing "  m=..." in single mode but not
# in batch; strip anything past the first whitespace just in case.
names(cons) <- sub("\\s.*", "", names(cons))

message("TRAs with per-TRA consensus: ", length(cons))
message("diagnostic rows: ", nrow(diag))
stopifnot(length(cons) == nrow(diag))

# -------- 1. per-TRA consensus vs TAREAN consensus --------------------
tarean_tsv <- read.table(file.path(DRAPA, "drapa_tarean_report.tsv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE, quote = "\"", comment.char = "")
# Extract TRC -> TAREAN consensus primary sequence + monomer length
tarean_cons <- character()
tarean_len  <- integer()
for (i in seq_len(nrow(tarean_tsv))) {
  trc <- trimws(tarean_tsv$TRC[i])
  raw <- tarean_tsv$Consensus[i]
  if (is.na(raw) || !nzchar(raw)) next
  seq <- gsub("<[^>]+>", "", raw)
  seq <- gsub("[^ACGT]", "", toupper(seq))
  if (nchar(seq) == 0) next
  tarean_cons[trc] <- seq
  tarean_len[trc]  <- suppressWarnings(as.integer(trimws(tarean_tsv$monomer_length[i])))
}
message("TAREAN consensuses loaded: ", length(tarean_cons))

rc_string <- function(s) {
  as.character(reverseComplement(DNAString(s)))
}

# Minimum Hamming distance across all rotations AND the reverse-
# complement rotation set. Tandem repeat arrays within one TRC can be
# on either strand, so rotation-only comparison underestimates
# similarity for half of them on average.
min_rot_hamm <- function(truth, candidate) {
  m <- nchar(truth)
  if (nchar(candidate) != m) return(NA_integer_)
  t_chars <- strsplit(truth, "", fixed = TRUE)[[1]]
  c_chars <- strsplit(candidate, "", fixed = TRUE)[[1]]
  r_chars <- strsplit(rc_string(candidate), "", fixed = TRUE)[[1]]
  best <- m
  for (r in 0:(m - 1)) {
    rot <- c(t_chars[(r + 1):m], if (r > 0) t_chars[1:r] else character())
    d_f <- sum(rot != c_chars)
    d_r <- sum(rot != r_chars)
    if (d_f < best) best <- d_f
    if (d_r < best) best <- d_r
  }
  best
}

# For every TRA whose KITE m1 equals the TAREAN monomer length of its
# TRC, compare the per-TRA consensus to the TAREAN consensus.
edits <- integer()
edits_info <- list()
for (i in seq_len(nrow(diag))) {
  trc <- diag$TRC_ID[i]
  m   <- as.integer(diag$m[i])
  if (!trc %in% names(tarean_cons)) next
  t_len <- tarean_len[[trc]]
  if (is.na(t_len) || t_len != m) next
  cand <- as.character(cons[[diag$id[i]]])
  truth <- tarean_cons[[trc]]
  d <- min_rot_hamm(truth, cand)
  if (!is.na(d)) {
    edits <- c(edits, d)
    edits_info[[length(edits_info) + 1]] <- list(trc = trc, m = m, d = d)
  }
}
message("\n=== Q1: per-TRA consensus vs TAREAN consensus ===")
message("matched TRAs: ", length(edits))
if (length(edits) > 0) {
  qs <- quantile(edits, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
  message("Hamming distance quantiles:")
  print(qs)
  # as fraction of m per case (approximate, each edit_info has its own m)
  frac <- vapply(edits_info, function(e) e$d / e$m, numeric(1))
  message("fraction-of-m quantiles:")
  print(quantile(frac, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1)))
  message(sprintf("median fraction of m: %.4f  (pass threshold <= 0.05)",
                  median(frac)))
}

# -------- 2. within-family variation on TRC_1 -------------------------
trc1_idx <- which(diag$TRC_ID == "TRC_1")
message("\n=== Q2: within-family variation on TRC_1 ===")
message("TRC_1 TRAs: ", length(trc1_idx))
if (length(trc1_idx) >= 2) {
  # Compute all-vs-all min-rotation Hamming within TRC_1 (same m).
  ids <- diag$id[trc1_idx]
  ms  <- as.integer(diag$m[trc1_idx])
  # Use only the dominant m to keep comparisons meaningful
  m_mode <- as.integer(names(sort(table(ms), decreasing = TRUE))[1])
  message("TRC_1 dominant m: ", m_mode, " (",
          sum(ms == m_mode), " / ", length(ms), " TRAs)")
  keep <- which(ms == m_mode)
  n <- length(keep)
  if (n >= 2) {
    dist_mat <- matrix(NA_integer_, n, n)
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      a <- as.character(cons[[ids[keep[i]]]])
      b <- as.character(cons[[ids[keep[j]]]])
      d <- min_rot_hamm(a, b)
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
    offs <- dist_mat[upper.tri(dist_mat)]
    message("pairwise Hamming quantiles:")
    print(quantile(offs, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1)))
    fracs <- offs / m_mode
    message("as fraction of m:")
    print(quantile(fracs, probs = c(0.5, 0.9)))
  }
}

# -------- 3. HOR-dominant / HOR-strong harmonic diagnostic ----------
message("\n=== Q3: HOR harmonic diagnostic ===")
kite_top3 <- read.table(
  file.path(DRAPA, "drapa_kite", "monomer_size_top3_estimats.csv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
  comment.char = "")
keys <- sprintf("%s__%s__%d-%d",
                kite_top3$TRC_ID, kite_top3$seqid,
                kite_top3$start, kite_top3$end)
hor_status <- setNames(as.character(kite_top3$HOR_status), keys)
diag$hor_status <- hor_status[diag$id]

for (cat in c("HOR strong", "HOR moderate", "HOR weak", "No HOR")) {
  sub <- diag[!is.na(diag$hor_status) & diag$hor_status == cat, , drop = FALSE]
  if (nrow(sub) == 0) { message(cat, ": 0 TRAs"); next }
  k2_mean <- mean(as.numeric(sub$harm_k2))
  k3_mean <- mean(as.numeric(sub$harm_k3))
  k2_frac_high <- mean(as.numeric(sub$harm_k2) > 0.15)
  message(sprintf("%-13s  n=%4d  mean(k2)=%.3f  mean(k3)=%.3f  frac(k2>0.15)=%.3f",
                  cat, nrow(sub), k2_mean, k3_mean, k2_frac_high))
}

# -------- 4. runtime / resource ---------------------------------------
message("\n=== Q4: runtime / memory ===")
summary_log <- readLines(file.path(OUT, "summary.log"))
cat(paste(summary_log, collapse = "\n"), "\n")

message("\ndone")
