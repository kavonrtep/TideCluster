#!/usr/bin/env Rscript
# Synthetic short-smoke fixtures + a driver that runs the prototype
# consensus_pfm.R in single mode on each one and asserts pass
# criteria from the plan §5.1.
#
# Run from repo root:
#   Rscript tarean/consensus_prototype/tests/make_synthetic.R

suppressPackageStartupMessages({
  library(Biostrings)
})

.script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  fa <- a[grep("--file=", a)]
  if (length(fa) > 0) normalizePath(dirname(sub("--file=", "", fa))) else getwd()
})()
repo_root <- normalizePath(file.path(.script_dir, "..", "..", ".."))
cli       <- file.path(repo_root, "tarean", "consensus_prototype", "consensus_pfm.R")
work_root <- file.path(repo_root, "tmp", "consensus_prototype_smoke")
dir.create(work_root, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
ACGT <- c("A", "C", "G", "T")
random_monomer <- function(m) paste(sample(ACGT, m, replace = TRUE), collapse = "")

mutate_seq <- function(s, rate) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  hit <- runif(length(chars)) < rate
  if (any(hit)) {
    chars[hit] <- vapply(chars[hit], function(c) {
      alt <- setdiff(ACGT, c)
      sample(alt, 1L)
    }, character(1))
  }
  paste(chars, collapse = "")
}

insert_indels <- function(s, rate) {
  if (rate <= 0) return(s)
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  n <- length(chars)
  hit <- runif(n) < rate
  out <- character()
  for (i in seq_along(chars)) {
    if (!hit[i]) { out <- c(out, chars[i]); next }
    act <- sample(c("ins", "del"), 1L)
    if (act == "ins") {
      out <- c(out, chars[i], sample(ACGT, 1L))
    }
    # del: drop this base (just skip appending)
  }
  paste(out, collapse = "")
}

build_tandem <- function(monomer, L_target, mutation = 0, indel = 0) {
  # Concatenate enough copies to exceed L_target, with per-copy
  # mutation + indel, then truncate to L_target.
  m <- nchar(monomer)
  n_copies <- as.integer(ceiling(L_target / m)) + 4L
  copies <- replicate(n_copies, insert_indels(mutate_seq(monomer, mutation), indel))
  s <- paste(copies, collapse = "")
  substr(s, 1L, L_target)
}

write_fixture <- function(name, sequence) {
  fa <- file.path(work_root, paste0(name, ".fasta"))
  hdr <- paste0("synthetic_", name, "_1_", nchar(sequence))
  ds <- DNAStringSet(sequence); names(ds) <- hdr
  writeXStringSet(ds, fa)
  fa
}

hamming <- function(a, b) {
  # strict-length Hamming distance between two strings of equal length
  n <- min(nchar(a), nchar(b))
  sum(strsplit(substr(a, 1, n), "", fixed = TRUE)[[1]] !=
      strsplit(substr(b, 1, n), "", fixed = TRUE)[[1]])
}

# For cyclic consensuses the tool may emit a rotated version of the
# true monomer. Report the minimum Hamming distance over all rotations.
min_rotation_hamming <- function(truth, candidate) {
  m <- nchar(truth)
  if (nchar(candidate) != m) return(m)        # length mismatch, worst case
  best <- m
  for (r in 0:(m - 1)) {
    rot <- paste0(substr(truth, r + 1, m), substr(truth, 1, r))
    d <- hamming(rot, candidate)
    if (d < best) best <- d
  }
  best
}

# Fixture definitions.
fx <- list(
  list(name = "clean-short",   m = 50L,  L = 50000L,
       mutation = 0,    indel = 0,   max_hamming = 1L),
  list(name = "noisy-short",   m = 50L,  L = 50000L,
       mutation = 0.02, indel = 0,   max_hamming = 5L),
  # Option 2 intentionally does not handle phase drift from indels
  # (noted in docs/kite_consensus_review.md §4 Option 2 cons). The
  # pass criterion here is just "consensus is produced and is better
  # than random (half of m)", not "consensus is exact".
  list(name = "indel-short",   m = 50L,  L = 50000L,
       mutation = 0.02, indel = 0.001, max_hamming = 30L),
  list(name = "clean-large",   m = 500L, L = 2000000L,
       mutation = 0,    indel = 0,   max_hamming = 10L),
  list(name = "HOR-2x",        m = 171L, L = 500000L,
       mutation = 0.01, indel = 0,   max_hamming = 18L, hor_2x = TRUE),
  list(name = "n-only",        m = 100L, L = 10000L, n_only = TRUE)
)

results <- list()
for (f in fx) {
  cat("\n=== fixture:", f$name, "===\n")
  if (isTRUE(f$n_only)) {
    seq <- paste(rep("N", f$L), collapse = "")
    truth_monomer <- NA
    run_m <- f$m
  } else if (isTRUE(f$hor_2x)) {
    monomer_short  <- random_monomer(f$m)
    monomer_shift  <- mutate_seq(monomer_short, 0.10)
    hor_monomer    <- paste0(monomer_short, monomer_shift)   # length = 2m
    seq <- build_tandem(hor_monomer, f$L,
                        mutation = f$mutation, indel = f$indel)
    truth_monomer <- hor_monomer
    run_m <- nchar(hor_monomer)
  } else {
    truth_monomer <- random_monomer(f$m)
    seq <- build_tandem(truth_monomer, f$L,
                        mutation = f$mutation, indel = f$indel)
    run_m <- f$m
  }
  fa_path <- write_fixture(f$name, seq)
  out_prefix <- file.path(work_root, f$name, f$name)
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  cmd <- sprintf("Rscript %s --mode single --fasta %s --m %d --out %s",
                 cli, fa_path, run_m, out_prefix)
  ret <- system(cmd)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cons_file <- paste0(out_prefix, "_consensus.fasta")
  cons <- if (file.exists(cons_file))
    as.character(readDNAStringSet(cons_file)[[1]])
  else NA
  hd <- if (!is.na(truth_monomer) && !is.na(cons))
    min_rotation_hamming(truth_monomer, cons)
  else NA
  passed <- if (isTRUE(f$n_only)) {
    # The tool does not refuse explicitly; the consensus will be
    # all zeroes (rowSums counts all zero -> pseudo fills argmax).
    # We just require it runs without error and produces something.
    ret == 0
  } else {
    ret == 0 && !is.na(hd) && hd <= f$max_hamming
  }
  results[[f$name]] <- list(
    name = f$name, m = run_m, truth_len = if (!is.na(truth_monomer)) nchar(truth_monomer) else NA,
    hamming = hd, max = f$max_hamming, wall = dt, passed = passed
  )
  cat(sprintf("  wall=%.1fs  hamming=%s  max=%s  pass=%s\n",
              dt,
              if (is.na(hd)) "NA" else as.character(hd),
              as.character(f$max_hamming),
              as.character(passed)))
}

cat("\n=== Summary ===\n")
for (r in results) {
  cat(sprintf("  %-14s  pass=%-5s  hamming=%s/%s  wall=%.2fs\n",
              r$name, as.character(r$passed),
              if (is.na(r$hamming)) "NA" else as.character(r$hamming),
              as.character(r$max), r$wall))
}
all_pass <- all(vapply(results, function(r) isTRUE(r$passed), logical(1)))
cat(if (all_pass) "\nALL PASSED\n" else "\nFAILURES PRESENT\n")
if (!all_pass) quit(status = 1)
