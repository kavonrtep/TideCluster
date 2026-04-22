#!/usr/bin/env Rscript
library(optparse)

# input parameter:
#   - directory with fasta files
#   - number of cpu
# output
#   - directory with results

# parse command line arguments
option_list <-  list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Directory with fasta files", metavar="character"),
  make_option(c("-c", "--cpu"), type="integer", default=1,
              help="Number of cpu", metavar="integer"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="Output prefix", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
file.arg.name <- "--file="

# get current script path
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub(file.arg.name, "",
                   initial.options[grep(file.arg.name, initial.options)])
script.dir <- normalizePath(dirname(script.name))




# check if all required arguments are provided - input directory and output directory
if (is.null(opt$dir) || is.null(opt$prefix)) {
  stop("Error: Input directory and prefix must be provided.")
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
  library(R2HTML)
  library(hwriter)
})

# creat output directory if it does not exist
output_dir <- paste0(opt$prefix, "_kite")
output_dir_base <- basename(output_dir)
prefix <- opt$prefix
dir.create(output_dir, showWarnings = FALSE)


## HOR CLASSIFICATION CONFIG
# Principles and rationale: see docs/hor_classification.md.
# All tunables live here so they are discoverable; edit here and the
# per-array and per-TRC outputs reflect the new calibration.
HOR_TOL              <- 0.10   # fractional tolerance around integer ratios
HOR_HARMONIC_BONUS   <- 0.5    # weight per extra distinct harmonic k >= 2
HOR_BIN_WEAK         <- 0.10   # confidence thresholds; see §3.4 of docs
HOR_BIN_MODERATE     <- 0.20
HOR_BIN_STRONG       <- 0.40
HOR_MAX_K            <- 5L     # enumerate m_i/k candidates up to k = MAX_K
HOR_GRID_STEP        <- 1L     # bp resolution of the fine candidate grid

## FUNCTIONS
# Score a single candidate base monomer m* against the top-3 peaks
# (ms = c(m1, m2, m3), ss = c(s1, s2, s3)). Returns a list with the
# fit info + confidence, or NULL if m* cannot be interpreted as a
# base-with-harmonic explanation of these peaks.
score_m_star <- function(ms, ss, m_star) {
  if (is.na(m_star) || m_star <= 0) return(NULL)
  total <- sum(ss, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NULL)
  n_peaks <- length(ms)
  ks    <- rep(NA_integer_, n_peaks)
  errs  <- rep(NA_real_, n_peaks)
  close <- rep(0, n_peaks)
  for (i in seq_len(n_peaks)) {
    if (is.na(ms[i])) next
    k <- as.integer(round(ms[i] / m_star))
    if (k < 1L) return(NULL)           # m_star larger than this peak
    ks[i]    <- k
    errs[i]  <- abs(ms[i] - k * m_star) / m_star
    close[i] <- max(0, 1 - errs[i] / HOR_TOL)
  }
  # Require at least one k == 1 (base present) AND at least one k >= 2
  if (!any(ks == 1L, na.rm = TRUE))  return(NULL)
  if (!any(ks >= 2L, na.rm = TRUE)) return(NULL)

  base_w <- sum((ss * close)[!is.na(ks) & ks == 1L], na.rm = TRUE)
  harm_w <- sum((ss * close)[!is.na(ks) & ks >= 2L], na.rm = TRUE)
  f_base <- base_w / total
  f_harm <- harm_w / total
  distinct_h <- length(unique(ks[!is.na(ks) & ks >= 2L & close > 0]))
  bonus   <- 1 + HOR_HARMONIC_BONUS * max(0, distinct_h - 1)
  support <- sqrt(max(0, f_base) * max(0, f_harm))
  conf    <- support * bonus
  contributing_ks <- ks[!is.na(ks) & ks >= 2L & close > 0]
  k_max <- if (length(contributing_ks) > 0) max(contributing_ks) else NA_integer_
  list(m_star = m_star, ks = ks, errs = errs, close = close,
       f_base = f_base, f_harm = f_harm, distinct = distinct_h,
       bonus = bonus, confidence = conf, k_max = k_max, total_s = total)
}

# Build the candidate m* set for an array: observed peaks, their k-th
# fractions for k in 2..HOR_MAX_K, and a 1 bp grid between min(m)/2
# and max(m)/2 (the range where a real base monomer is plausible).
hor_candidates <- function(ms) {
  valid <- ms[!is.na(ms) & ms > 0]
  if (length(valid) == 0) return(integer(0))
  cand <- as.integer(round(valid))
  for (k in 2:HOR_MAX_K) {
    cand <- c(cand, as.integer(round(valid / k)))
  }
  lo <- max(5L, as.integer(min(valid) * 0.5))
  hi <- as.integer(max(valid) * 0.55)
  if (hi > lo) cand <- c(cand, seq.int(lo, hi, by = HOR_GRID_STEP))
  unique(cand[cand > 0])
}

# Find the best-fitting m* for one array (highest confidence).
# Returns the score_m_star() result for that m*, or NULL if no
# candidate satisfied the base + harmonic requirement.
best_hor_fit <- function(m1, m2, m3, s1, s2, s3) {
  ms <- c(m1, m2, m3); ss <- c(s1, s2, s3)
  if (all(is.na(ms))) return(NULL)
  best <- NULL
  for (m_star in hor_candidates(ms)) {
    fit <- score_m_star(ms, ss, m_star)
    if (is.null(fit)) next
    if (is.null(best) || fit$confidence > best$confidence) best <- fit
  }
  best
}

# Map a continuous confidence to a categorical label.
classify_hor_confidence <- function(conf) {
  if (is.na(conf) || conf < HOR_BIN_WEAK)     return("No HOR")
  if (conf < HOR_BIN_MODERATE)                 return("HOR weak")
  if (conf < HOR_BIN_STRONG)                   return("HOR moderate")
  "HOR strong"
}

# Wrap a HOR status string in a coloured HTML badge for use in HTML tables.
hor_status_badge <- function(status) {
  cls <- switch(status,
                "HOR strong"   = "hor-strong",
                "HOR moderate" = "hor-mod",
                "HOR weak"     = "hor-weak",
                "No HOR"       = "hor-none",
                "hor-none")
  sprintf('<span class="hor-badge %s">%s</span>', cls, status)
}

# Wrap a per-TRC count cell with a tinted background matching its category.
# kind is one of "strong", "mod", "weak", "none".
hor_count_cell <- function(value, kind) {
  cls <- switch(kind,
                "strong" = "hor-cnt-strong",
                "mod"    = "hor-cnt-mod",
                "weak"   = "hor-cnt-weak",
                "none"   = "hor-cnt-none",
                "hor-cnt-none")
  sprintf('<span class="%s" style="display:block; padding:2px 8px; text-align:center;">%d</span>',
          cls, as.integer(value))
}

calculate_neighbor_distances <- function(dna_sequence, k) {
  dna_str <- as.character(dna_sequence)
  # Extract all k-mers from the sequence
  start_positions <- 1:(nchar(dna_str) - k + 1)
  kmers <- substring(dna_str, start_positions, start_positions + k - 1)
  # Identify unique k-mers
  # Calculate distances between neighboring occurrences for each unique k-mer
  # and name the list by the k-mer
  kmers_positions <- split(seq_along(kmers), kmers)
  # Calculate distances between neighboring occurrences for each k-mer
  neighbor_distances <- lapply(kmers_positions, function(positions) {
    distances <- diff(positions)
    return(distances)
  })
  # Filter out empty elements and k-mers that occur only once
  neighbor_distances <- neighbor_distances[sapply(neighbor_distances, length) > 0]
  return(neighbor_distances)
}

simplified_findpeaks <- function(x, minpeakdistance = 1, npeaks = 0) {
  stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) == 0)
  xc <- paste(as.character(sign(diff(x))), collapse = "")
  xc <- gsub("1", "+", gsub("-1", "-", xc))
  peakpat <- "[+][-]"
  rc <- gregexpr(peakpat, xc)[[1]]
  if (rc[1] < 0) return(NULL)
  x1 <- rc
  x2 <- rc + attr(rc, "match.length")
  n <- length(x1)
  xv <- xp <- numeric(n)
  for (i in 1:n) {
    xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
    xv[i] <- x[xp[i]]
  }
  X <- cbind(xv, xp, x1, x2)
  sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
  X <- X[sl, , drop = FALSE]
  if (length(X) == 0) return(c())
  if (minpeakdistance > 1) {
    no_peaks <- nrow(X)
    badpeaks <- rep(FALSE, no_peaks)
    for (i in 1:no_peaks) {
      ipos <- X[i, 2]
      if (!badpeaks[i]) {
        dpos <- abs(ipos - X[, 2])
        badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
      }
    }
    X <- X[!badpeaks, , drop = FALSE]
  }
  colnames(X) <- c("peak", "position", "left", "right")
  return(X)
}


random_dna_sequence <- function(L, prob = c(0.25, 0.25, 0.25, 0.25)) {
  # L - length of the sequence
  # return - random DNA sequence of length L
  dna_letters <- c("A", "C", "G", "T")
  dna_sequence <- paste(sample(dna_letters, L, replace = TRUE, prob = prob), collapse = "")
  return(dna_sequence)
}


get_neighgor_distances_background <- function(S, k, N = 10){
  # s - DNA sequence
  # k - length of k-mer
  # N - number of random sequences
  # return - list of neighbor distances for each k-mer
  # get bakgroud distribution on random sequences
  # get probability of each nucleotide from the sequence
  L <- nchar(S)
  actg_letters <- c("A", "C", "T", "G")
  ACTGprob <- table(factor(unlist(strsplit(as.character(S), "")), levels = actg_letters))/L
  pks <- numeric(L)
  for (j in 1:N){
    pks <- pmax(tabulate(unlist(calculate_neighbor_distances(
      random_dna_sequence(L, ACTGprob), k)), nbins = L),
                pks)
  }

  fit <- smooth.spline(x = seq_along(pks), y = pks + 1)
  y_predicted <- predict(fit, x = seq_along(pks))$y
  multiple <- max((pks/y_predicted)[1:round(L/2)], na.rm = TRUE)
  y_predicted <- y_predicted * multiple
  return(y_predicted)
}

get_peaks_from_neighbor_distances <- function(nd,L,  minpeakdistance = 1, threshold = 0.001, background = 0){
  # get backgroud distribution on random sequences

  km <- nchar(names(nd)[1])
  y <- tabulate(unlist(nd),nbins = L)
  # y[y<=background] <- 0
  x <- seq_along(y)
  # get peaks
  peaks <- data.frame(simplified_findpeaks(y, minpeakdistance = minpeakdistance))
  peaks$score <- peaks$peak/(L - peaks$position - km + 1)
  score2 <- peaks$score * log2(peaks$position)
  peaks$score2 <- score2/sum(score2)
  peaks <- peaks[peaks$score2 > threshold,, drop = FALSE]
  # positions of peaks should be larger than half of the sequence length L
  peaks <- peaks[peaks$position < L/2,, drop = FALSE]
  peaks$L <- L
  peaks$background <- background[peaks$position]
  # keep only peaks with peak height larger than background
  if (any(peaks$peak > peaks$background)) {
    peaks <- peaks[peaks$peak > peaks$background,, drop = FALSE]
  } else {
    # get peak closest to background
    ind <- which.max(peaks$peak/peaks$background)
    peaks <- peaks[ind,, drop = FALSE]
  }
  return(peaks)
}

get_peaks_from_seq <- function (S, k = 6, N = 10, minpeakdistance = 1, threshold = 0.001){
  nd <- calculate_neighbor_distances(S, k)
  profile <- tabulate(unlist(nd), nbins = nchar(S))
  nd_background <- get_neighgor_distances_background(S, k, N = 10)
  peaks <- get_peaks_from_neighbor_distances(nd, L = nchar(S), minpeakdistance = minpeakdistance,
                                             threshold = threshold, background = nd_background)
  return(list(peaks=peaks, profile=profile))
}


# END OF FUNCTIONS


# for testing, remove later
if (FALSE) {
  opt$dir <- "/mnt/raid/454_data/TideCluster/analysis/athaliana/tidecluster/tc1_tarean/fasta/"
  opt$output <- "/mnt/raid/454_data/TideCluster/analysis/athaliana/tidecluster/tc1_tarean/peaks/"
  opt$cpu <- 10
}
tcr_files <- dir(opt$dir, pattern = ".fasta", full.names = TRUE)
cpu <- opt$cpu



# Prepare sequences:
s <- sapply(tcr_files, readDNAStringSet)
trc_names <- gsub(".fasta", "", basename(tcr_files))
# the sequence are dimmers, take only the first half using subseq
names(s) <- trc_names
s <- sapply(names(s), function(x) {
  names(s[[x]]) <- paste0(x, ":", names(s[[x]]))
  return(s[[x]])
})

# unlist sequences
s <- sapply(s, function (x) subseq(x, start = 1, end = nchar(x)/2))
s <- Reduce(c, s)

peaks_list <- mclapply(s, get_peaks_from_seq, k = 6, N = 10, minpeakdistance = 1, threshold = 0.001, mc.cores = cpu)

peaks_best_list <- lapply(peaks_list, function(p) {
  x <- p$peaks
  if (nrow(x) > 0) {
    x <- x[which.max(x$score), , drop = FALSE]
  } else {
    x <- data.frame(peak = NA, position = NA, left = NA, right = NA, score = NA, score2 = NA, L = NA)
  }
  return(x)
})
profile_list <- lapply(peaks_list, function(p) {
  x <- p$profile
  return(x)
})
peaks_best <- do.call(rbind, peaks_best_list)
peaks_max_position <- sapply(peaks_list, function(p) {
  x <- p$peaks
  if (nrow(x) > 0) {
    return(max(x$position))
  } else {
    return(NA)
  }
})

# all_peaks in format position1, position2, position3 ; score1, score2, score3
# show best 3 peaks
top3_peaks_list <- lapply(peaks_list, function(p) {
  x <- p$peaks
  if (nrow(x) > 0) {
    x <- x[order(x$score, decreasing = TRUE), , drop = FALSE]
    x <- x[1:3, , drop = FALSE]
    out <- data.frame(l1 = x$position[1],
                      s1 = x$score[1],
                      l2 = x$position[2],
                      s2 = x$score[2],
                      l3 = x$position[3],
                      s3 = x$score[3],
                      L = x$L[1])
  } else {
    out <- data.frame(l1 = NA,
                      s1 = NA,
                      l2 = NA,
                      s2 = NA,
                      l3 = NA,
                      s3 = NA,
                      L = NA)
  }
  return(out)
})
top3_peaks <- do.call(rbind, top3_peaks_list)
peaks_best$ID <- rownames(peaks_best)

peaks_best$TRC_ID <- sapply(strsplit(peaks_best$ID, ":"), "[", 1)
ID_split <- strsplit(peaks_best$ID, "_|:")
peaks_best$seqid <- sapply(ID_split, function(x){N=length(x); paste0(x[3:(N-2)], collapse = "_")})
peaks_best$start <- as.numeric(sapply(ID_split, function(x){N=length(x); x[N-1]}))
peaks_best$end <- as.numeric(sapply(ID_split, function(x){N=length(x); x[N]}))
peaks_best$TRC_index <- as.numeric(sapply(strsplit(peaks_best$TRC_ID, "_"), "[", 2))
# index for export, order by TRC_ID
ord_index <- order(peaks_best$TRC_index)

# make concise table with best peaks
best_peaks_concise <- cbind(
  peaks_best[,c("TRC_ID", "seqid", "start", "end", "position", "score", "L")],
    top3_peaks[,c("l2", "s2", "l3", "s3")])

colnames(best_peaks_concise) <- c("TRC_ID", "seqid", "start", "end", "monomer_size",
                                  "score", "array_length", "monomer_size_2", "score_2",
                                  "monomer_size_3", "score_3")

# Per-array HOR classification via confidence search. Five columns are
# added: HOR_status (4-bin label), HOR_confidence (continuous score),
# HOR_base_monomer (fitted m*, bp), HOR_hor_period (k_max * m*, bp),
# HOR_n_harmonics (distinct k >= 2 that contributed). See
# docs/hor_classification.md for the algorithm and rationale.
hor_fits <- lapply(seq_len(nrow(best_peaks_concise)), function(i) {
  best_hor_fit(best_peaks_concise$monomer_size[i],
               best_peaks_concise$monomer_size_2[i],
               best_peaks_concise$monomer_size_3[i],
               best_peaks_concise$score[i],
               best_peaks_concise$score_2[i],
               best_peaks_concise$score_3[i])
})
best_peaks_concise$HOR_confidence <- vapply(hor_fits,
  function(f) if (is.null(f)) 0 else f$confidence, numeric(1))
best_peaks_concise$HOR_status <- vapply(best_peaks_concise$HOR_confidence,
  classify_hor_confidence, character(1))
best_peaks_concise$HOR_base_monomer <- vapply(hor_fits,
  function(f) if (is.null(f)) NA_integer_ else as.integer(round(f$m_star)),
  integer(1))
best_peaks_concise$HOR_hor_period <- vapply(hor_fits, function(f) {
  if (is.null(f) || is.na(f$k_max)) NA_integer_
  else as.integer(round(f$k_max * f$m_star))
}, integer(1))
best_peaks_concise$HOR_n_harmonics <- vapply(hor_fits,
  function(f) if (is.null(f)) 0L else as.integer(f$distinct), integer(1))

peaks_best$HOR_status       <- best_peaks_concise$HOR_status
peaks_best$HOR_confidence   <- best_peaks_concise$HOR_confidence
peaks_best$HOR_base_monomer <- best_peaks_concise$HOR_base_monomer
peaks_best$HOR_hor_period   <- best_peaks_concise$HOR_hor_period
peaks_best$HOR_n_harmonics  <- best_peaks_concise$HOR_n_harmonics


save.image(file = paste0(output_dir, "/tcr.RData"))
write.table(peaks_best[ord_index,], file = paste0(output_dir, "/monomer_size_best_estimate_stat.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(best_peaks_concise[ord_index,], file = paste0(output_dir, "/monomer_size_top3_estimats.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(profile_list, file = paste0(output_dir, "/peaks_list.RDS"))


xmax <- 10000

dir.create(paste0(output_dir, "/profile_plots/"), showWarnings = FALSE)

for (trc in unique(best_peaks_concise$TRC_ID)) {
  ind <- which(best_peaks_concise$TRC_ID == trc)
  prof_trc_list <- profile_list[ind]
  xmax <- max(peaks_max_position[ind]) * 1.3
  outpng <- paste0(output_dir, "/profile_plots/profile_", trc, ".png")
  prof_trc_matrix <- sapply(prof_trc_list, function(x) {
    x <- x[1:xmax]
    x[is.na(x)] <- 0
    x <- x/max(x)
    return(x)
  })

  # calculate size of the png, it should be proportional to the number of regions + margin for labels
  png(outpng, width = 1000, height = 50 * log2(ncol(prof_trc_matrix) + 1) + 300)
  par(xaxs = "i", yaxs = "i", mar = c(5, 5, 5, 5))

  n_regions <- ncol(prof_trc_matrix)
  plot(NA, xlim = c(1,xmax), ylim = rev(c(1, ncol(prof_trc_matrix)+1)),
       main = paste0(trc, " (Number of arays: ", n_regions, ")"), type='n',
       axes = FALSE, frame.plot = TRUE, xlab="Monomer length [bp] ", ylab = "TRC array index",
       log = "x"
  )

  # convert values to rgb colors, low is white, high is black
  prof_trc_matrix_color <- rgb(0, 0, 0, alpha = prof_trc_matrix/max(prof_trc_matrix))
  ax <- c(1,1.2, 1.5, 2,2.5, 3,5, 4, 7, 8, 10)
  log_ticks <- unique(round(c(ax, ax * 10, ax * 100, ax * 1000, ax * 10000)))

  #points(rep(1:xmax, ncol(prof_trc_matrix)), rep(1:ncol(prof_trc_matrix), each = xmax), col = prof_trc_matrix_color, pch = "|", cex = 2)
  # try rect instead of points
  rect(rep(1:xmax, ncol(prof_trc_matrix)),
       rep(1:ncol(prof_trc_matrix), each = xmax),
       rep(1:xmax, ncol(prof_trc_matrix)),
       rep(1:ncol(prof_trc_matrix), each = xmax) + 1,
       col = prof_trc_matrix_color, border = prof_trc_matrix_color, lwd = 2)


  ax <- pretty(c(1, xmax), n = 20)
  axis(1, las=2, at=log_ticks, labels=log_ticks)
  axis(2, at = 1:ncol(prof_trc_matrix) + 0.5, labels = 1:ncol(prof_trc_matrix))
  dev.off()

  # make plot from top3 peaks, one plot for each TRC_ID

  # store par setting

  pngout <- paste0(output_dir, "/profile_plots/profile_top3_", trc, ".png")
  png(pngout, width = 600, height = 100)
  par(mar = c(3,0,0,0), cex=0.8)
  prf <- rowSums(prof_trc_matrix)/max(rowSums(prof_trc_matrix))
  ylim <- c(0, max(prf))
  xlim <- c(1, xmax)
  plot(NA, xlim = xlim, ylim = ylim,
         type='n',
         axes = FALSE, frame.plot = FALSE, xlab="", ylab = "",
         log = "x"
   )
  points(prf, col = "#00000090", pch = 19, cex = 2, type='h')
  axis(1, las=2, at=log_ticks, labels=log_ticks)
  dev.off()



}

# make html report:
# html_out <- paste0(output_dir, "/kite_report.html")
html_out <- paste0(output_dir, "_report.html")
# first list all unique TRC_ID, sort it by numerical suffix
trc_list <- unique(best_peaks_concise$TRC_ID)
trc_list <- trc_list[order(as.numeric(gsub(".*_", "", trc_list)))]
# start html report
trc_df <- data.frame(TRC_ID = trc_list, stringsAsFactors = FALSE)
# add info about the most frequent monomer size
trc_df$monomer_size <- sapply(trc_df$TRC_ID, function(x) {
  x <- best_peaks_concise[best_peaks_concise$TRC_ID == x, "monomer_size"]
  x <- x[!is.na(x)]
  if (length(x) > 0) {
    x <- names(sort(table(x), decreasing = TRUE))[1]
  } else {
    x <- NA
  }
  return(x)
})

# number of regions
trc_df$number_of_regions <- sapply(trc_df$TRC_ID, function(x) {
  x <- best_peaks_concise[best_peaks_concise$TRC_ID == x, "monomer_size"]
  x <- x[!is.na(x)]
  return(length(x))
})

# Per-TRC HOR rollup: 4-bin counts + median confidence. No single TRC-
# level HOR label (TRCs are often heterogeneous).
hor_count <- function(trc, status) {
  sum(best_peaks_concise$TRC_ID == trc & best_peaks_concise$HOR_status == status)
}
trc_df$N_no_HOR       <- vapply(trc_df$TRC_ID, hor_count, integer(1), "No HOR")
trc_df$N_HOR_weak     <- vapply(trc_df$TRC_ID, hor_count, integer(1), "HOR weak")
trc_df$N_HOR_moderate <- vapply(trc_df$TRC_ID, hor_count, integer(1), "HOR moderate")
trc_df$N_HOR_strong   <- vapply(trc_df$TRC_ID, hor_count, integer(1), "HOR strong")
trc_df$HOR_median_conf <- vapply(trc_df$TRC_ID, function(t) {
  v <- best_peaks_concise$HOR_confidence[best_peaks_concise$TRC_ID == t]
  if (length(v) == 0) NA_real_ else median(v, na.rm = TRUE)
}, numeric(1))

# link image <img> to monomer size profile of top3 peaks for each TRC_ID
trc_df$monomer_size_profile <- sapply(trc_df$TRC_ID, function(x) {
  x <- paste0("<img src=\"",output_dir_base, "/profile_plots/profile_top3_", x, ".png\" width=\"600\">")
  return(x)
})
trc_df$TRC_ID <- paste0("<a href=\"",output_dir_base, "/trc_", trc_df$TRC_ID, ".html\">", trc_df$TRC_ID, "</a>")
# trc_df$TRC_ID <- paste0("<a href=\"#", trc_df$TRC_ID, "\">", trc_df$TRC_ID, "</a>")

source(paste0(script.dir, "/htmlheader.R")) # set htmlheader variable
# add header
cat(htmlheader, file = html_out)
HTML.title("K-Mer Interval Tandem Repeat Estimation (KITE)", file = html_out, HR = 1)
HTML("
The Monomer Size Estimate Score Plot displays weighted scores for
 estimated monomer sizes in tandem repeat clusters (TRC).
  Points represent monomer size estimates, where the score reflects the estimate's
  significance derived from k-mer interval analysis.
  Peaks indicate probable monomer sizes, with higher scores
   denoting stronger evidence.
   For scores for each tandem repeat array within the TRC,
    follow the provided link to see detailed report.
", file = html_out)
HTML(sprintf("
<h3>Higher-order repeat (HOR) classification</h3>
<p>Each tandem repeat array (TRA) is assigned a continuous HOR
confidence score by searching for a base monomer m* that best explains
the three observed peaks as a harmonic series, and measuring how much
real score mass sits on the base and on the harmonics. See
<code>docs/hor_classification.md</code> in the TideCluster source for
the full algorithm. The score is binned into four categories:</p>
<ul style=\"max-width:760px;\">
<li>%s &mdash; confidence &ge; %.2f. Clean base + one or more
well-supported harmonics at integer multiples within %.0f%%.</li>
<li>%s &mdash; %.2f &le; confidence &lt; %.2f. Base and harmonic
both visible but one of them carries modest score mass or sits
slightly outside the tolerance band.</li>
<li>%s &mdash; %.2f &le; confidence &lt; %.2f. Weak evidence: a
harmonic is present but with low score, or the integer-multiple fit
is noisy.</li>
<li>%s &mdash; confidence &lt; %.2f. No supported HOR structure.</li>
</ul>
<p class=\"hor-legend\">Per-TRC counts below show how many of that TRC's
arrays fall into each category, plus the median confidence.</p>",
hor_status_badge("HOR strong"),   HOR_BIN_STRONG, HOR_TOL * 100,
hor_status_badge("HOR moderate"), HOR_BIN_MODERATE, HOR_BIN_STRONG,
hor_status_badge("HOR weak"),     HOR_BIN_WEAK,     HOR_BIN_MODERATE,
hor_status_badge("No HOR"),       HOR_BIN_WEAK), file = html_out)


# add table with TRC_ID and link to subsections within the report
# for each TRC_ID add link to the corresponding section (in HTML.title below)
# Tint the per-TRC HOR count cells before renaming columns for display.
trc_df$N_no_HOR       <- vapply(trc_df$N_no_HOR,       hor_count_cell, character(1), "none")
trc_df$N_HOR_weak     <- vapply(trc_df$N_HOR_weak,     hor_count_cell, character(1), "weak")
trc_df$N_HOR_moderate <- vapply(trc_df$N_HOR_moderate, hor_count_cell, character(1), "mod")
trc_df$N_HOR_strong   <- vapply(trc_df$N_HOR_strong,   hor_count_cell, character(1), "strong")
trc_df$HOR_median_conf <- ifelse(is.na(trc_df$HOR_median_conf), "",
                                 sprintf("%.3f", trc_df$HOR_median_conf))

# adjust column names for html output
colnames(trc_df)[colnames(trc_df) == "monomer_size"] <- "Monomer size <br> primary estimate"
colnames(trc_df)[colnames(trc_df) == "number_of_regions"] <- "Number of arrays"
colnames(trc_df)[colnames(trc_df) == "monomer_size_profile"] <- "Monomer Size Estimate Score Plot"
colnames(trc_df)[colnames(trc_df) == "N_no_HOR"]       <- "No HOR"
colnames(trc_df)[colnames(trc_df) == "N_HOR_weak"]     <- "HOR weak"
colnames(trc_df)[colnames(trc_df) == "N_HOR_moderate"] <- "HOR moderate"
colnames(trc_df)[colnames(trc_df) == "N_HOR_strong"]   <- "HOR strong"
colnames(trc_df)[colnames(trc_df) == "HOR_median_conf"] <- "Median<br>confidence"


HTML(trc_df, header = c("TRC_ID"), rownames = FALSE, align = "c", file = html_out)
HTML("<br>", file = html_out)
# export individual TRC_ID tables:
# get script path, it will run as command line script
# this must be split into parts beacase single page will be too large


for (trc in trc_list){
  html_out_trc <- paste0(output_dir, "/trc_", trc, ".html")
  # start html report
  cat(htmlheader, file = html_out_trc)
  HTML.title(trc, file = html_out_trc, HR = 1)
  # make section for each TRC_ID
  # add link to the main page, (relative)
  # HTML(paste0("<a href=\"../",prefix, "_kite_report.html\">Back to main K-Mer Interval Tandem Repeat Estimation (KITE)</a>"), file = html_out_trc)
  # include image with profile
  HTML("
  This heatmap displays the estimate scores for monomer sizes
  across individual tandem repeat arrays in a TRC, with each
  array represented by a row. Shades from white to grey signify
   score intensity, where darker shades denote higher scores. Table
   below provide three best monomer size estimates for each array.

  ", file = html_out_trc)
  png_rel_path <- paste0("profile_plots/profile_", trc, ".png")
  # use original size of the image
  HTMLInsertGraph(file = html_out_trc, GraphFileName = png_rel_path, Align="left", Width = 1000)
  df1 <- best_peaks_concise[best_peaks_concise$TRC_ID == trc,, drop = FALSE]
  df1$"TRC array index" <- 1:nrow(df1)
  # Render HOR_status as a coloured badge. Base / HOR period / confidence
  # become displayable columns; HOR_n_harmonics is kept in the TSV only.
  df1$HOR_status <- vapply(df1$HOR_status, hor_status_badge, character(1))
  df1$HOR_confidence <- sprintf("%.3f", df1$HOR_confidence)
  df1$HOR_base_monomer <- ifelse(is.na(df1$HOR_base_monomer), "",
                                 as.character(df1$HOR_base_monomer))
  df1$HOR_hor_period <- ifelse(is.na(df1$HOR_hor_period), "",
                               as.character(df1$HOR_hor_period))
  # drop the TSV-only column from the HTML view
  df1$HOR_n_harmonics <- NULL
  colnames(df1)[colnames(df1) == "monomer_size"]   <- "Monomer size<br>(primary estimate)"
  colnames(df1)[colnames(df1) == "monomer_size_2"] <- "Monomer size 2<br>(alternative estimate)"
  colnames(df1)[colnames(df1) == "monomer_size_3"] <- "Monomer size 3<br>(alternative estimate)"
  colnames(df1)[colnames(df1) == "array_length"]   <- "Array length [nt]"
  colnames(df1)[colnames(df1) == "HOR_status"]       <- "HOR status"
  colnames(df1)[colnames(df1) == "HOR_confidence"]   <- "HOR confidence"
  colnames(df1)[colnames(df1) == "HOR_base_monomer"] <- "Base monomer (bp)"
  colnames(df1)[colnames(df1) == "HOR_hor_period"]   <- "HOR period (bp)"
  # columns start and end are numerical, the number should be printent completelly without scientific notation
  df1$start <- format(df1$start, scientific = FALSE)
  df1$end <- format(df1$end, scientific = FALSE)
  HTML(df1, align = "c", rownames = FALSE, file = html_out_trc)
  HTML("<br>", file = html_out_trc)
  HTMLEndFile(file = html_out_trc)
}


HTMLEndFile(file = html_out)


save.image(paste0(output_dir, "/all_data.RData"))


