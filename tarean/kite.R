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
prefix <- opt$prefix
dir.create(output_dir, showWarnings = FALSE)


## FUNCTIONS
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

# link image <img> to monomer size profile of top3 peaks for each TRC_ID
trc_df$monomer_size_profile <- sapply(trc_df$TRC_ID, function(x) {
  x <- paste0("<img src=\"",output_dir, "/profile_plots/profile_top3_", x, ".png\" width=\"600\">")
  return(x)
})
trc_df$TRC_ID <- paste0("<a href=\"",output_dir, "/trc_", trc_df$TRC_ID, ".html\">", trc_df$TRC_ID, "</a>")
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


# add table with TRC_ID and link to subsections within the report
# for each TRC_ID add link to the corresponding section (in HTML.title below)
# adjust column names for html output
colnames(trc_df)[colnames(trc_df) == "monomer_size"] <- "Monomer size <br> primary estimate"
colnames(trc_df)[colnames(trc_df) == "number_of_regions"] <- "Number of arrays"
colnames(trc_df)[colnames(trc_df) == "monomer_size_profile"] <- "Monomer Size Estimate Score Plot"


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
  colnames(df1)[colnames(df1) == "monomer_size"] <- "Monomer size<br>(primary estimate)"
  colnames(df1)[colnames(df1) == "monomer_size_2"] <- "Monomer size 2<br>(alternative estimate)"
  colnames(df1)[colnames(df1) == "monomer_size_3"] <- "Monomer size 3<br>(alternative estimate)"
  colnames(df1)[colnames(df1) == "array_length"] <- "Array length [nt]"
  # columns start and end are numerical, the number should be printent completelly without scientific notation
  df1$start <- format(df1$start, scientific = FALSE)
  df1$end <- format(df1$end, scientific = FALSE)
  HTML(df1, align = "c", rownames = FALSE, file = html_out_trc)
  HTML("<br>", file = html_out_trc)
  HTMLEndFile(file = html_out_trc)
}


HTMLEndFile(file = html_out)


save.image(paste0(output_dir, "/all_data.RData"))


