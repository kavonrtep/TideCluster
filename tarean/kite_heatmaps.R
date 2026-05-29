#!/usr/bin/env Rscript
# Render per-TRC profile heatmaps from kitehor's --periodogram output.
#
# Replaces the per-TRC PNG generation that used to live in the deleted
# tarean/kite.R. Reads the FASTA-like bundle produced by
# `kitehor analyze --periodogram <PATH>` (alternating `>case_id|H` and
# `>case_id|bg` records, each followed by a single line of
# space-separated numbers), groups arrays by their TRC prefix
# (`<TRC>:<seqid>_<start>_<end>` headers), and emits two PNGs per TRC:
#   - profile_<TRC>.png       per-array H[d] heatmap (rows = arrays)
#   - profile_top3_<TRC>.png  per-TRC summed H[d] spectrum
#
# No new conda dependencies — uses base R + scales (already in
# conda-deps.txt for the rest of the TAREAN / comparative-analysis
# plots).

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--periodogram", type = "character", default = NULL,
              help = "kitehor --periodogram bundle (FASTA-like)."),
  make_option("--top3-csv",    type = "character", default = NULL,
              help = "monomer_size_top3_estimats.csv (joined kitehor output)."),
  make_option("--out-dir",     type = "character", default = NULL,
              help = "Output directory for profile_<TRC>.png files.")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (f in c("periodogram", "top3-csv", "out-dir")) {
  if (is.null(opt[[f]])) stop("missing required --", f)
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

# ---- read periodogram bundle ----------------------------------------
# Streaming parse: alternating `>case_id|H` and `>case_id|bg` headers,
# each followed by exactly one numeric line. We keep only the H rows
# (heatmap intensity); bg is available if a future enhancement wants
# to subtract the noise envelope before plotting.
read_periodogram <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con))
  prof <- list()
  current_id  <- NULL
  current_kind <- NULL
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (startsWith(line, ">")) {
      # Header: >TRC_x:seqid_start_end|H length=L kmer=k
      hdr <- substr(line, 2L, nchar(line))
      bar <- regexpr("\\|", hdr)
      if (bar < 1L) next
      current_id   <- substr(hdr, 1L, bar - 1L)
      kind_field   <- substr(hdr, bar + 1L, nchar(hdr))
      current_kind <- sub(" .*$", "", kind_field)
      next
    }
    if (is.null(current_id)) next
    if (current_kind == "H") {
      vals <- as.numeric(strsplit(line, "\\s+", fixed = FALSE)[[1]])
      vals <- vals[!is.na(vals)]
      prof[[current_id]] <- vals
    }
    current_id <- NULL; current_kind <- NULL
  }
  prof
}

# ---- read top-3 CSV to know per-TRC array order ---------------------
read_top3 <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE, check.names = FALSE,
                   comment.char = "")
  df$record_id <- sprintf("%s:%s_%s_%s",
                          df$TRC_ID, df$seqid, df$start, df$end)
  df
}

prof <- read_periodogram(opt$periodogram)
top3 <- read_top3(opt$`top3-csv`)

if (nrow(top3) == 0L || length(prof) == 0L) {
  message("nothing to render (top3 rows: ", nrow(top3),
          ", periodogram records: ", length(prof), ")")
  quit(status = 0)
}

trc_list <- unique(top3$TRC_ID)
trc_list <- trc_list[order(as.numeric(gsub(".*_", "", trc_list)))]

log_ticks <- (function() {
  ax <- c(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 7, 8, 10)
  unique(round(c(ax, ax * 10, ax * 100, ax * 1000, ax * 10000)))
})()

# Find per-array max-position to size the x-axis (mirrors kite.R).
top3$max_position <- pmax(
  suppressWarnings(as.numeric(top3$monomer_size)),
  suppressWarnings(as.numeric(top3$monomer_size_2)),
  suppressWarnings(as.numeric(top3$monomer_size_3)),
  na.rm = TRUE)

rendered <- 0L
for (trc in trc_list) {
  rows <- top3[top3$TRC_ID == trc, , drop = FALSE]
  ids  <- rows$record_id
  arrays <- prof[ids]
  arrays <- arrays[!vapply(arrays, is.null, logical(1))]
  if (length(arrays) == 0L) next

  pmax_vals <- rows$max_position[!is.na(rows$max_position) & rows$max_position > 0]
  if (length(pmax_vals) == 0L) next
  xmax <- max(round(max(pmax_vals) * 1.3), 100L)

  # ---- per-TRC heatmap -------------------------------------------
  png_out <- file.path(opt$`out-dir`, paste0("profile_", trc, ".png"))
  matrix_h <- sapply(arrays, function(v) {
    x <- numeric(xmax)
    n <- min(length(v), xmax)
    if (n > 0L) x[seq_len(n)] <- v[seq_len(n)]
    mx <- max(x)
    if (mx > 0) x / mx else x
  })
  if (max(matrix_h) <= 0) next

  png(png_out, width = 1000,
      height = 50 * log2(ncol(matrix_h) + 1) + 300)
  par(xaxs = "i", yaxs = "i", mar = c(5, 5, 5, 5))
  plot(NA, xlim = c(1, xmax), ylim = rev(c(1, ncol(matrix_h) + 1)),
       main = sprintf("%s (Number of arrays: %d)", trc, ncol(matrix_h)),
       type = "n", axes = FALSE, frame.plot = TRUE,
       xlab = "Monomer length [bp]", ylab = "TRC array index", log = "x")
  cols <- rgb(0, 0, 0, alpha = matrix_h / max(matrix_h))
  rect(rep(seq_len(xmax), ncol(matrix_h)),
       rep(seq_len(ncol(matrix_h)), each = xmax),
       rep(seq_len(xmax), ncol(matrix_h)),
       rep(seq_len(ncol(matrix_h)), each = xmax) + 1,
       col = cols, border = cols, lwd = 2)
  ticks <- log_ticks[log_ticks >= 1 & log_ticks <= xmax]
  axis(1, las = 2, at = ticks, labels = ticks)
  axis(2, at = seq_len(ncol(matrix_h)) + 0.5,
       labels = seq_len(ncol(matrix_h)))
  dev.off()

  # ---- summed-spectrum overlay ------------------------------------
  png_top3 <- file.path(opt$`out-dir`, paste0("profile_top3_", trc, ".png"))
  png(png_top3, width = 600, height = 100)
  par(mar = c(3, 0, 0, 0), cex = 0.8)
  summed <- rowSums(matrix_h)
  if (max(summed) > 0) summed <- summed / max(summed)
  plot(NA, xlim = c(1, xmax), ylim = c(0, max(summed)),
       type = "n", axes = FALSE, frame.plot = FALSE,
       xlab = "", ylab = "", log = "x")
  points(summed, col = "#00000090", pch = 19, cex = 2, type = "h")
  axis(1, las = 2, at = ticks, labels = ticks)
  dev.off()

  rendered <- rendered + 1L
}

message("Rendered heatmaps for ", rendered, " TRC(s) into ", opt$`out-dir`)
