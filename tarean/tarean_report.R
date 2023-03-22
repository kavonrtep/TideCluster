#!/usr/bin/env Rscript
library(optparse, quietly = TRUE)
library(R2HTML, quietly = TRUE)
library(hwriter, quietly = TRUE)
# Parse command line arguments
option_list <- list(
  make_option(c('-i', '--input_dir'), action = 'store', type = 'character', help
    = 'directory with TAREAN outputs', default = NA),
  make_option(c('-o', '--output'), action = 'store', type = 'character', help =
    'output base name for html and csv', default = "tarean_report.html")
)
description <- paste(strwrap("Summarize TAREAN reports"), collapse = "\n")
epilogue <- paste(strwrap(""), collapse = "\n")

parser <- OptionParser(
  option_list = option_list,
  epilogue = epilogue,
  description = description,
)

initial.options <- commandArgs(trailingOnly = FALSE)
args <- parse_args(parser, args = commandArgs(TRUE))
file.arg.name <- "--file="

# tarean dir is in format TCR_INDEX.fasta_tarean
tarean_dirs <- dir(args$input_dir, full.names = FALSE, pattern = "_tarean")
input_dir_base <- basename(args$input_dir)
INDEX <-  numeric()
kmer <-  numeric()
total_score <-  numeric()
monomer_length <- numeric()
TRC <-  character()
tarean_dir <- character()
graph_link <- character()
logo_link <- character()
n_gap50 <- numeric()
consensus <- character()

for (i in seq_along(tarean_dirs)) {
  td_full <- paste0(args$input_dir, "/", tarean_dirs[i])
  INDEX[i] <- as.numeric(gsub("TRC_", "",
                           gsub("[.]fasta_tarean",
                                "", tarean_dirs[i])))
  # read summary
  summary <- read.table(paste0(td_full, "/summary_table.csv"), header = TRUE,
                        as.is = TRUE, sep = ",")
  j <- 1 # take only first row
  # DF columns to report:
  TRC[i] <-  paste0("TRC_", INDEX[i])
  tarean_dir[i] <-  tarean_dirs[i]
  kmer[i] <-  summary[j, "kmer"]
  total_score[i] <-  summary[j, "total_score"]
  monomer_length[i] <- summary[j, "monomer_length"]
  gl <- paste0(input_dir_base, "/", tarean_dirs[i],"/", summary[j, "graph_link"])
  graph_link[i] <- hwriteImage(gl, table = FALSE, width = 200, height = 200)
  ll <- paste0(input_dir_base, "/", tarean_dirs[i],"/", summary[j, "logo_link"])
  logo_link[i] <- hwriteImage(ll, table = FALSE, width = 200, height = 200)
  n_gap50[i] <- summary[j, "n_gap50"]
  consensus[i] <- summary[j, "consensus"]

}


summary_df <- data.frame(INDEX, TRC,monomer_length, kmer, total_score, n_gap50,
                         tarean_dir,graph_link, logo_link)
# sort by INDEX
summary_df <- summary_df[order(summary_df$INDEX),]


write.table(summary_df, file = "summary_table.csv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# get current script path
script.name <- sub(file.arg.name, "",
                   initial.options[grep(file.arg.name, initial.options)])

script.dir <- normalizePath(dirname(script.name))
source(paste(script.dir, "/", "htmlheader.R", sep = ''))

# export to html
html_out <- paste0(args$output, ".html")
cat(htmlheader, file = html_out)
summary_df$Consensus = consensus
HTML(summary_df, file = html_out, title = "TAREAN report",
       caption = "TAREAN report", row.names = FALSE, sortableDF = TRUE)
HTMLEndFile(file = html_out)
