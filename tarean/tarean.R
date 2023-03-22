#!/usr/bin/env Rscript
library(optparse, quiet = TRUE)
library(parallel)

## get options from command line
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "",
                   initial.options[grep(file.arg.name, initial.options)])
script.dir <- normalizePath(dirname(script.name))
oridir <- getwd()
## parse arguments
option_list <- list(
  make_option(c('-i', '--input_sequences'), action = 'store', type = 'character', help
    = 'fasta file with input sequences', default = NA),
  make_option(c('-o', '--output_dir'), action = 'store', type = 'character', help =
    'output directory', default = "./kmer_analysis"),
  make_option(c('-m', '--min_kmer_length'), action = 'store', type = 'numeric', help =
    'min kmer length', default = 7),
  make_option(c('-x', '--max_kmer_length'), action = 'store', type = 'numeric', help =
    'min kmer length', default = 27),
  make_option(c('-n', '--cpu'), action = 'store', type = 'numeric', help = 'number of
  cpu to use', default = NULL),
  make_option(c('-s', '--sample_size'), action = 'store', type = 'numeric', help =
    'number of sequences to use for analysis, is set to 0 all sequences are used',
              default = 10000),
  make_option(c('-l', '--no_layout'), action = 'store_true', type = 'logical', help =
    'do not calculate graph layout', default = FALSE),
  make_option(c('-p', '--paired'), action = 'store_true', type = 'logical', help =
    'reads are paired', default = FALSE)
)

description <- paste(strwrap("Tandem Trapea Analyzer"), collapse = "\n")
epilogue <- paste(strwrap(""), collapse = "\n")
parser <- OptionParser(
  option_list = option_list,
  epilogue = epilogue,
  description = description,
)
opt <- parse_args(parser, args = commandArgs(TRUE))
## as Rscript
options(OGDF = paste0(script.dir, "/OGDF/runOGDFlayout"))
CPU <- ifelse(is.null(opt$cpu), detectCores(), opt$cpu)
source(paste(script.dir, "/", "methods.R", sep = ''))
source(paste(script.dir, "/", "logo_methods.R", sep = ''))
source(paste(script.dir, "/", "htmlheader.R", sep = ''))
## set number of CPU to use

## run tarean:
out <- tarean(
  opt$input_sequences,
  opt$output_dir,
  opt$min_kmer_length,
  opt$max_kmer_length,
  CPU,
  opt$sample_size,
  !opt$no_layout,
  paired = opt$paired
)

