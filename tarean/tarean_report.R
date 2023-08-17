#!/usr/bin/env Rscript
library(optparse, quietly = TRUE)
suppressPackageStartupMessages({
  library(R2HTML, quietly = TRUE)
  library(hwriter, quietly = TRUE)
  library(Biostrings, quietly = TRUE)
})
# Parse command line arguments
option_list <- list(
  make_option(c('-i', '--input_dir'), action = 'store', type = 'character', help
    = 'directory with TAREAN outputs', default = NA),
  make_option(c('-o', '--output'), action = 'store', type = 'character',
              help = 'output base name for html and csv', default = "tarean_report.html"),
  make_option(c('-g', '--gff_file'), action = 'store', type = 'character',
              help = 'gff file with TR annotation', default = NA)
)

read_gff3 <- function(file_path) {
  # Read the file with read.table
  gff3_data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", quote = "", na.strings = ".", col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

  # Process attributes
  attr_list <- lapply(gff3_data$attributes, process_attributes)
  gff3_data$attributes <- list_to_dataframe(attr_list)
  return(gff3_data)
}

list_to_dataframe <- function(input_list) {
  all_columns <- unique(unlist(lapply(input_list, names)))

  df <- do.call(rbind, lapply(input_list, function(x) {
    row_data <- setNames(vector("list", length(all_columns)), all_columns)
    for (name in names(x)) {
      row_data[[name]] <- x[[name]]
    }
    missing_names <- setdiff(all_columns, names(x))
    for (name in missing_names) {
      row_data[[name]] <- NA
    }
    as.data.frame(row_data, stringsAsFactors = FALSE)
  }))

  return(df)
}


# Process attributes
process_attributes <- function(attributes) {
  key_value_pairs <- strsplit(attributes, ";")[[1]]
  attr_list <- sapply(key_value_pairs, function(pair) {
    split_pair <- strsplit(pair, "=")[[1]]
    return(setNames(split_pair[2], split_pair[1]))
  }, simplify = "data.frame", USE.NAMES = FALSE)
  return(attr_list)
}


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
# sort by INDEX
INDEX <- as.numeric(gsub("TRC_", "",
                         gsub("[.]fasta_tarean",
                              "", tarean_dirs)))
tarean_dirs <- tarean_dirs[order(INDEX)]

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

consensus_library <- character()
consensus_cutoff <- 0.25

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

  include <- summary$total_score > consensus_cutoff
  if (sum(include) > 0){
    consensus_top <- summary$consensus[include]
    consensus_top <- gsub("\n", "", gsub("<pre>", "", gsub("</pre>", "", consensus_top)))
    consensus_top_dimer <- paste0(consensus_top, consensus_top)
    names(consensus_top_dimer) <- paste0("cons_", summary$kmer[include], "_",summary$variant[include],"#", TRC[i])
    consensus_library <- c(consensus_library, consensus_top_dimer)
  } else {
    message("No consensus sequences with score < ", consensus_cutoff, " for ", TRC[i])
  }
}
consensus_library <- DNAStringSet(consensus_library)
# rename sequences - keep syntax TRC_INDEX#TRC_INDEX
names(consensus_library) <- sapply(strsplit(names(consensus_library), "#"), function(x) paste0(x[2], "#", x[2]))


gff <- read_gff3(args$gff_file)

parse_ssr_string <- function(input_string) {
  # Split the string by comma and space
  parts <- unlist(strsplit(input_string, ", "))
  # Extract sequences using regular expressions
  sequences <- gsub(" \\(.*", "", parts)
  # Extract percentages using regular expressions
  percentages <- as.numeric(gsub(".*\\(|%25\\)", "", parts))
  include <- percentages > 10
  sequences <- sequences[include]
  sequences <- sapply(sequences, function(x)paste(rep(x, 200), collapse=""))
  return(sequences)
}


# get info about SSRS - if repeat type is SSR - these must be replaced in concentus library
# because tarean consensus is not usually correct for SSRs
ssr_seq_all <- c()

ssr_gff <- gff[gff$attributes$repeat_type == "SSR",]
if (nrow(ssr_gff)>0){
  # get just one line for each TRC
  ssr_gff <- ssr_gff[!duplicated(ssr_gff$attributes$Name),,drop=FALSE]
  for (ss in 1:nrow(ssr_gff)){
    ssr_seq <- parse_ssr_string(ssr_gff$attributes$ssr[ss])
    names(ssr_seq) <- rep(ssr_gff$attributes$Name[ss], length(ssr_seq))
    ssr_seq_all <- c(ssr_seq_all, ssr_seq)
  }
  ssr_seq_all <- DNAStringSet(ssr_seq_all)
  names(ssr_seq_all) <- paste0(names(ssr_seq_all), "#", names(ssr_seq_all))
  # replace sequences in consensus library, but do not add new ones
  include <- names(ssr_seq_all) %in% names(consensus_library)
  consensus_library <- consensus_library[!names(consensus_library) %in% names(ssr_seq_all)]
  consensus_library <- c(consensus_library, ssr_seq_all[include])
}

writeXStringSet(consensus_library,
                file = paste0(args$output, "_consensus_dimer_library.fasta"))


summary_df <- data.frame(INDEX, TRC,monomer_length, kmer, total_score, n_gap50,
                         tarean_dir,graph_link, logo_link)
# sort by INDEX
summary_df <- summary_df[order(summary_df$INDEX),]
summary_df$number_of_arrays <- table(gff$attributes$Name)[summary_df$TRC]
summary_df$min_array_length <- sapply(summary_df$TRC, function(x) {
  min(gff[gff$attributes$Name == x, ]$end - gff[gff$attributes$Name == x, ]$start)
})
summary_df$max_array_length <- sapply(summary_df$TRC, function(x) {
  max(gff[gff$attributes$Name == x, ]$end - gff[gff$attributes$Name == x, ]$start)
})
summary_df$median_array_length <- sapply(summary_df$TRC, function(x) {
  median(gff[gff$attributes$Name == x, ]$end - gff[gff$attributes$Name == x, ]$start)
})
summary_df$type <- sapply(summary_df$TRC, function(x) {
  unique(gff[gff$attributes$Name == x, ]$attributes$repeat_type)
})
# add SSR sequence
summary_df$SSRs <- sapply(summary_df$TRC, function(x) {
  gsub(",", "<br>", unique(gff[gff$attributes$Name == x, ]$attributes$ssr))
})

minmedmax <- paste0("min:    ",summary_df$min_array_length,"<br>",
                                "median: ",   round(summary_df$median_array_length), "<br>",
                                "max:    ", summary_df$max_array_length)
summary_df$size_of_arrays <- minmedmax
summary_df$Consensus <- consensus
summary_df$Annotation <- sapply(summary_df$TRC, function(x) {
  ann <- unique(gff[gff$attributes$Name == x, ]$attributes$annotation)
  if (is.null(ann)){
    ann <- ""
  }
  ann
})

Total_size <- sapply(summary_df$TRC, function(x) {
  sum(gff[gff$attributes$Name == x, ]$end - gff[gff$attributes$Name == x, ]$start)
})

summary_df$Total_size <-  Total_size


summary_df$TRC_with_link <- paste0("<a href='", input_dir_base, "/TRC_", summary_df$INDEX, ".fasta_tarean/report.html'>", summary_df$TRC, "</a>")

# get current script path
script.name <- sub(file.arg.name, "",
                   initial.options[grep(file.arg.name, initial.options)])

script.dir <- normalizePath(dirname(script.name))
source(paste(script.dir, "/", "htmlheader.R", sep = ''))

# export to html

include_cols <- c("TRC_with_link"="TRC",
                  "monomer_length" = "Monomer size",
                  "total_score" = "Score",
                  "Total_size" = "Total size",
                  "SSRs" = "SSRs",
                  "Annotation" = "Annotation",
                  "number_of_arrays" = "Number of arrays",
                  "size_of_arrays" = "Size of arrays",
                  "Consensus" = "Consensus",
                  "graph_link" = "Graph",
                  "logo_link" = "Logo")

summary_df_out <- summary_df[, names(include_cols)]

# rename columns
names(summary_df_out) <- include_cols[names(summary_df_out)]

html_out <- paste0(args$output, ".html")
csv_out <- paste0(args$output, ".tsv")
credits <- readLines(paste0(script.dir,"/../credits.html"))
cat(htmlheader, file = html_out)
cat(credits, file = html_out, append = TRUE)

is_ssrs <-  x <-  sapply(summary_df_out$SSRs, function(x) {
    if (length(x)==0){
      FALSE
    }else{
      if (is.na(x)){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  }
 )

# remove HTML encoding of "%" character in SSRs and Annotation columns - this encoding
# was there to avoid problems in GFF3 output
summary_df_out$SSRs <- gsub("%25", "%", summary_df_out$SSRs)
summary_df_out$Annotation <- gsub("%25", "%", summary_df_out$Annotation)


if ("character(0)" %in% summary_df_out$SSRs){
  summary_df_out$SSRs[summary_df_out$SSRs == "character(0)"] <- ""
}

if (sum(!is_ssrs) > 0){
  summary_df_out_tr <- summary_df_out[!is_ssrs,]
  HTML.title("Tandem Repeats TAREAN Summary", file = html_out)
  HTML(summary_df_out_tr, file = html_out, title = "TAREAN report",
       caption = "TAREAN report", row.names = FALSE, classfirstline="sticky-header",
       align = 'left'
  )
}

if (sum(is_ssrs) > 0){
  summary_df_out_ssrs <- summary_df_out[is_ssrs,]
  # exclude TAREAN logo and graph,..., it is misleading
  summary_df_out_ssrs$Graph <- NULL
  summary_df_out_ssrs$Logo <- NULL
  summary_df_out_ssrs$Score <- NULL
  summary_df_out_ssrs$"Monomer size" <- NULL
  summary_df_out_ssrs$Consensus <- NULL
  HTML.title("Simple Sequence Repeats Summary", file = html_out)
  HTML(summary_df_out_ssrs, file = html_out, title = "TAREAN report",
       caption = "TAREAN report - SSRs", row.names = FALSE, classfirstline="sticky-header",
       align = 'left')

}


HTMLEndFile(file = html_out)

# reformat consensus sequence for tsv output
# the html tags and end of line characters must be removed
consensus_best <- gsub("\n", "",  gsub("<pre>", "", summary_df_out$Consensus))

write.table(apply(summary_df, 2, as.character),
            file = csv_out, sep = "\t",
            row.names = FALSE, quote = FALSE)
