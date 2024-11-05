#!/usr/bin/env Rscript
library(optparse)

style <-
  "<style>
    .sticky-header {
        position: sticky;
        top: 0;
        background-color: #f1f1f1;
        z-index: 100;
    }
    table {
        background: #FFFFFF;
        border: 1px solid gray;
        border-collapse: collapse;
        color: #fff;
        font: normal 10px verdana, arial, helvetica, sans-serif;
        max-width: 400pt;
        width: 100%; /* Optional: Makes the table width responsive */
    }
    p {
        max-width: 700px;
        text-align: left;
    }
    caption {
        border: 1px solid #5C443A;
        color: #5C443A;
        font-weight: bold;
        font-size: 20pt;
        padding: 16px 14px 18px 10px;
        text-align: center;
    }
    td, th {
        color: #363636;
        padding: 1.0em;
    }
    tr {
        border: 1px dotted gray;
    }
    thead th, tfoot th {
        background: #9C443A;
        color: #FFFFFF;
        padding: 10px;
        text-align: left;
        text-transform: uppercase;
    }
    tbody td a {
        color: #3636FF;
        text-decoration: underline;
    }
    tbody td a:visited {
        color: gray;
        text-decoration: line-through;
    }
    tbody td a:hover {
        text-decoration: underline;
    }
    tbody th a {
        color: #3636FF;
        font-weight: normal;
        text-decoration: none;
    }
    tbody th a:hover {
        color: #363636;
    }
    tbody td + td + td + td a {
        background-image: url('bullet_blue.png');
        background-position: left center;
        background-repeat: no-repeat;
        color: #FFFFFF;
        padding-left: 15px;
    }
    tbody td + td + td + td a:visited {
        background-image: url('bullet_white.png');
        background-position: left center;
        background-repeat: no-repeat;
    }
    tbody th, tbody td {
        text-align: left;
        vertical-align: top;
    }
    tfoot td {
        background: #5C443A;
        color: #FFFFFF;
        padding-top: 3px;
    }
    .odd {
        background: #fff;
    }
    tbody tr:hover {
        background: #EEEEEE;
        border: 1px solid #03476F;
        color: #000000;
    }
</style>

"
# FUNCTIONS:
create_table_data <- function(cls) {
  n <- length(cls)

  # Generate Superfamily names and Anchors
  Superfamily <- paste("TRC Superfamily", seq_len(n))
  Anchor <- paste0("superfamily", seq_len(n))

  # Generate group_names
  group_names <- sapply(cls, function(trcs) {
    trc_numbers <- gsub("TRC_", "", trcs, fixed = TRUE)
    trc_numbers <- as.numeric(trc_numbers)
    trc_numbers <- sort(trc_numbers)
    group_name <- paste0("TRCS_", paste0(trc_numbers, collapse = "_"))
    return(group_name)
  })

  # Generate TRC lists
  TRCs <- sapply(cls, function(trcs) {
    trc_numbers <- gsub("TRC_", "", trcs, fixed = TRUE)
    trc_numbers <- as.numeric(trc_numbers)
    trc_numbers <- sort(trc_numbers)
    trc_list <- paste0("TRC_", trc_numbers, collapse = ", ")
    return(trc_list)
  })

  # Create hyperlinks for the Superfamily column
  SuperfamilyLink <- paste0(
    '<a href="#', Anchor, '">', Superfamily, '</a>'
  )

  # Create the final data frame
  table_data_for_output <- data.frame(
    Superfamily = SuperfamilyLink,
    TRCs = TRCs,
    Anchor = Anchor,
    GroupName = group_names,
    stringsAsFactors = FALSE
  )

  return(table_data_for_output)
}

create_dotplot <- function(sequences) {
  # Load necessary libraries
  # Check if sequences is a DNAStringSet
  if (!inherits(sequences, "DNAStringSet")) {
    stop("Input must be a DNAStringSet")
  }

  # Create temporary files for query and database
  query_fasta <- tempfile(fileext = ".fasta")
  db_fasta <- tempfile(fileext = ".fasta")

  # Write sequences to fasta files
  writeXStringSet(sequences, filepath = query_fasta)
  writeXStringSet(sequences, filepath = db_fasta)

  # Make BLAST database
  system(paste("makeblastdb -in", db_fasta, "-dbtype nucl"))

  # Run BLASTN with tabular output format
  blast_output <- tempfile(fileext = ".txt")
  blast_cmd <- paste("blastn -task blastn -query", query_fasta, "-db", db_fasta,
                     "-outfmt 6 -out", blast_output, "-evalue 10",
                     "-num_alignments 10000000 -num_threads 4 -word_size 4",
                     "-ungapped -perc_identity 60 -dust no"
  )
  system(blast_cmd)

  # Read BLAST output
  blast_results <- read.table(blast_output, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(blast_results) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                               "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  # Determine the orientation of the alignment
  blast_results$orientation <- ifelse(blast_results$sstart <= blast_results$send, "forward", "reverse")

  # Normalize percentage identity for color scaling
  min_pident <- min(blast_results$pident)
  min_pident <- min(60, min_pident)
  max_pident <- max(blast_results$pident)
  pident_range <- max_pident - min_pident
  if (pident_range == 0) {
    blast_results$pident_norm <- rep(0.5, nrow(blast_results))
  } else {
    blast_results$pident_norm <- (blast_results$pident - min_pident) / pident_range
  }

  # Define color ramps for forward and reverse alignments
  forward_col_ramp <- colorRamp(c("lightblue", "blue"))
  reverse_col_ramp <- colorRamp(c("pink", "red"))

  # Function to get color based on pident_norm and orientation
  get_color <- function(pident_norm, orientation) {
    if (orientation == "forward") {
      col_rgb <- forward_col_ramp(pident_norm)
    } else {
      col_rgb <- reverse_col_ramp(pident_norm)
    }
    col_hex <- rgb(col_rgb[1]/255, col_rgb[2]/255, col_rgb[3]/255)
    return(col_hex)
  }

  # Apply the color function to each alignment
  blast_results$color <- mapply(get_color, blast_results$pident_norm, blast_results$orientation)


  # Get sequence lengths and cumulative positions
  seq_lengths <- width(sequences)
  seq_names <- names(sequences)
  cum_lengths <- cumsum(c(0, seq_lengths))
  names(cum_lengths) <- c(seq_names, "end")

  # Map sequence names to cumulative positions
  seq_positions <- data.frame(seqid = seq_names,
                              start = cum_lengths[-length(cum_lengths)],
                              end = cum_lengths[-1])

  # Function to get cumulative positions
  get_cum_pos <- function(seqid, pos) {
    offset <- seq_positions$start[seq_positions$seqid == seqid]
    return(offset + pos)
  }

  # Map BLAST positions to cumulative positions
  blast_results$q_cum_start <- mapply(get_cum_pos, blast_results$qseqid, blast_results$qstart)
  blast_results$q_cum_end <- mapply(get_cum_pos, blast_results$qseqid, blast_results$qend)
  blast_results$s_cum_start <- mapply(get_cum_pos, blast_results$sseqid, blast_results$sstart)
  blast_results$s_cum_end <- mapply(get_cum_pos, blast_results$sseqid, blast_results$send)

  # Create the plot
  plot(NULL, xlim = c(0, max(cum_lengths)), ylim = c(0, max(cum_lengths)),
       xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')


  # Add sequence labels
  N <- length(seq_lengths)
  label_size <- ifelse(N < 5, 1,ifelse(N < 10, 0.8, ifelse(N < 20, 0.6, 0.4)))

  axis(1, at = (cum_lengths[-length(cum_lengths)] + cum_lengths[-1]) / 2, labels = seq_names, las = 2, cex.axis=label_size)
  axis(2, at = (cum_lengths[-length(cum_lengths)] + cum_lengths[-1]) / 2, labels = seq_names, las = 2, cex.axis=label_size)
  axis(3)
  axis(4)

  # what is the density of lines in the plot - length of all segments/(max(cum_lengths)^2)
  density <- sum(blast_results$length)/max(cum_lengths)
  print("Density of lines in the plot")
    print(density)
  print(sequences)
  print("------------------")
  # for density <5 -> LWD=3
  # for density <10 -> LWD=2
  # for density <50 -> LWD=1.5
  # for density >=50 -> LWD=1
  if (density < 5) {
    lwd <- 4
  } else if (density < 10) {
    lwd <- 2
  } else if (density < 50) {
    lwd <- 1.5
  } else if (density < 100){
    lwd <- 1
  } else {
    lwd <- 0.5
  }

  # Plot the alignments
  segments(blast_results$q_cum_start, blast_results$s_cum_start,
           blast_results$q_cum_end, blast_results$s_cum_end,
           col = blast_results$color, lwd = lwd)
  # Add vertical and horizontal lines to separate sequences
  abline(v = cum_lengths, col = "darkgreen", lwd=3)
  abline(h = cum_lengths, col = "darkgreen", lwd=3)


  box()
  # Clean up temporary files
  unlink(c(query_fasta, db_fasta, blast_output))
}

create_blast_db <- function(fasta_file, db_name) {
  cmd <- paste("makeblastdb -in", fasta_file, "-dbtype nucl -parse_seqids -out", db_name)
  if (system(cmd) != 0) {
    stop("Error creating BLAST database")
  }
}

run_blast <- function(query, db, threads, evalue) {
  outfmt <- '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
  blast_out <- tempfile(fileext = ".txt")
  cmd <- paste("blastn -task blastn -db", db, "-query", query,
               "-num_alignments 10000 -num_threads", threads,
               "-out", blast_out, "-outfmt", outfmt, "-evalue", evalue)
  if (system(cmd) != 0) {
    stop("Error running BLAST")
  }
  bl <- read.table(blast_out, header = FALSE, sep = "\t", as.is = TRUE, comment.char = "")
  colnames(bl) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")
  return(bl)
}

make_empty_outputs <- function(prefix){
  # create empty html file
  html_out <- paste(prefix, "_trc_superfamilies.html", sep = "")
  page <- openPage(html_out, title = "TRC Superfamilies", css = style)
  hwrite("No TRC superfamilies found", page = page)
  # create empty file
  cat("", file = paste(prefix, "_superfamilies.csv", sep = ""))
  closePage(page)
}


## get input arguments - FASTA file with TRC dimers
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input FASTA file with TRC dimers"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
                help="Prefix for output files"),
  make_option(c("-t", "--threads"), type="numeric", default=4,
              help="Number of threads"),
  make_option(c("-d", "--debug"), type="logical", default=FALSE,
              help="Save debug information")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# validate input arguments input and prefix must be provided
if (is.null(opt$input) || is.null(opt$prefix)) {
  stop("Input and prefix must be provided")
}

suppressPackageStartupMessages({
  library(igraph)
  library(Biostrings)
  library(R2HTML)
  library(hwriter)
})

dir.create(dirname(opt$prefix), showWarnings = FALSE)
# Load input sequences
fasta <- readDNAStringSet(opt$input)
IDs <- names(fasta)
NewIDs <- paste0(IDs, "_", seq_along(IDs))
fasta_repr <- fasta[!duplicated(IDs)]
TRC_id <- gsub("#.+", "", names(fasta_repr))
names(fasta_repr) <- TRC_id
names(fasta) <- NewIDs

# Create temporary directory and files
workdir <- tempdir()
db <- paste(workdir, "/db", sep="")
blast_out <- paste(workdir, "/blast.out", sep="")
on.exit(unlink(workdir, recursive = TRUE))



# Write sequences to file and run BLAST database
writeXStringSet(fasta, db)
create_blast_db(db, db)
bl <- run_blast(db, db, opt$threads, 1e-5)
bl$ID1 <- gsub("#.+", "", bl$qseqid)
bl$ID2 <- gsub("#.+", "", bl$sseqid)
bl <- bl[bl$ID1 != bl$ID2,]
bl <- bl[order(bl$bitscore, decreasing=TRUE),]
bl2 <- bl[!duplicated((bl[,c("ID1", "ID2")])),]
bl2$score <- (bl2$length * bl2$pident - bl2$gapopen) / ifelse(bl2$qlen > bl2$slen, bl2$qlen, bl2$slen)
score_threshold <-  20
bl3 <- bl2[bl2$score >= score_threshold,]
score_dict <- bl2$score
names(score_dict) <- paste(bl2$ID1, bl2$ID2, sep="|")

# if blast is empty, exit
if (nrow(bl3) == 0) {
  message("No significant similarity found")
  make_empty_outputs(opt$prefix)
  q()
}

# Create graph and identify connected components
g <- graph_from_data_frame(bl3[,c("ID1", "ID2")], directed=FALSE)


gcls <- components(g)
cls <- split(names(gcls$membership), gcls$membership)
# if all components are singletons, exit, there are no TRC superfamilies

if (all(sapply(cls, length) == 1)) {
  message("No TRC superfamilies found")
  make_empty_outputs(opt$prefix)
  q()
}

ord <- order(sapply(cls, length), decreasing=TRUE)
cls <- cls[ord]
fasta_repr_cls <- split(fasta_repr, gcls$membership[TRC_id])[ord]
table_data_for_output <- create_table_data(cls)


#Create directory for images
img_dir <- paste(dirname(opt$prefix), "/dotplots", sep = "")
dir.create(img_dir, showWarnings = FALSE)

# Define HTML output file
html_out <- paste(opt$prefix, "_trc_superfamilies.html", sep = "")

# Open the HTML file for writing
page <- openPage(html_out, title = "TRC Superfamilies", css = style)
hwrite("<h1> TRC Superfamilies</h1>", page = page)
hwrite("<p>  A TRC Superfamily is a group of Tandem Repeat Clusters (TRCs)
 that share significant sequence similarity. The dotplots below show the
 pairwise similarity between TRCs within each superfamily. The color intensity
 of the lines represents the similarity score, with darker colors
  indicating higher similarity. </p>", page = page)
hwrite(table_data_for_output[, c("Superfamily", "TRCs")], page = page, row.names = FALSE, border = 1)


# Loop over groups
for (i in seq_along(cls)) {
  group_name <- gsub("TRC_", "", cls[[i]], fixed = TRUE) %>%
    as.numeric() %>% sort() %>% paste0(collapse = "_") %>% paste0("TRCS_", .)
  N <- sum(nchar(fasta_repr_cls[[i]]))
  size <- log(N) * 200
  min_size <- 500
  size <- max(size, min_size)
  print(size)

  # Generate the image
  img_filename <- paste0(group_name, ".png")
  img_filepath <- file.path(img_dir, img_filename)
  png(img_filepath, width = size, height = size, pointsize = 36)
  par(mar = c(5, 5, 5, 5))
  # some sequenes are too short - if sequence is bellow 70 - multiplay
  L <- nchar(fasta_repr_cls[[i]])
  M <- round(140/L)
  M <- ifelse(M < 1, 1, M)
  fasta_extended <- DNAStringSet(mapply(function(x, y) paste(rep(x, y),collapse = ""), fasta_repr_cls[[i]], M))
  create_dotplot(fasta_extended)
  dev.off()

  # Prepare the TRC list for the section header
  trc_list <- paste("TRC_", strsplit(group_name, "_")[[1]][-1], sep = "", collapse = ", ")

  # Include section header in the HTML
  # Include the section header with anchor
  anchor_name <- table_data_for_output$Anchor[i]
  superfamily_name <- paste("TRC Superfamily", i)
  hwrite(paste0('<h2 id="', anchor_name, '">', superfamily_name, '</h2>'), page = page)
  hwrite("TRCs in superfamily: ", page = page)
  hwrite(paste("<p> ",trc_list,"</p>") , page = page)
  img_src <- file.path(basename(img_dir), img_filename)
  hwriteImage(img_src, page = page, br = TRUE, width=500, link=img_src)
  hwrite('<p><a href="#top">Top</a></p>', page = page)
  hwrite("<hr>", page = page)
}

# Close the HTML file
closePage(page)
# export superfamily assignments to csv table
tbl_out <- data.frame(Superfamily=rep(seq_along(cls), sapply(cls, length), stringsAsFactors = FALSE),
                      TRC=unlist(cls, use.names = FALSE), stringsAsFactors = FALSE)
write.csv(tbl_out, paste(opt$prefix, "_trc_superfamilies.csv", sep = ""), row.names = FALSE)

if (opt$debug) {
save.image(paste(opt$prefix, "_debug.RData", sep = ""))
}

