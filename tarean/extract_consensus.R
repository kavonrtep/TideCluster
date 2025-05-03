#!/usr/bin/env Rscript
library(optparse)
# input arguments - table with genomic regions
#                 - fasta file with sequences

option_list <- list(
  make_option(c("-t", "--table"),
    type = "character",
    default = NULL, help = "Table with genomic regions"
  ),
  make_option(c("-f", "--fasta"),
    type = "character",
    default = NULL, help = "Fasta file with sequences"
  ),
  make_option(c("-o", "--output"),
    type = "character",
    default = NULL, help = "Output file with sequences"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check if input and output are provided
if (is.null(opt$table) | is.null(opt$fasta) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All input arguments are mandatory")
}

# Load the library
suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome)
})

# for testings
# opt <- list(table="/mnt/ceph/454_data/Vicia_faba_pangenome/assemblies/Peru14/ver_240815_split_chr1/analysis/repeats/Peru14_v1_split_chr1_241016/TideCluster/default/TideCluster_kite/monomer_size_top3_estimats.csv", fasta="/mnt/ceph/454_data/Vicia_faba_pangenome/assemblies/Peru14/ver_240815_split_chr1/Peru14_pseudomolecules_v1+unanchored_contigs_split_chr1.fasta")
df <- read.table(opt$table, header = TRUE, sep = "\t")
df$max_monomer_size <- apply(
  df[, c("monomer_size", "monomer_size_2", "monomer_size_3")],
  1,
  function(x) max(x, na.rm = TRUE)
)
# load fasta file
sequence <- readDNAStringSet(opt$fasta)

df_sub <- df[df$seqid %in% names(sequence), ]


# Verify that ranges and sequence have matching seqlevels






# Create GRanges object from df
ranges <- makeGRangesFromDataFrame(df_sub, keep.extra.columns = FALSE)
if (!all(seqlevels(ranges) %in% names(sequence))) {
  stop(
    "Sequence levels in ranges do not match sequence names in fasta file. ",
    "Missing sequences: ",
    paste(setdiff(seqlevels(ranges), names(sequence)), collapse = ", ")
  )
}

# Extract all sequences at once using GRanges
extracted_sequences <- getSeq(sequence, ranges)


# Write sequences to output file
writeXStringSet(extracted_sequences, filepath = opt$output)




# ---------------------------
# Set Parameters
# ---------------------------

# k-mer length will be set based on the max monomer size
# the frequency of this k-mer whould be approximatelly
# proportional to array_length/max_monomer_size

n_monomers <- df_sub$array_length / df_sub$max_monomer_size
k <- ceiling(log2(df_sub$max_monomer_size) + 5)
# k cannot be NA, verify
if (any(is.na(k))) {
  stop("k cannot be NA")
}



L <- 330 # Length of sequences to extract (adjust as needed)
L <- 330 # Length of sequences to extract (adjust as needed)

# ---------------------------
# Load the DNA Sequence
# ---------------------------


# Alternatively, load from a FASTA file
# sequence <- readDNAStringSet("path_to_your_sequence.fasta")[[1]]

# ---------------------------
# Extract All K-mers and Their Positions
# ---------------------------

# Convert the sequence to character

# Function to analyze sequence and generate consensus
get_consensus <- function(s, kl, ml) {
  # Convert sequence to character and find most frequent kmers
  sequence_chars <- as.character(s)
  sequence_length <- nchar(sequence_chars)

  # Generate all k-mers in order
  all_kmers <- substring(
    sequence_chars, 1:(sequence_length - kl + 1),
    kl:sequence_length
  )

  # Get the most frequent k-mer(s)
  most_freq_kmers <- names(sort(table(all_kmers), decreasing = TRUE)[1])

  # Find positions where the most frequent k-mers occur
  positions <- which(all_kmers %in% most_freq_kmers)

  # Calculate start and end positions for each sequence
  start_positions <- positions
  end_positions <- pmin(positions + ml - 1, sequence_length)

  # Extract sequences using vectorized operations
  extracted_sequences <- DNAStringSet(Views(s,
    start = start_positions, end = end_positions
  ))

  # Write sequences to temporary FASTA file
  tmp_fasta <- tempfile(fileext = ".fasta")
  writeXStringSet(extracted_sequences, filepath = tmp_fasta)

  # Perform multiple sequence alignment with MAFFT
  tmp_aligned <- tempfile(fileext = ".fasta")
  mafft_command <- paste(
    "/home/petr/data/miniforge3/envs/tidecluster_1.5/bin/mafft --auto",
    tmp_fasta, ">", tmp_aligned
  )
  system(mafft_command)

  # Load aligned sequences and generate consensus
  aligned_sequences <- readDNAStringSet(tmp_aligned)
  cm <- consensusMatrix(aligned_sequences)

  # Get most frequent base at each position and remove gaps
  consensus_seq <- apply(cm, 2, function(x) {
    names(sort(x, decreasing = TRUE)[1])
  })
  consensus_seq <- consensus_seq[consensus_seq != "-"] |> paste(collapse = "")

  # Clean up temporary files
  unlink(c(tmp_fasta, tmp_aligned))

  return(consensus_seq)
}

get_consensus(sequence[1], k[1  ])

sequence_chars <- as.character(sequence)
sequence_length <- nchar(sequence_chars)

# Generate all k-mers in order
all_kmers <- substring(sequence_chars, 1:(sequence_length - k + 1), k:sequence_length)


# Get the most frequent k-mer(s)
most_freq_kmers <- names(sort(table(all_kmers), decreasing = TRUE)[1])


# ---------------------------
# Find Positions of the Most Frequent K-mers
# ---------------------------

# Find positions where the most frequent k-mers occur
positions <- which(all_kmers %in% most_freq_kmers)

# ---------------------------
# Extract Sequences Starting with the Most Frequent K-mers
# ---------------------------

# Calculate start and end positions for each sequence
start_positions <- positions
end_positions <- pmin(positions + L - 1, sequence_length)

# Extract sequences using vectorized operations
extracted_sequences <- DNAStringSet(Views(sequence, start = start_positions, end = end_positions))




# Print the number of sequences extracted
cat("Number of sequences extracted:", length(extracted_sequences), "\n")

# ---------------------------
# Write Sequences to FASTA File
# ---------------------------

# Define the output FASTA file name
fasta_file <- "tmp/extracted_sequences.fasta"

# Write the sequences to the FASTA file
writeXStringSet(extracted_sequences, filepath = fasta_file, format = "fasta")

# ---------------------------
# Perform Multiple Sequence Alignment with MAFFT
# ---------------------------

# Define the output alignment file name
alignment_file <- "tmp/aligned_sequences.fasta"

# Construct the MAFFT command
# Adjust the path to MAFFT if necessary
mafft_command <- paste("/home/petr/data/miniforge3/envs/tidecluster_1.5/bin/mafft --auto", fasta_file, ">", alignment_file)

# Execute the MAFFT command
system(mafft_command)

# ---------------------------
# Load the MSA and Extract Consensus Sequence
# ---------------------------

# Load the aligned sequences
aligned_sequences <- readDNAStringSet(alignment_file, format = "fasta")

# Generate the consensus sequence
CM <- consensusMatrix(aligned_sequences)

# get the most frequent base at each position
consensus_seq <- apply(CM, 2, function(x) names(sort(x, decreasing = TRUE)[1]))

# remove all gaps
consensus_seq <- consensus_seq[consensus_seq != "-"] |> paste(collapse = "")


# Print the consensus sequence
cat("Consensus Sequence:\n")
cat(consensus_seq, "\n")

# ---------------------------
# Save the Consensus Sequence to a File
# ---------------------------

# Define the output consensus file name
consensus_file <- "consensus_sequence.fasta"

# Create a DNAStringSet for the consensus
consensus_DNAStringSet <- DNAStringSet(consensus_seq)
names(consensus_DNAStringSet) <- "Consensus_Sequence"

# Write the consensus sequence to a FASTA file
writeXStringSet(consensus_DNAStringSet, filepath = consensus_file, format = "fasta")
