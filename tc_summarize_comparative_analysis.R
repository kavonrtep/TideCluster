#!/usr/bin/env Rscript

# TC Comparative Analysis Summary and Visualization Script
# Generates HTML report with interactive visualizations of TRC comparative analysis results

library(optparse)
library(jsonlite)

# Get the directory where this script is located
script_path <- getwd()

# Function to parse satellite families data from directory
parse_satellite_families <- function(input_dir) {
  cat("Reading satellite families data from directory:", input_dir, "\n")
  
  # Find TSV file in directory
  tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
  tsv_files <- tsv_files[grepl("satellite_families", basename(tsv_files))]
  
  if (length(tsv_files) == 0) {
    stop("No satellite families TSV file found in directory: ", input_dir)
  }
  
  input_file <- tsv_files[1]
  cat("Using TSV file:", input_file, "\n")
  
  # Read the TSV file
  data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     check.names = FALSE, na.strings = "", quote = "")
  
  # Identify sample columns (those without suffixes)
  all_cols <- colnames(data)
  sample_cols <- all_cols[!grepl("_length$|_annot$|_prevalent_annot$|^Satellite_family$|^prevalent_annot$|^group_id$", all_cols)]
  samples <- sample_cols
  
  cat("Detected samples:", paste(samples, collapse = ", "), "\n")
  cat("Total families:", nrow(data), "\n")
  
  # Check for GFF3 directory
  gff3_dir <- file.path(input_dir, "gff3")
  if (!dir.exists(gff3_dir)) {
    cat("Warning: No gff3 directory found. Karyotype visualization will be disabled.\n")
    gff3_dir <- NULL
  }
  
  return(list(
    data = data,
    samples = samples,
    n_families = nrow(data),
    gff3_dir = gff3_dir,
    input_dir = input_dir
  ))
}

# Function to calculate overview statistics
calculate_overview_stats <- function(parsed_data) {
  data <- parsed_data$data
  samples <- parsed_data$samples
  
  overview_stats <- data.frame(
    Sample = samples,
    stringsAsFactors = FALSE
  )
  
  for (sample in samples) {
    # Count TRCs (non-empty cells with TRC content)
    trc_col <- sample
    length_col <- paste0(sample, "_length")
    
    # Count non-empty TRC entries
    trc_entries <- data[[trc_col]]
    trc_entries <- trc_entries[!is.na(trc_entries) & trc_entries != ""]
    
    # Count total TRCs (split by comma)
    total_trcs <- sum(sapply(trc_entries, function(x) length(strsplit(x, ", ")[[1]])))
    
    # Count families present in this sample
    families_present <- sum(!is.na(trc_entries) & trc_entries != "")
    
    # Count unique families (families present only in this sample)
    unique_families <- 0
    for (i in 1:nrow(data)) {
      present_samples <- 0
      for (s in samples) {
        if (!is.na(data[i, s]) && data[i, s] != "") {
          present_samples <- present_samples + 1
        }
      }
      if (present_samples == 1 && !is.na(data[i, sample]) && data[i, sample] != "") {
        unique_families <- unique_families + 1
      }
    }
    
    # Calculate total length
    if (length_col %in% colnames(data)) {
      total_length <- sum(data[[length_col]], na.rm = TRUE)
      # Convert to appropriate units (Mbp if > 1M, kbp if > 1K)
      if (total_length > 1000000) {
        length_display <- sprintf("%.2f Mbp", total_length / 1000000)
      } else if (total_length > 1000) {
        length_display <- sprintf("%.2f kbp", total_length / 1000)
      } else {
        length_display <- sprintf("%d bp", total_length)
      }
    } else {
      total_length <- 0
      length_display <- "0 bp"
    }
    
    # Add to overview stats
    overview_stats[overview_stats$Sample == sample, "Number_of_TRCs"] <- total_trcs
    overview_stats[overview_stats$Sample == sample, "Number_of_Families"] <- families_present
    overview_stats[overview_stats$Sample == sample, "Number_of_Unique_Families"] <- unique_families
    overview_stats[overview_stats$Sample == sample, "Total_Length"] <- length_display
    overview_stats[overview_stats$Sample == sample, "Total_Length_bp"] <- total_length
  }
  
  return(overview_stats)
}

# Function to calculate shared families matrix
calculate_shared_families <- function(parsed_data) {
  data <- parsed_data$data
  samples <- parsed_data$samples
  
  # Create matrix for shared families
  shared_matrix <- matrix(0, nrow = length(samples), ncol = length(samples))
  rownames(shared_matrix) <- samples
  colnames(shared_matrix) <- samples
  
  # Calculate shared families between each pair of samples
  for (i in 1:length(samples)) {
    for (j in 1:length(samples)) {
      sample1 <- samples[i]
      sample2 <- samples[j]
      
      if (i == j) {
        # Diagonal: count families present in this sample
        families_present <- sum(!is.na(data[[sample1]]) & data[[sample1]] != "")
        shared_matrix[i, j] <- families_present
      } else {
        # Off-diagonal: count shared families
        shared_count <- 0
        for (k in 1:nrow(data)) {
          present_in_1 <- !is.na(data[k, sample1]) && data[k, sample1] != ""
          present_in_2 <- !is.na(data[k, sample2]) && data[k, sample2] != ""
          if (present_in_1 && present_in_2) {
            shared_count <- shared_count + 1
          }
        }
        shared_matrix[i, j] <- shared_count
      }
    }
  }
  
  return(shared_matrix)
}

# Function to prepare detailed family data
prepare_detailed_family_data <- function(parsed_data, karyotype_data = NULL) {
  data <- parsed_data$data
  samples <- parsed_data$samples
  
  # Create detailed family matrix
  detailed_data <- list()
  
  for (i in 1:nrow(data)) {
    family_id <- data[i, "Satellite_family"]
    prevalent_annot <- ifelse("prevalent_annot" %in% colnames(data), 
                             data[i, "prevalent_annot"], "")
    if (is.na(prevalent_annot)) prevalent_annot <- ""
    
    family_row <- list(
      family_id = family_id,
      prevalent_annot = prevalent_annot
    )
    
    # Add data for each sample
    for (sample in samples) {
      trc_col <- sample
      length_col <- paste0(sample, "_length")
      
      # Get TRC count
      trc_entry <- data[i, trc_col]
      if (is.na(trc_entry) || trc_entry == "") {
        trc_count <- 0
        total_length <- 0
      } else {
        trc_count <- length(strsplit(trc_entry, ", ")[[1]])
        total_length <- ifelse(length_col %in% colnames(data), 
                              data[i, length_col], 0)
        if (is.na(total_length)) total_length <- 0
      }
      
      # Count arrays from GFF3 data (number of annotated features for this family in this sample)
      array_count <- 0
      if (!is.null(karyotype_data) && sample %in% names(karyotype_data)) {
        sample_karyotype <- karyotype_data[[sample]]
        if (!is.null(sample_karyotype$contigs)) {
          family_id_str <- paste0("SF_", sprintf("%04d", family_id))
          for (contig_name in names(sample_karyotype$contigs)) {
            contig_data <- sample_karyotype$contigs[[contig_name]]
            if (!is.null(contig_data$satellites)) {
              array_count <- array_count + sum(contig_data$satellites$family == family_id_str)
            }
          }
        }
      }
      
      family_row[[paste0(sample, "_trc_count")]] <- trc_count
      family_row[[paste0(sample, "_length")]] <- total_length
      family_row[[paste0(sample, "_array_count")]] <- array_count
    }
    
    detailed_data[[i]] <- family_row
  }
  
  return(detailed_data)
}

# Function to perform hierarchical clustering of samples and return family order by index
cluster_families <- function(parsed_data) {
  data <- parsed_data$data
  samples <- parsed_data$samples
  
  # Create matrix of total lengths for clustering (samples as rows, families as columns)
  length_matrix <- matrix(0, nrow = length(samples), ncol = nrow(data))
  rownames(length_matrix) <- samples
  colnames(length_matrix) <- data$Satellite_family
  
  for (i in 1:length(samples)) {
    for (j in 1:nrow(data)) {
      sample <- samples[i]
      length_col <- paste0(sample, "_length")
      if (length_col %in% colnames(data)) {
        length_val <- data[j, length_col]
        length_matrix[i, j] <- ifelse(is.na(length_val), 0, length_val)
      }
    }
  }
  
  # Perform hierarchical clustering on samples if we have more than 1 sample
  if (length(samples) > 1) {
    # Use Euclidean distance and ward.D2 method
    dist_matrix <- dist(length_matrix, method = "euclidean")
    hc <- hclust(dist_matrix, method = "ward.D2")
    sample_order <- hc$order
  } else {
    sample_order <- 1
  }
  
  # Sort families by their Satellite_family index (ascending order)
  family_order <- order(data$Satellite_family)
  
  return(list(
    sample_order = sample_order,
    family_order = family_order
  ))
}

# Function to read GFF3 files and extract karyotype data
read_gff3_karyotype_data <- function(parsed_data) {
  if (is.null(parsed_data$gff3_dir)) {
    cat("No GFF3 directory found. Skipping karyotype data extraction.\n")
    return(NULL)
  }
  
  cat("Reading GFF3 files for karyotype data...\n")
  
  samples <- parsed_data$samples
  gff3_files <- list.files(parsed_data$gff3_dir, pattern = "\\.gff3$", full.names = TRUE)
  
  if (length(gff3_files) == 0) {
    cat("No GFF3 files found in directory:", parsed_data$gff3_dir, "\n")
    return(NULL)
  }
  
  karyotype_data <- list()
  
  for (gff_file in gff3_files) {
    # Extract sample name from filename (assuming format: SAMPLE_tc_annotated.gff3)
    sample_name <- gsub("_tc_annotated\\.gff3$", "", basename(gff_file))
    
    # Check if this sample is in our samples list
    if (!sample_name %in% samples) {
      next
    }
    
    cat("Processing GFF3 file for sample:", sample_name, "\n")
    
    # Read GFF3 file
    gff_data <- read.table(gff_file, sep = "\t", stringsAsFactors = FALSE, 
                          comment.char = "#", header = FALSE,
                          col.names = c("seqname", "source", "feature", "start", "end", 
                                      "score", "strand", "frame", "attribute"))
    
    # Extract chromosome/contig information and satellite family data
    contigs <- list()
    
    for (i in 1:nrow(gff_data)) {
      contig <- gff_data$seqname[i]
      start_pos <- gff_data$start[i]
      end_pos <- gff_data$end[i]
      attributes <- gff_data$attribute[i]
      
      # Extract Satellite_family from attributes
      sf_match <- regmatches(attributes, regexpr("Satellite_family=SF_[0-9]+", attributes))
      if (length(sf_match) > 0) {
        sf_id <- gsub("Satellite_family=", "", sf_match)
        
        # Initialize contig data if not exists
        if (!contig %in% names(contigs)) {
          contigs[[contig]] <- list(
            satellites = data.frame(
              start = integer(0),
              end = integer(0),
              family = character(0),
              stringsAsFactors = FALSE
            )
          )
        }
        
        # Add satellite family position
        contigs[[contig]]$satellites <- rbind(contigs[[contig]]$satellites,
                                            data.frame(start = start_pos, 
                                                     end = end_pos,
                                                     family = sf_id,
                                                     stringsAsFactors = FALSE))
      }
    }
    
    # Calculate contig lengths and sort contigs by order in GFF3 file
    contig_info <- data.frame(
      contig = character(0),
      length = integer(0),
      order = integer(0),
      stringsAsFactors = FALSE
    )
    
    contig_order <- 1
    for (contig_name in unique(gff_data$seqname)) {
      max_pos <- max(gff_data$end[gff_data$seqname == contig_name])
      contig_info <- rbind(contig_info, 
                          data.frame(contig = contig_name, 
                                   length = max_pos,
                                   order = contig_order,
                                   stringsAsFactors = FALSE))
      contig_order <- contig_order + 1
    }
    
    karyotype_data[[sample_name]] <- list(
      contigs = contigs,
      contig_info = contig_info
    )
  }
  
  cat("Karyotype data extracted for", length(karyotype_data), "samples\n")
  return(karyotype_data)
}

# Function to generate HTML report
generate_html_report <- function(parsed_data, overview_stats, shared_matrix, 
                                detailed_data, clustering_result, karyotype_data, output_file) {
  
  # Create output directory and js subdirectory if they don't exist
  output_dir <- dirname(output_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  js_dir <- file.path(output_dir, "js")
  dir.create(js_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Reorder detailed data according to family index (ascending)
  detailed_data_ordered <- detailed_data[clustering_result$family_order]
  
  # Reorder samples according to hierarchical clustering
  samples_ordered <- parsed_data$samples[clustering_result$sample_order]
  
  # Reorder shared matrix according to sample clustering
  shared_matrix_ordered <- shared_matrix[clustering_result$sample_order, clustering_result$sample_order]
  
  # Prepare data for JavaScript
  js_data <- list(
    samples = samples_ordered,
    overview_stats = overview_stats,
    shared_matrix = shared_matrix_ordered,
    detailed_families = detailed_data_ordered,
    n_families = parsed_data$n_families,
    sample_order = clustering_result$sample_order,
    karyotype_data = karyotype_data
  )
  
  # Convert to JSON
  json_data <- toJSON(js_data, auto_unbox = FALSE, pretty = TRUE)
  
  # Read HTML template
  template_path <- file.path(dirname(getwd()), "TideCluster", "html", "report_template.html")
  if (!file.exists(template_path)) {
    # Try relative path
    template_path <- file.path("html", "report_template.html")
  }
  
  if (file.exists(template_path)) {
    html_content <- readLines(template_path, warn = FALSE)
    html_content <- paste(html_content, collapse = "\n")
    
    # Replace placeholders
    html_content <- gsub("\\{\\{TIMESTAMP\\}\\}", Sys.time(), html_content)
    html_content <- gsub("\\{\\{DATA_JSON\\}\\}", json_data, html_content)
  } else {
    stop("HTML template file not found: ", template_path)
  }
  
  # Write HTML file
  writeLines(html_content, output_file)
  cat("HTML report written to:", output_file, "\n")
}

# Function to create JavaScript files
create_javascript_files <- function(output_dir) {
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  js_dir <- file.path(output_dir, "js")
  dir.create(js_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Copy JavaScript and CSS files from templates
  template_dir <- file.path(script_path, "html")
  
  if (!dir.exists(template_dir)) {
    stop("Template directory not found: ", template_dir)
  }
  
  # List of files to copy
  template_files <- c(
    "main.js",
    "utils.js", 
    "overview-table.js",
    "shared-matrix.js",
    "detailed-matrix.js",
    "karyotype.js",
    "styles.css"
  )
  
  for (file in template_files) {
    src_file <- file.path(template_dir, file)
    if (file.exists(src_file)) {
      if (file == "styles.css") {
        dest_file <- file.path(js_dir, file)
      } else {
        dest_file <- file.path(js_dir, file)
      }
      file.copy(src_file, dest_file, overwrite = TRUE)
      cat("Copied template file:", file, "\n")
    } else {
      cat("Warning: Template file not found:", src_file, "\n")
    }
  }
  
  cat("JavaScript and CSS files created in:", js_dir, "\n")
}

# Main function
main <- function(opt) {
  cat("Starting TRC comparative analysis summarization...\n")
  
  # Parse input data
  parsed_data <- parse_satellite_families(opt$input)
  
  # Calculate overview statistics
  cat("Calculating overview statistics...\n")
  overview_stats <- calculate_overview_stats(parsed_data)
  
  # Calculate shared families matrix
  cat("Calculating shared families matrix...\n")
  shared_matrix <- calculate_shared_families(parsed_data)
  
  # Read karyotype data from GFF3 files first
  karyotype_data <- read_gff3_karyotype_data(parsed_data)
  
  # Prepare detailed family data (now with karyotype data for array counting)
  cat("Preparing detailed family data...\n")
  detailed_data <- prepare_detailed_family_data(parsed_data, karyotype_data)
  
  # Perform hierarchical clustering of samples and sort families by index
  cat("Performing hierarchical clustering of samples and sorting families by index...\n")
  clustering_result <- cluster_families(parsed_data)
  
  # Create output directory
  output_dir <- dirname(opt$output)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create JavaScript and CSS files
  cat("Creating JavaScript and CSS files...\n")
  create_javascript_files(output_dir)
  
  # Generate HTML report
  cat("Generating HTML report...\n")
  generate_html_report(parsed_data, overview_stats, shared_matrix, 
                      detailed_data, clustering_result, karyotype_data, opt$output)
  
  cat("Report generation completed successfully!\n")
  cat("Open", opt$output, "in your web browser to view the interactive report.\n")
}

# Parse command line arguments
option_list <- list(
  make_option(
    c("-i", "--input"), type="character", default=NULL,
    help="Input directory containing TSV file and GFF3 subdirectory from tc_comparative_analysis.R script (required)"),
  make_option(
    c("-o", "--output"), type="character", default="trc_comparative_report.html",
    help="Output HTML file path [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input directory is mandatory. Use -i or --input to specify it.")
}

if (!dir.exists(opt$input)) {
  stop("Input directory does not exist: ", opt$input)
}

# Run main analysis
main(opt)