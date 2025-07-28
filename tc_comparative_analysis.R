#!/usr/bin/env Rscript
# MMseqs2 similarity search wrapper function
mmseqs2_search <- function(sequences,
                           search_type = 3,
                           sensitivity = 7.5,
                           threads = 1,
                           format_output = "query,target,pident,qcov,tcov",
                           additional_params = "",
                           mmseqs2_path = "mmseqs",
                           min_coverage = 0.2,
                           min_identity = 70,
                           awk_filter = NULL
) {

  # Create temporary directory for all operations
  temp_dir <- tempfile(pattern = "mmseqs2_search_")
  dir.create(temp_dir, recursive = TRUE)

  # Define file paths
  input_fasta <- file.path(temp_dir, "input.fasta")
  output_file <- file.path(temp_dir, "pairs.m8")
  tmp_dir <- file.path(temp_dir, "tmp")

  # Create tmp directory
  dir.create(tmp_dir, recursive = TRUE)

  tryCatch({
    # Write input sequences to temporary file
    writeXStringSet(sequences, input_fasta, format = "fasta")

    # Build MMseqs2 command with filtering options
    mmseqs_cmd <- paste(mmseqs2_path, " easy-search",
                        input_fasta,
                        input_fasta,  # Search against itself
                        output_file,
                        tmp_dir,
                        "--search-type", search_type,
                        "-a",  # Include amino acid sequences in output
                        paste0("--format-output \"", format_output, "\""),
                        "-s", sensitivity, "-v 0", # No verbose output
                        "--threads", threads,
                        "--min-seq-id", min_identity / 100,  # Convert percentage to fraction
                        additional_params)

    cat("Running MMseqs2 search command:\n", mmseqs_cmd, "\n")

    # Execute MMseqs2
    result_code <- system(mmseqs_cmd)

    if (result_code != 0) {
      stop("MMseqs2 search command failed with exit code: ", result_code)
    }

    # Apply additional AWK filtering if specified
    if (!is.null(awk_filter)) {
      cat("Applying additional AWK filtering...\n")
      temp_filtered_file <- file.path(temp_dir, "pairs_filtered.m8")
      
      # Build AWK command for additional filtering
      awk_cmd <- paste("awk", paste0("'", awk_filter, "'"), 
                      output_file, ">", temp_filtered_file)
      
      cat("Running AWK filter command:\n", awk_cmd, "\n")
      
      awk_result_code <- system(awk_cmd)
      
      if (awk_result_code != 0) {
        stop("AWK filtering command failed with exit code: ", awk_result_code)
      }
      
      # Replace original output file with filtered version
      output_file <- temp_filtered_file
    }

    # Read results
    if (!file.exists(output_file)) {
      stop("Output file not created: ", output_file)
    }

    if (file.size(output_file) == 0) {
      warning("Output file is empty - no matches found")
      # Return empty data frame with correct column names
      column_names <- strsplit(format_output, ",")[[1]]
      empty_df <- data.frame(matrix(ncol = length(column_names), nrow = 0))
      colnames(empty_df) <- column_names
      return(empty_df)
    }

    # Read the results file
    search_results <- read.table(output_file,
                                header = FALSE,
                                stringsAsFactors = FALSE,
                                sep = "\t",
                                comment.char = "",
                                quote = "")

    # Set column names based on format_output
    column_names <- strsplit(format_output, ",")[[1]]

    if (ncol(search_results) != length(column_names)) {
      warning("Number of columns in output (", ncol(search_results),
              ") doesn't match expected format (", length(column_names), ")")
      # Use generic column names
      colnames(search_results) <- paste0("col", 1:ncol(search_results))
    } else {
      colnames(search_results) <- column_names
    }

    cat("Read", nrow(search_results), "similarity pairs\n")

    # Convert numeric columns to appropriate types
    numeric_cols <- c("pident", "alnlen", "qcov", "tcov", "evalue", "bits")
    for (col in numeric_cols) {
      if (col %in% colnames(search_results)) {
        search_results[[col]] <- as.numeric(search_results[[col]])
      }
    }

    return(search_results)

  }, finally = {
    # Clean up temporary directory
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      cat("Cleaned up temporary directory:", temp_dir, "\n")
    }
  })
}


cluster_trc_sequences <- function(tc_seq, th_seq, mmseqs2_path = NULL,
                                  min_coverage = 0.2, min_identity = 70,
                                  output_directory = NULL, ncpu = 5
) {

  # Input validation
  if (!inherits(tc_seq, "DNAStringSet") || !inherits(th_seq, "DNAStringSet")) {
    stop("Both tc_seq and th_seq must be DNAStringSet objects")
  }
  if (length(tc_seq) == 0 || length(th_seq) == 0) {
    stop("Input sequences cannot be empty")
  }
  if (!is.numeric(min_coverage) || min_coverage < 0 || min_coverage > 1) {
    stop("min_coverage must be numeric between 0 and 1")
  }
  if (!is.numeric(min_identity) || min_identity < 0 || min_identity > 100) {
    stop("min_identity must be numeric between 0 and 100")
  }
  # Load required libraries
  required_packages <- c("Biostrings", "igraph")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Required packages not found:", paste(missing_packages, collapse = ", ")))
  }
  # Ensure unique names for TC sequences
  names(tc_seq) <- paste(names(tc_seq), seq_along(tc_seq), sep = "_")

  # Extract TRC IDs
  tc_trc_id <- gsub("#.+", "", names(tc_seq))
  th_trc_id <- gsub("_rep.+", "", names(th_seq))

  # Combine all unique TRC IDs
  trc_all <- unique(c(tc_trc_id, th_trc_id))

  # Combine sequences, avoiding duplicates based on transcript IDs
  th_tc_seq <- c(tc_seq, th_seq[!th_trc_id %in% tc_trc_id])

  message("Combined ", length(th_tc_seq), " sequences for clustering")

  # add also reverse complement sequences
  th_tc_seq <- c(th_tc_seq, reverseComplement(th_tc_seq))

  # Build AWK filter for coverage filtering (since MMseqs2 doesn't have --min-cov)
  # format_output = "query,target,pident,qcov,tcov"
  # Column indices: qcov=4, tcov=5, pident=3
  awk_coverage_filter <- sprintf("(($4 >= %f) || ($5 >= %f))", min_coverage, min_coverage)
  
  # Run MMseqs2 search with filtering parameters
  tryCatch({
    if (is.null(mmseqs2_path)) {
      df <- mmseqs2_search(th_tc_seq, threads = ncpu, 
                          min_identity = min_identity,
                          awk_filter = awk_coverage_filter)
    } else {
      df <- mmseqs2_search(th_tc_seq, threads = ncpu, 
                          mmseqs2_path = mmseqs2_path,
                          min_identity = min_identity,
                          awk_filter = awk_coverage_filter)
    }
  }, error = function(e) {
    stop("MMseqs2 execution failed. Please check that MMseqs2 is properly installed and accessible: ", e$message)
  })

  # Calculate maximum coverage
  df$max_cov <- pmax(df$qcov, df$tcov)

  # Extract  IDs from query and target
  df$trc_q <- gsub("_rep.+", "", gsub("#.+", "", df$query))
  df$trc_t <- gsub("_rep.+", "", gsub("#.+", "", df$target))



  # Filter based on coverage and identity thresholds, or same transcript ID
  message("MMseqs2 output before filtering: ", nrow(df), " rows")
  df_pass <- df[(df$max_cov >= min_coverage & df$pident >= min_identity) |
                df$trc_q == df$trc_t, ]
  message("MMseqs2 output after filtering: ", nrow(df_pass), " rows")

  # Calculate edge weights
  df_pass$weight <- df_pass$max_cov * df_pass$pident

  # Build graph and perform community detection
  g <- igraph::simplify(
    igraph::graph_from_data_frame(df_pass[, c("query", "target", "weight")],
                                  directed = FALSE)
  )



  fg <- cluster_fast_greedy(g)
  fg_groups <- igraph::groups(fg)

  # add fastgreedy groups to graph as vertex attribute
  V(g)$group <- as.integer(factor(igraph::membership(fg)))
  # Add species code attribute - this is the prefix of the sequence name before the first colon
  V(g)$code <- gsub(":.+", "", names(V(g)))
  # add TRC ID attribute - this is the transcript ID without prefix and the #rep or _rep suffix
  V(g)$trc_id <- gsub(".+:", "",
                      gsub("_rep.+", "",
                           gsub("#.+", "", names(V(g)))))

  if (!is.null(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    write.table(df, file = file.path(output_directory, "mmseqs2_results.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    saveRDS(g, file = file.path(output_directory, "trc_graph.rds"))
    write_graph(g, file = file.path(output_directory, "trc_graph.graphml"),
              format = "graphml")
  }


  # Extract transcript IDs for each group
  fg_groups_trc_id <- lapply(fg_groups, function(x) {
    unique(gsub("_rep.+", "", gsub("#.+", "", x)))
  })

  # Convert to data frame format
  trc_groups <- do.call(rbind, lapply(seq_along(fg_groups_trc_id), function(i) {
    data.frame(trc_id = fg_groups_trc_id[[i]],
               group_id = i,
               stringsAsFactors = FALSE)
  }))

  message("Initial clustering produced ", max(trc_groups$group_id), " groups")

  # Handle transcript IDs appearing in multiple groups
  trc_duplicated <- trc_groups$trc_id[duplicated(trc_groups$trc_id)]

  if (length(trc_duplicated) > 0) {
    message("Merging groups for ", length(trc_duplicated), " duplicated TRC IDs")

    for (trc in trc_duplicated) {
      # Find all groups for this transcript ID
      grps <- unique(trc_groups$group_id[trc_groups$trc_id == trc])

      if (length(grps) > 1) {
        # Merge into the smallest group ID
        new_group_id <- min(grps)
        trc_groups$group_id[trc_groups$group_id %in% grps] <- new_group_id
      }
    }
  }
  # Remove duplicated rows and renumber groups starting from 1
  trc_groups <- trc_groups[!duplicated(trc_groups), ]
  trc_groups$group_id <- as.integer(factor(trc_groups$group_id))

  message("Final clustering produced ", max(trc_groups$group_id), " groups covering ",
          nrow(trc_groups), " TRC IDs")

  return(trc_groups)
}
pivot_trc_data <- function(df, prefix) {
  df2 <- data.frame(
    prefix = gsub(":.+", "", df$trc_id),
    trc_id = gsub(".+:", "", df$trc_id),
    group_id = df$group_id,
    stringsAsFactors = FALSE
  )
  # make wide format - prefix as columns, group_id as rows, values are concatenated trc_id
  result_df <- reshape2::dcast(df2, group_id ~ prefix, value.var = "trc_id", fun.aggregate = function(x) paste(unique(x), collapse = ", "))
  # make also lists
  result_list <- lapply(seq_along(colnames(result_df)[-1]), function(i) {
    res <- result_df[[i + 1]]
    trc <- strsplit(res, ", ")
    # remove empty strings within lists with NULL
    names(trc) <- result_df$group_id
    return(trc)
  })
  # set same order as prefixes
  correct_order <- match(prefix, colnames(result_df)[-1])
  result_list <- result_list[correct_order]
  result_df <- result_df[, c("group_id", colnames(result_df)[-1][correct_order])]
  return(list(
    df = result_df,
    list = result_list
  ))

}
read_keyvalue_file <- function(file_path, sep = "\t", comment = "#",
                              skip_empty = TRUE, auto_convert = TRUE) {

  # Input validation
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  if (!is.character(sep) || length(sep) != 1) {
    stop("sep must be a single character string")
  }

  # Read all lines from file
  lines <- readLines(file_path, warn = FALSE)

  # Remove empty lines if requested
  if (skip_empty) {
    lines <- lines[nzchar(trimws(lines))]
  }

  # Remove comment lines
  if (!is.null(comment) && nchar(comment) > 0) {
    lines <- lines[!grepl(paste0("^\\s*", comment), lines)]
  }

  if (length(lines) == 0) {
    warning("No data lines found in file")
    return(list())
  }

  # Initialize result list
  result <- list()

  # Process each line
  for (i in seq_along(lines)) {
    line <- lines[i]

    # Skip empty lines
    if (!nzchar(trimws(line))) {
      next
    }

    # Split line by separator
    parts <- strsplit(line, sep, fixed = TRUE)[[1]]

    # Remove leading/trailing whitespace
    parts <- trimws(parts)

    # Must have at least identifier
    if (length(parts) < 1) {
      warning("Line ", i, " is empty or invalid, skipping")
      next
    }

    # First part is the identifier
    identifier <- parts[1]

    # Check if identifier is empty
    if (!nzchar(identifier)) {
      warning("Line ", i, " has empty identifier, skipping")
      next
    }

    # Remaining parts are key-value pairs
    kv_parts <- parts[-1]

    # Initialize sub-list
    sub_list <- list()

    # Process key-value pairs
    if (length(kv_parts) > 0) {
      # Check if we have an even number of key-value elements
      if (length(kv_parts) %% 2 != 0) {
        warning("Line ", i, " has odd number of key-value elements. Last element ignored.")
        kv_parts <- kv_parts[-length(kv_parts)]
      }

      # Process pairs
      if (length(kv_parts) > 0) {
        # Create sequence of indices for keys (odd positions) and values (even positions)
        key_indices <- seq(1, length(kv_parts), by = 2)
        value_indices <- seq(2, length(kv_parts), by = 2)

        keys <- kv_parts[key_indices]
        values <- kv_parts[value_indices]

        # Check for duplicate keys within the same row
        if (any(duplicated(keys))) {
          dup_keys <- keys[duplicated(keys)]
          warning("Line ", i, " has duplicate keys: ", paste(dup_keys, collapse = ", "),
                  ". Only the last value will be kept.")
        }

        # Auto-convert values if requested
        if (auto_convert) {
          values <- lapply(values, function(x) {
            # Try numeric conversion
            if (!is.na(suppressWarnings(as.numeric(x)))) {
              return(as.numeric(x))
            }
            # Try logical conversion
            if (tolower(x) %in% c("true", "false", "t", "f")) {
              return(as.logical(toupper(substr(x, 1, 1))))
            }
            # Keep as character
            return(x)
          })
        }

        # Create named list
        names(values) <- keys
        sub_list <- as.list(values)
      }
    }

    # Add to result (handle duplicate identifiers)
    if (identifier %in% names(result)) {
      warning("Duplicate identifier '", identifier, "' found. Previous entry will be overwritten.")
    }

    result[[identifier]] <- sub_list
  }

  message("Successfully read ", length(result), " entries from ", file_path)
  return(result)
}

make_multimer <- function(seq_set, k = 2){
  N <-  names(seq_set)
  seq_char <- as.character(seq_set)
  seq_char_multy <- sapply(seq_char, function(x, k) {
    rep_seq <- paste(rep(x, k), collapse = "")
  }, k = k)
  seq_set_multy <- DNAStringSet(seq_char_multy)
  names(seq_set_multy) <- N
  seq_set_multy <- unique(seq_set_multy)
  return(seq_set_multy)
}
# read and concatenate all input files from input_tc_dirs
get_seq_files <- function(input_dirs, prefix){
  # TODO - handle tc_ prefixes
  tarean_consensus_path <- "tc_consensus_dimer_library.fasta"  # take first TRC
  consensus_groups_path <- "tc_consensus/consensus_sequences_all.fasta"
  fasta_tc <- paste0(input_dirs, "/", tarean_consensus_path)
  fasta_th <- paste0(input_dirs, "/", consensus_groups_path)
  s_tc <- sapply(fasta_tc, readDNAStringSet)
  s_th <- sapply(fasta_th, readDNAStringSet)
  s_tc <- lapply(s_tc, function(x) {
    if (length(x) == 0) {
      return(DNAStringSet())
    } else {
      return(make_multimer(x, k = 1))  # make dimers
    }
  })
    s_th <- lapply(s_th, function(x) {
        if (length(x) == 0) {
        return(DNAStringSet())
        } else {
        return(make_multimer(x, k = 2))  # make tetramers
        }
    })


  # add prefixes to sequence names:
  s_tc <- lapply(seq_along(prefix), function(i) {
    names(s_tc[[i]]) <- paste(prefix[i], names(s_tc[[i]]), sep = ":")
    s_tc[[i]]
  })
  s_th <- lapply(seq_along(prefix), function(i) {
    names(s_th[[i]]) <- paste(prefix[i], names(s_th[[i]]), sep = ":")
    s_th[[i]]
  })
  # concatenate all sequences into one DNAStringSet
  all_tc <- Reduce(c, s_tc)
  all_th <- Reduce(c, s_th)
  return(list(tc = all_tc, th = all_th))
}

cluster_ssrs_sequences <- function(data_list, min_percentage = 10) {

  # Helper function to extract SSR types above threshold
  extract_ssrs_above_threshold <- function(ssrs_string, threshold = min_percentage) {
    # Split by comma and trim whitespace
    ssrs_parts <- trimws(strsplit(ssrs_string, ",")[[1]])

    # Extract SSR type and percentage for each part
    valid_ssrs <- c()
    for (part in ssrs_parts) {
      # Extract percentage using regex
      percentage_match <- regmatches(part, regexpr("\\(([0-9.]+)%", part))
      if (length(percentage_match) > 0) {
        # Extract the numeric percentage
        percentage <- as.numeric(gsub("[(%)]", "", percentage_match))
        if (percentage >= threshold) {
          # Extract SSR type (everything before the percentage)
          ssrs_type <- trimws(gsub("\\s*\\([0-9.]+%.*\\)", "", part))
          valid_ssrs <- c(valid_ssrs, ssrs_type)
        }
      }
    }

    # Return as comma-separated string in original order (by percentage)
    if (length(valid_ssrs) > 0) {
      return(paste(valid_ssrs, collapse = ", "))
    } else {
      return(NA)
    }
  }

  # Process each sample and create mapping
  all_samples <- list()

  for (sample_name in names(data_list)) {
    df <- data_list[[sample_name]]

    # Create the merging column with SSRs above threshold
    df$ssrs_merge <- sapply(df$ssrs_type, extract_ssrs_above_threshold)

    # Keep only rows with valid SSR types (non-NA)
    df <- df[!is.na(df$ssrs_merge), ]

    # Store processed data
    all_samples[[sample_name]] <- df
  }

  # Get all unique SSR patterns (clusters)
  all_patterns <- unique(unlist(lapply(all_samples, function(x) x$ssrs_merge)))
  all_patterns <- all_patterns[!is.na(all_patterns)]

  # Get all sample names
  sample_names <- names(data_list)

  # Create result dataframe
  result_df <- data.frame(
    cluster_index = seq_along(all_patterns),
    stringsAsFactors = FALSE
  )

  # Add columns for each sample's TRC ID and SSR type
  for (sample_name in sample_names) {
    # Initialize columns with NA
    result_df[[paste0(sample_name, "_trc_id")]] <- NA
    result_df[[paste0(sample_name, "_ssrs_type")]] <- NA

    # Fill in the data
    sample_data <- all_samples[[sample_name]]

    for (i in seq_along(all_patterns)) {
      pattern <- all_patterns[i]

      # Find matching TRC in this sample
      matching_rows <- which(sample_data$ssrs_merge == pattern)

      if (length(matching_rows) > 0) {
        # Take the first match if multiple (shouldn't happen with proper data)
        match_idx <- matching_rows[1]
        result_df[i, paste0(sample_name, "_trc_id")] <- sample_data$trc_id[match_idx]
        result_df[i, paste0(sample_name, "_ssrs_type")] <- sample_data$ssrs_type[match_idx]
      }
    }
  }

  # Add the cluster pattern for reference
  result_df$major_pattern <- all_patterns

  # Reorder columns: cluster_index, then alternating trc_id and ssrs_type for each sample, then cluster_pattern
  trc_cols <- paste0(sample_names, "_trc_id")
  ssrs_cols <- paste0(sample_names, "_ssrs_type")

  ordered_columns <- c("cluster_index", trc_cols, ssrs_cols, "major_pattern")
  result_df <- result_df[, ordered_columns, drop = FALSE]
  return(result_df)
}

# Function to calculate total lengths from GFF3 files
calculate_trc_lengths <- function(input_tc_dirs, prefix, tc_clust_path = "tc_clustering.gff3") {
  total_trc_length <- list()

  for (i in seq_along(input_tc_dirs)) {
    gff_path <- file.path(input_tc_dirs[i], tc_clust_path[i])
    if (!file.exists(gff_path)) {
      warning("GFF3 file not found: ", gff_path)
      total_trc_length[[prefix[i]]] <- numeric(0)
      next
    }

    gff <- rtracklayer::import.gff3(gff_path)
    total_trc_length[[prefix[i]]] <- sapply(split(width(gff), gff$Name), sum)
  }

  return(total_trc_length)
}

# Function to calculate group lengths
calculate_group_lengths <- function(grps_pivoted, total_trc_length, prefix) {
  total_grps_length <- list()

  for (i in seq_along(prefix)) {
    total_grps_length[[prefix[i]]] <- sapply(grps_pivoted$list[[i]], function(x) {
      if (length(x) == 0) {
        return(0)
      } else {
        return(sum(total_trc_length[[prefix[i]]][x], na.rm = TRUE))
      }
    })
  }

  total_grps_length <- do.call(cbind, total_grps_length)
  colnames(total_grps_length) <- paste0(colnames(total_grps_length), "_length")

  return(total_grps_length)
}

# Function to process annotations for groups
process_group_annotations <- function(input_tc_dirs, prefix, grps_pivoted, total_trc_length,
                                    total_grps_length, tc_annot_path = "tc_annotation.tsv") {
  annot_prefix_groups <- list()
  annot_prefix_groups_fraction <- list()
  prevalent_prefix_groups_annot <- list()
  if (length(tc_annot_path) == 1){
    tc_annot_path <- rep(tc_annot_path, length(prefix))
  }
  for (i in seq_along(prefix)) {
    annot_file <- file.path(input_tc_dirs[i], tc_annot_path[i])

    if (file.exists(annot_file)) {
      message("Adding annotation from ", tc_annot_path[i], " for prefix ", prefix[i])
      annot <- read_keyvalue_file(annot_file)
    } else {
      message("Annotation file not found: ", annot_file)
      annot <- list()
    }

    annot_prefix_groups[[prefix[i]]] <- list()
    annot_prefix_groups_fraction[[prefix[i]]] <- list()
    prevalent_prefix_groups_annot[[prefix[i]]] <- list()

    # Process each group
    for (sf in seq_along(grps_pivoted$list[[i]])) {
      trc_composition <- grps_pivoted$list[[i]][[sf]]

      if (length(trc_composition) == 0) {
        annot_prefix_groups[[prefix[i]]][[sf]] <- NA
        annot_prefix_groups_fraction[[prefix[i]]][[sf]] <- NA
        prevalent_prefix_groups_annot[[prefix[i]]][[sf]] <- NA
        next
      }

      # Check if any TRC in this group has annotation
      if (any(trc_composition %in% names(annot))) {
        annot_values <- annot[trc_composition[trc_composition %in% names(annot)]]

        # Calculate weighted annotation lengths
        annot_length_recalculated <- numeric()
        for (n in names(annot_values)) {
          if (n %in% names(total_trc_length[[prefix[i]]])) {
            summed_length <- unlist(annot_values[[n]]) * total_trc_length[[prefix[i]]][n]

            for (ann_type in names(summed_length)) {
              if (ann_type %in% names(annot_length_recalculated)) {
                annot_length_recalculated[[ann_type]] <- annot_length_recalculated[[ann_type]] + summed_length[[ann_type]]
              } else {
                annot_length_recalculated[[ann_type]] <- summed_length[[ann_type]]
              }
            }
          }
        }

        annot_prefix_groups[[prefix[i]]][[sf]] <- annot_length_recalculated

        # Calculate fractions
        if (total_grps_length[sf, i] > 0) {
          annot_prefix_groups_fraction[[prefix[i]]][[sf]] <- annot_length_recalculated / total_grps_length[sf, i]

          # Find prevalent annotations (>50%)
          prevalent_annot <- annot_prefix_groups_fraction[[prefix[i]]][[sf]] > 0.5
          if (any(prevalent_annot)) {
            prevalent_prefix_groups_annot[[prefix[i]]][[sf]] <- paste(names(annot_length_recalculated)[prevalent_annot], collapse = ", ")
          } else {
            prevalent_prefix_groups_annot[[prefix[i]]][[sf]] <- NA
          }
        } else {
          annot_prefix_groups_fraction[[prefix[i]]][[sf]] <- NA
          prevalent_prefix_groups_annot[[prefix[i]]][[sf]] <- NA
        }
      } else {
        annot_prefix_groups[[prefix[i]]][[sf]] <- NA
        annot_prefix_groups_fraction[[prefix[i]]][[sf]] <- NA
        prevalent_prefix_groups_annot[[prefix[i]]][[sf]] <- NA
      }
    }
  }

  return(list(
    annot_groups = annot_prefix_groups,
    annot_fractions = annot_prefix_groups_fraction,
    prevalent_annot = prevalent_prefix_groups_annot
  ))
}

# Function to create annotation data frames
create_annotation_dataframes <- function(grps_pivoted, annotation_results, prefix) {

  # Create annotation details data frame
  annot_prefix_groups_df <- data.frame(
    group_id = grps_pivoted$df$group_id,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(prefix)) {
    cln <- paste0(prefix[i], "_annot")
    annot_prefix_groups_df[[cln]] <- sapply(annotation_results$annot_fractions[[i]], function(x) {
      if (is.na(x)[1]) {
        return("")
      } else {
        return(paste(names(x), " (", round(unlist(x) * 100, 1), "%)", sep = "", collapse = ", "))
      }
    })
  }

  # Create prevalent annotation data frame
  prevalent_annot_df <- data.frame(
    group_id = grps_pivoted$df$group_id,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(prefix)) {
    cln <- paste0(prefix[i], "_prevalent_annot")
    prevalent_annot_df[[cln]] <- sapply(annotation_results$prevalent_annot[[i]], function(x) {
      if (is.na(x)) {
        return("")
      } else {
        return(x)
      }
    })
  }

  # Combine prevalent annotations across prefixes
  prevalent_annot_df$prevalent_annot <- apply(prevalent_annot_df[, -1, drop = FALSE], 1, function(x) {
    non_empty <- x[x != ""]
    if (length(non_empty) > 0) {
      return(paste(unique(non_empty), collapse = ", "))
    } else {
      return("")
    }
  })

  return(list(
    annot_df = annot_prefix_groups_df,
    prevalent_df = prevalent_annot_df
  ))
}

# Function to renumber groups by size and finalize data frame
finalize_groups_dataframe <- function(grps_pivoted, total_grps_length, annotation_dfs) {

  # Calculate total length across all prefixes
  total_grps_length_in_all_prefixes <- rowSums(total_grps_length)

  # Add length columns to main data frame
  grps_pivoted$df <- cbind(grps_pivoted$df, total_grps_length)

  # Add annotation columns
  grps_pivoted$df <- cbind(grps_pivoted$df,
                          annotation_dfs$annot_df[,-1],
                          annotation_dfs$prevalent_df[,-1])

  # Renumber groups by size (largest = 1)
  groups_id_renumbered <- rank(-total_grps_length_in_all_prefixes, ties.method = "first")
  grps_pivoted$df$group_id_new <- groups_id_renumbered

  # Sort by new group ID
  new_order <- order(grps_pivoted$df$group_id_new)
  grps_pivoted$df <- grps_pivoted$df[new_order, ]

  return(grps_pivoted)
}

# Main workflow function that combines everything
process_trc_analysis <- function(input_tc_dirs, prefix,tc_code = "tc",
                                 mmseqs2_path = NULL,
                                 output_directory = "tc_comparative_analysis",
                                 ncpu = opt$cpu
) {

  tc_clust_path <- paste0(tc_code, "_clustering.gff3")
  tc_annot_path <- paste0(tc_code, "_annotation.tsv")
  ssrs_table_path <- paste0(tc_code, "_tarean/SSRS_summary.csv")

  # Get sequences and perform clustering
  message("Loading sequences...")
  all_seq <- get_seq_files(input_tc_dirs, prefix)

  message("Clustering sequences...")
  grps <- cluster_trc_sequences(all_seq$tc, all_seq$th, mmseqs2_path = mmseqs2_path,
                                output_directory = output_directory, ncpu = ncpu)
  grps_pivoted <- pivot_trc_data(grps, prefix)


  # Load and process SSRS data
  message("Processing SSRS data...")

  ssrs_tables <- mapply(
  function(dir, tbl_path) {
    ssrs_path <- file.path(dir, tbl_path)
    if (file.exists(ssrs_path)) {
      message("Loading SSRS table from ", ssrs_path)
      read.table(
        ssrs_path,
        header           = FALSE,
        sep              = "\t",
        stringsAsFactors = FALSE,
        col.names        = c("trc_id", "ssrs_type")
      )
    } else {
      message("SSRS table not found at ", ssrs_path)
      data.frame()
    }
  },
  input_tc_dirs,
  ssrs_table_path,
  SIMPLIFY = FALSE
  )





  names(ssrs_tables) <- prefix

  ssrs_groups <- cluster_ssrs_sequences(ssrs_tables, min_percentage = 10)

  # Calculate lengths
  message("Calculating TRC lengths...")
  total_trc_length <- calculate_trc_lengths(input_tc_dirs, prefix, tc_clust_path)
  total_grps_length <- calculate_group_lengths(grps_pivoted, total_trc_length, prefix)

  # Process annotations
  message("Processing annotations...")
  annotation_results <- process_group_annotations(input_tc_dirs, prefix, grps_pivoted,
                                                 total_trc_length, total_grps_length, tc_annot_path)

  # Create annotation data frames
  annotation_dfs <- create_annotation_dataframes(grps_pivoted, annotation_results, prefix)

  # Finalize the main data frame
  message("Finalizing results...")
  grps_pivoted_final <- finalize_groups_dataframe(grps_pivoted, total_grps_length, annotation_dfs)

  return(list(
    groups = grps_pivoted_final,
    ssrs_groups = ssrs_groups,
    total_trc_length = total_trc_length,
    annotation_results = annotation_results
  ))
}
# Function to validate and read input configuration
validate_and_read_input <- function(input_file) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  input_table <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Validate required columns
  required_cols <- c("input_dir", "sample_code", "tidecluster_prefix")
  missing_cols <- setdiff(required_cols, colnames(input_table))

  if (length(missing_cols) > 0) {
    stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
  }

  # Validate that input directories exist
  missing_dirs <- !dir.exists(input_table$input_dir)
  if (any(missing_dirs)) {
    warning("Following input directories do not exist: ",
            paste(input_table$input_dir[missing_dirs], collapse = ", "))
  }

  return(input_table)
}

# Function to format and clean results dataframe
format_results_dataframe <- function(grps_pivoted) {
  # Rename group_id_new to Satellite_family
  colnames(grps_pivoted)[colnames(grps_pivoted) == "group_id_new"] <- "Satellite_family"

  # Remove group_id column
  grps_pivoted$group_id <- NULL

  # Make Satellite_family the first column
  grps_pivoted <- grps_pivoted[, c("Satellite_family",
                                   setdiff(colnames(grps_pivoted), "Satellite_family"))]

  return(grps_pivoted)
}

# Function to generate annotation statistics
generate_annotation_report <- function(grps_pivoted) {
  # Calculate annotation counts per family
  annot_counts <- sapply(grps_pivoted$prevalent_annot,
                        function(x) length(unlist(strsplit(x, ", "))))

  # Calculate annotation distribution across families
  annot_in_family_counts <- table(grps_pivoted$prevalent_annot[grps_pivoted$prevalent_annot != ""])

  report <- list(
    no_annotation = sum(annot_counts == 0),
    single_annotation = sum(annot_counts == 1),
    conflicting_annotation = sum(annot_counts > 1),
    one_to_one_mapping = sum(annot_in_family_counts == 1),
    one_to_many_mapping = sum(annot_in_family_counts > 1),
    total_families = nrow(grps_pivoted)
  )

  # Print report
  cat("\n=== Annotation Statistics ===\n")
  cat(sprintf("Total satellite families: %d\n", report$total_families))
  cat(sprintf("Families with no similarity-based annotation: %d (%.1f%%)\n",
              report$no_annotation, 100 * report$no_annotation / report$total_families))
  cat(sprintf("Families with single annotation: %d (%.1f%%)\n",
              report$single_annotation, 100 * report$single_annotation / report$total_families))
  cat(sprintf("Families with conflicting annotations: %d (%.1f%%)\n",
              report$conflicting_annotation, 100 * report$conflicting_annotation / report$total_families))
  cat(sprintf("Annotations with one-to-one family mapping: %d\n", report$one_to_one_mapping))
  cat(sprintf("Annotations with one-to-many family mapping: %d\n", report$one_to_many_mapping))
  cat("===========================\n\n")

  return(report)
}

# Function to construct file paths for a given prefix
construct_file_paths <- function(tc_prefix) {
  list(
    annotation_gff = paste0(tc_prefix, "_annotation.gff3"),
    clustering_gff = paste0(tc_prefix, "_clustering.gff3"),
    annotation_tsv = paste0(tc_prefix, "_annotation.tsv"),
    ssrs_summary = paste0(tc_prefix, "_tarean/SSRS_summary.csv")
  )
}

# Function to process and export annotated GFF files
export_annotated_gff_files <- function(input_table, grps_pivoted, ssrs_groups, output_directory) {
  # Create output directory for GFF files
  gff_output_dir <- file.path(output_directory, "gff3")
  dir.create(gff_output_dir, recursive = TRUE, showWarnings = FALSE)

  # Create satellite family name mapping
  sf_names <- sprintf("SF_%04d", grps_pivoted$Satellite_family)

  for (i in seq_len(nrow(input_table))) {
    prefix <- input_table$sample_code[i]
    tc_prefix <- input_table$tidecluster_prefix[i]
    input_dir <- input_table$input_dir[i]

    # Get file paths
    paths <- construct_file_paths(tc_prefix)

    # Find available GFF file (prefer annotation over clustering)
    gff_path <- file.path(input_dir, paths$annotation_gff)
    if (!file.exists(gff_path)) {
      gff_path <- file.path(input_dir, paths$clustering_gff)
      if (!file.exists(gff_path)) {
        warning(sprintf("No GFF3 file found for prefix '%s', skipping", prefix))
        next
      }
      message(sprintf("Using clustering GFF3 file for prefix: %s", prefix))
    }

    # Create TRC to Satellite family mapping
    trc_to_sf <- create_trc_to_sf_mapping(grps_pivoted, prefix, sf_names)

    # Create SSRS type mapping
    ssrs_mapping <- create_ssrs_mapping(ssrs_groups, prefix)

    # Read, annotate, and export GFF
    gff <- rtracklayer::import(gff_path)

    # Add Satellite_family annotation
    gff$Satellite_family <- trc_to_sf[gff$Name]

    # Add SSRS type annotation
    gff$ssrs_type <- ssrs_mapping[gff$Name]

    # Export annotated GFF
    output_path <- file.path(gff_output_dir, paste0(prefix, "_tc_annotated.gff3"))
    rtracklayer::export(gff, output_path, format = "gff3")

    message(sprintf("Exported annotated GFF for %s to %s", prefix, output_path))
  }
}

# Helper function to create TRC to Satellite Family mapping
create_trc_to_sf_mapping <- function(grps_pivoted, prefix, sf_names) {
  sf_list <- setNames(grps_pivoted[[prefix]], sf_names)

  # Remove empty entries
  sf_list <- sf_list[!is.na(sf_list) & sf_list != ""]

  # Convert to TRC -> SF dictionary
  trc_list <- strsplit(sf_list, ", ")
  trc_to_sf <- setNames(rep(names(trc_list), lengths(trc_list)), unlist(trc_list))

  return(trc_to_sf)
}

# Helper function to create SSRS mapping
create_ssrs_mapping <- function(ssrs_groups, prefix) {
  ssrs_trc_col <- paste0(prefix, "_trc_id")
  ssrs_trc <- ssrs_groups[[ssrs_trc_col]]

  # Create mapping, replacing commas with underscores in pattern names
  ssrs_dict <- setNames(
    gsub(", ", "_", ssrs_groups$major_pattern[!is.na(ssrs_trc)]),
    ssrs_trc[!is.na(ssrs_trc)]
  )

  return(ssrs_dict)
}

# Function to export all results
export_results <- function(grps_pivoted, ssrs_groups, annotation_report, output_directory) {
  # Export satellite families table
  families_file <- file.path(output_directory, "trc_satellite_families.tsv")
  write.table(grps_pivoted, file = families_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
  message(sprintf("Exported satellite families to %s", families_file))

  # Export SSRS groups table
  ssrs_file <- file.path(output_directory, "ssrs_groups.tsv")
  write.table(ssrs_groups, file = ssrs_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
  message(sprintf("Exported SSRS groups to %s", ssrs_file))

  # Export annotation report
  report_file <- file.path(output_directory, "annotation_report.txt")
  report_lines <- c(
    "TideCluster Comparative Analysis - Annotation Report",
    "====================================================",
    sprintf("Analysis date: %s", Sys.Date()),
    sprintf("Total satellite families: %d", annotation_report$total_families),
    sprintf("Families with no annotation: %d", annotation_report$no_annotation),
    sprintf("Families with single annotation: %d", annotation_report$single_annotation),
    sprintf("Families with conflicting annotations: %d", annotation_report$conflicting_annotation),
    sprintf("One-to-one annotation mappings: %d", annotation_report$one_to_one_mapping),
    sprintf("One-to-many annotation mappings: %d", annotation_report$one_to_many_mapping)
  )
  writeLines(report_lines, report_file)
  message(sprintf("Exported annotation report to %s", report_file))
}

# Main execution function
main <- function(opt) {
  # Load required libraries
  suppressPackageStartupMessages({
    library(Biostrings)
    library(igraph)
    library(rtracklayer)
    library(reshape2)
  })

  # Validate and read input
  message("Reading and validating input file...")
  input_table <- validate_and_read_input(opt$input)

  # Extract configuration
  prefix <- input_table$sample_code
  tc_code <- input_table$tidecluster_prefix
  input_tc_dirs <- input_table$input_dir

  # Run main analysis
  message("\nStarting TRC comparative analysis...")
  result <- process_trc_analysis(
    input_tc_dirs = input_tc_dirs,
    prefix = prefix,
    tc_code = tc_code,
    mmseqs2_path = opt$mmseqs2_path,
    output_directory = opt$output_directory,
    ncpu = opt$cpu
  )

  # Format results
  message("\nFormatting results...")
  grps_pivoted <- format_results_dataframe(result$groups$df)

  # Generate annotation report
  annotation_report <- generate_annotation_report(grps_pivoted)

  # Export annotated GFF files
  message("\nExporting annotated GFF files...")
  export_annotated_gff_files(input_table, grps_pivoted, result$ssrs_groups,
                            opt$output_directory)

  # Export all results
  message("\nExporting result tables...")
  export_results(grps_pivoted, result$ssrs_groups, annotation_report,
                opt$output_directory)

  message("\nAnalysis complete! Results saved to: ", opt$output_directory)
}

# Parse command line arguments
library(optparse)

option_list <- list(
  make_option(
    c("-i", "--input"), type="character", default=NULL,
    help=("Input tab delimited file specifying input TideCluster data, the file should have three columns: 'input_dir', 'sample_code' and 'tidecluster_prefix'.")),
  make_option(
    c("-c", "--cpu"), type="numeric", default=5,
    help="Number of CPU threads to use for MMseqs2 search [default: %default]"),
  make_option(
    c("-m", "--mmseqs2_path"), type="character", default='mmseqs',
    help="Path to MMseqs2 executable [default: %default]"),
  make_option(
    c("-o", "--output_directory"), type="character", default="tc_comparative_analysis",
    help="Output directory for results [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is mandatory. Use -i or --input to specify it.")
}

# Run main analysis
main(opt)


