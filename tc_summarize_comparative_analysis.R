#!/usr/bin/env Rscript

# TC Comparative Analysis Summary and Visualization Script
# Generates HTML report with interactive visualizations of TRC comparative analysis results

library(optparse)
library(jsonlite)

# Function to parse satellite families data
parse_satellite_families <- function(input_file) {
  cat("Reading satellite families data from:", input_file, "\n")
  
  # Read the TSV file
  data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     check.names = FALSE, na.strings = "", quote = "")
  
  # Identify sample columns (those without suffixes)
  all_cols <- colnames(data)
  sample_cols <- all_cols[!grepl("_length$|_annot$|_prevalent_annot$|^Satellite_family$|^prevalent_annot$|^group_id$", all_cols)]
  samples <- sample_cols
  
  cat("Detected samples:", paste(samples, collapse = ", "), "\n")
  cat("Total families:", nrow(data), "\n")
  
  return(list(
    data = data,
    samples = samples,
    n_families = nrow(data)
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
prepare_detailed_family_data <- function(parsed_data) {
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
      
      family_row[[paste0(sample, "_trc_count")]] <- trc_count
      family_row[[paste0(sample, "_length")]] <- total_length
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

# Function to generate HTML report
generate_html_report <- function(parsed_data, overview_stats, shared_matrix, 
                                detailed_data, clustering_result, output_file) {
  
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
    sample_order = clustering_result$sample_order
  )
  
  # Convert to JSON
  json_data <- toJSON(js_data, auto_unbox = FALSE, pretty = TRUE)
  
  # Create HTML content
  html_content <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>TRC Comparative Analysis Report</title>
    <link rel="stylesheet" href="js/styles.css">
</head>
<body>
    <div class="container">
        <header>
            <h1>TRC Comparative Analysis Report</h1>
            <p class="subtitle">Interactive visualization of tandem repeat cluster families across samples</p>
            <p class="timestamp">Generated on: %s</p>
        </header>
        
        <nav class="tab-nav">
            <button class="tab-button active" onclick="showTab(\'overview\')">Overview</button>
            <button class="tab-button" onclick="showTab(\'shared-families\')">Shared Families</button>
            <button class="tab-button" onclick="showTab(\'detailed-families\')">Detailed Families</button>
        </nav>
        
        <main>
            <div id="overview" class="tab-content active">
                <h2>Sample Overview</h2>
                <div id="overview-table-container"></div>
            </div>
            
            <div id="shared-families" class="tab-content">
                <h2>Family Presence Patterns - Shared Families Matrix</h2>
                <div id="shared-matrix-container"></div>
            </div>
            
            <div id="detailed-families" class="tab-content">
                <h2>Family Presence Patterns - Detailed View</h2>
                <div class="controls">
                    <label>
                        <input type="radio" name="view-mode" value="trc-count" checked> TRC Count
                    </label>
                    <label>
                        <input type="radio" name="view-mode" value="length"> Total Length
                    </label>
                </div>
                <div id="detailed-matrix-container"></div>
            </div>
        </main>
    </div>
    
    <script>
        // Embed data
        const data = %s;
    </script>
    <script src="js/utils.js"></script>
    <script src="js/overview-table.js"></script>
    <script src="js/shared-matrix.js"></script>
    <script src="js/detailed-matrix.js"></script>
    <script src="js/main.js"></script>
</body>
</html>', 
    Sys.time(),
    json_data
  )
  
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
  
  # Create main.js
  main_js <- '// Main JavaScript file for TRC Comparative Analysis Report

// Tab switching functionality
function showTab(tabName) {
    // Hide all tab contents
    const tabContents = document.querySelectorAll(".tab-content");
    tabContents.forEach(content => {
        content.classList.remove("active");
    });
    
    // Remove active class from all tab buttons
    const tabButtons = document.querySelectorAll(".tab-button");
    tabButtons.forEach(button => {
        button.classList.remove("active");
    });
    
    // Show selected tab content
    document.getElementById(tabName).classList.add("active");
    
    // Add active class to clicked button
    event.target.classList.add("active");
    
    // Initialize the appropriate visualization
    switch(tabName) {
        case "overview":
            initOverviewTable();
            break;
        case "shared-families":
            initSharedMatrix();
            break;
        case "detailed-families":
            initDetailedMatrix();
            break;
    }
}

// Initialize the report
document.addEventListener("DOMContentLoaded", function() {
    initOverviewTable();
    
    // Add event listeners for view mode radio buttons
    const viewModeInputs = document.querySelectorAll(\'input[name="view-mode"]\');
    viewModeInputs.forEach(input => {
        input.addEventListener("change", function() {
            if (document.getElementById("detailed-families").classList.contains("active")) {
                initDetailedMatrix();
            }
        });
    });
});'
  
  writeLines(main_js, file.path(js_dir, "main.js"))
  
  # Create utils.js
  utils_js <- '// Utility functions for TRC Comparative Analysis Report

// Format numbers with appropriate units
function formatLength(length) {
    if (length >= 1000000) {
        return (length / 1000000).toFixed(2) + " Mbp";
    } else if (length >= 1000) {
        return (length / 1000).toFixed(2) + " kbp";
    } else {
        return length + " bp";
    }
}

// Create color scale for heatmaps
function getColorScale(values, colorScheme = "Blues") {
    const minVal = Math.min(...values);
    const maxVal = Math.max(...values);
    
    if (minVal === maxVal) {
        return (val) => "#f0f0f0";
    }
    
    return function(val) {
        const intensity = (val - minVal) / (maxVal - minVal);
        
        if (colorScheme === "Blues") {
            const blue = Math.floor(255 - intensity * 200);
            return `rgb(${blue}, ${blue}, 255)`;
        } else if (colorScheme === "Reds") {
            const red = Math.floor(255 - intensity * 200);
            return `rgb(255, ${red}, ${red})`;
        } else {
            // Default grayscale
            const gray = Math.floor(255 - intensity * 200);
            return `rgb(${gray}, ${gray}, ${gray})`;
        }
    };
}

// Create tooltip
function createTooltip() {
    const tooltip = document.createElement("div");
    tooltip.id = "tooltip";
    tooltip.style.cssText = `
        position: absolute;
        background: rgba(0, 0, 0, 0.8);
        color: white;
        padding: 8px;
        border-radius: 4px;
        font-size: 12px;
        pointer-events: none;
        z-index: 1000;
        display: none;
    `;
    document.body.appendChild(tooltip);
    return tooltip;
}

// Show tooltip
function showTooltip(event, text) {
    const tooltip = document.getElementById("tooltip") || createTooltip();
    tooltip.innerHTML = text;
    tooltip.style.display = "block";
    tooltip.style.left = (event.pageX + 10) + "px";
    tooltip.style.top = (event.pageY - 10) + "px";
}

// Hide tooltip
function hideTooltip() {
    const tooltip = document.getElementById("tooltip");
    if (tooltip) {
        tooltip.style.display = "none";
    }
}'
  
  writeLines(utils_js, file.path(js_dir, "utils.js"))
  
  # Create overview-table.js
  overview_js <- '// Overview table functionality

function initOverviewTable() {
    const container = document.getElementById("overview-table-container");
    const stats = data.overview_stats;
    
    let html = `
        <table class="overview-table">
            <thead>
                <tr>
                    <th>Sample Name</th>
                    <th>Number of TRCs</th>
                    <th>Number of Families</th>
                    <th>Number of Unique Families</th>
                    <th>Total Length</th>
                </tr>
            </thead>
            <tbody>
    `;
    
    stats.forEach(row => {
        html += `
            <tr>
                <td class="sample-name">${row.Sample}</td>
                <td class="numeric">${row.Number_of_TRCs.toLocaleString()}</td>
                <td class="numeric">${row.Number_of_Families.toLocaleString()}</td>
                <td class="numeric">${row.Number_of_Unique_Families.toLocaleString()}</td>
                <td class="numeric">${row.Total_Length}</td>
            </tr>
        `;
    });
    
    html += `
            </tbody>
        </table>
    `;
    
    container.innerHTML = html;
}'
  
  writeLines(overview_js, file.path(js_dir, "overview-table.js"))
  
  # Create shared-matrix.js
  shared_js <- '// Shared families matrix functionality

function initSharedMatrix() {
    const container = document.getElementById("shared-matrix-container");
    const matrix = data.shared_matrix;
    const samples = data.samples;
    
    // Get all values for color scaling
    const allValues = [];
    for (let i = 0; i < samples.length; i++) {
        for (let j = 0; j < samples.length; j++) {
            allValues.push(matrix[i][j]);
        }
    }
    
    const colorScale = getColorScale(allValues, "Blues");
    
    let html = `
        <div class="matrix-container">
            <table class="shared-matrix">
                <thead>
                    <tr>
                        <th class="corner-cell"></th>
    `;
    
    // Column headers
    samples.forEach(sample => {
        html += `<th class="sample-header">${sample}</th>`;
    });
    
    html += `
                    </tr>
                </thead>
                <tbody>
    `;
    
    // Matrix rows
    for (let i = 0; i < samples.length; i++) {
        html += `<tr>`;
        html += `<th class="sample-header">${samples[i]}</th>`;
        
        for (let j = 0; j < samples.length; j++) {
            const value = matrix[i][j];
            const bgColor = colorScale(value);
            const isDiagonal = i === j;
            
            html += `
                <td class="matrix-cell ${isDiagonal ? \'diagonal\' : \'\'}" 
                    style="background-color: ${bgColor}"
                    onmouseover="showTooltip(event, \'${samples[i]} vs ${samples[j]}: ${value} shared families\')"
                    onmouseout="hideTooltip()">
                    ${value}
                </td>
            `;
        }
        
        html += `</tr>`;
    }
    
    html += `
                </tbody>
            </table>
        </div>
        <div class="matrix-legend">
            <p><strong>Interpretation:</strong></p>
            <ul>
                <li>Rows and columns are samples ordered by hierarchical clustering</li>
                <li>Diagonal cells show total families present in each sample</li>
                <li>Off-diagonal cells show shared families between sample pairs</li>
                <li>Darker colors indicate higher numbers of shared families</li>
            </ul>
        </div>
    `;
    
    container.innerHTML = html;
}'
  
  writeLines(shared_js, file.path(js_dir, "shared-matrix.js"))
  
  # Create detailed-matrix.js
  detailed_js <- '// Detailed families matrix functionality

function initDetailedMatrix() {
    const container = document.getElementById("detailed-matrix-container");
    const families = data.detailed_families;
    const samples = data.samples;
    
    // Get selected view mode
    const viewMode = document.querySelector(\'input[name="view-mode"]:checked\').value;
    const isLengthView = viewMode === "length";
    
    // Get all values for color scaling
    const allValues = [];
    families.forEach(family => {
        samples.forEach(sample => {
            const key = isLengthView ? `${sample}_length` : `${sample}_trc_count`;
            allValues.push(family[key] || 0);
        });
    });
    
    const colorScale = getColorScale(allValues.filter(v => v > 0), "Reds");
    
    let html = `
        <div class="matrix-container">
            <table class="detailed-matrix">
                <thead>
                    <tr>
                        <th class="corner-cell">Family ID</th>
    `;
    
    // Sample headers
    samples.forEach(sample => {
        html += `<th class="sample-header">${sample}</th>`;
    });
    
    html += `<th class="annotation-header">Prevalent Annotation</th>`;
    html += `</tr></thead><tbody>`;
    
    // Family rows
    families.forEach(family => {
        html += `<tr>`;
        html += `<td class="family-id">SF_${String(family.family_id).padStart(4, "0")}</td>`;
        
        samples.forEach(sample => {
            const trcKey = `${sample}_trc_count`;
            const lengthKey = `${sample}_length`;
            const trcCount = family[trcKey] || 0;
            const length = family[lengthKey] || 0;
            
            const displayValue = isLengthView ? formatLength(length) : trcCount;
            const colorValue = isLengthView ? length : trcCount;
            const bgColor = colorValue > 0 ? colorScale(colorValue) : "#f9f9f9";
            
            const tooltipText = `Family ${family.family_id} in ${sample}<br>` +
                              `TRCs: ${trcCount}<br>` +
                              `Length: ${formatLength(length)}`;
            
            html += `
                <td class="matrix-cell ${colorValue > 0 ? \'has-value\' : \'empty\'}" 
                    style="background-color: ${bgColor}"
                    onmouseover="showTooltip(event, \'${tooltipText}\')"
                    onmouseout="hideTooltip()">
                    ${colorValue > 0 ? displayValue : ""}
                </td>
            `;
        });
        
        // Prevalent annotation column
        const annotation = family.prevalent_annot || "";
        html += `<td class="annotation-cell">${annotation}</td>`;
        html += `</tr>`;
    });
    
    html += `
                </tbody>
            </table>
        </div>
        <div class="matrix-legend">
            <p><strong>Interpretation:</strong></p>
            <ul>
                <li>Rows are satellite families sorted by family index</li>
                <li>Columns are samples ordered by hierarchical clustering</li>
                <li>Cells show ${isLengthView ? "total genomic length" : "number of TRCs"} for each family-sample combination</li>
                <li>Empty cells indicate the family is absent in that sample</li>
                <li>Darker colors indicate ${isLengthView ? "longer total length" : "more TRCs"}</li>
            </ul>
        </div>
    `;
    
    container.innerHTML = html;
}'
  
  writeLines(detailed_js, file.path(js_dir, "detailed-matrix.js"))
  
  # Create styles.css
  styles_css <- '/* Styles for TRC Comparative Analysis Report */

* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    line-height: 1.6;
    color: #333;
    background-color: #f5f5f5;
}

.container {
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
    background-color: white;
    box-shadow: 0 0 10px rgba(0,0,0,0.1);
    min-height: 100vh;
}

header {
    text-align: center;
    margin-bottom: 30px;
    padding-bottom: 20px;
    border-bottom: 2px solid #e0e0e0;
}

header h1 {
    color: #2c3e50;
    font-size: 2.5em;
    margin-bottom: 10px;
}

header .subtitle {
    color: #7f8c8d;
    font-size: 1.2em;
    margin-bottom: 5px;
}

header .timestamp {
    color: #95a5a6;
    font-size: 0.9em;
}

.tab-nav {
    display: flex;
    justify-content: center;
    margin-bottom: 30px;
    border-bottom: 1px solid #ddd;
}

.tab-button {
    background: none;
    border: none;
    padding: 12px 24px;
    cursor: pointer;
    font-size: 16px;
    color: #666;
    border-bottom: 3px solid transparent;
    transition: all 0.3s ease;
}

.tab-button:hover {
    color: #2c3e50;
    background-color: #f8f9fa;
}

.tab-button.active {
    color: #2c3e50;
    border-bottom-color: #3498db;
    background-color: #f8f9fa;
}

.tab-content {
    display: none;
}

.tab-content.active {
    display: block;
}

.tab-content h2 {
    color: #2c3e50;
    margin-bottom: 20px;
    font-size: 1.8em;
}

/* Overview Table Styles */
.overview-table {
    width: 100%;
    border-collapse: collapse;
    margin-bottom: 20px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

.overview-table th,
.overview-table td {
    padding: 12px 15px;
    text-align: center;
    border-bottom: 1px solid #ddd;
}

.overview-table th {
    background-color: #3498db;
    color: white;
    font-weight: 600;
    text-transform: uppercase;
    font-size: 0.9em;
}

.overview-table tbody tr:hover {
    background-color: #f5f5f5;
}

.overview-table .sample-name {
    font-weight: 600;
    color: #2c3e50;
}

/* Matrix Styles */
.matrix-container {
    overflow-x: auto;
    margin-bottom: 20px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

.shared-matrix,
.detailed-matrix {
    border-collapse: collapse;
    font-size: 12px;
    background-color: white;
}

.shared-matrix th,
.shared-matrix td,
.detailed-matrix th,
.detailed-matrix td {
    border: 1px solid #ddd;
    text-align: center;
    min-width: 80px;
    position: relative;
}

.corner-cell {
    background-color: #ecf0f1 !important;
}

.sample-header {
    background-color: #3498db !important;
    color: white !important;
    font-weight: 600;
    padding: 8px 4px;
    writing-mode: vertical-rl;
    text-orientation: mixed;
    min-height: 100px;
}

.annotation-header {
    background-color: #e74c3c !important;
    color: white !important;
    font-weight: 600;
    padding: 8px;
    min-width: 200px;
}

.matrix-cell {
    padding: 6px 4px;
    cursor: pointer;
    transition: all 0.2s ease;
    font-weight: 500;
}

.matrix-cell:hover {
    border: 2px solid #2c3e50;
    z-index: 10;
}

.matrix-cell.diagonal {
    font-weight: bold;
    border: 2px solid #2c3e50;
}

.matrix-cell.empty {
    color: #ccc;
}

.matrix-cell.has-value {
    color: #2c3e50;
}

.family-id {
    background-color: #f8f9fa !important;
    font-weight: 600;
    color: #2c3e50;
    padding: 6px 8px;
    text-align: center;
    min-width: 80px;
}

.annotation-cell {
    padding: 6px 8px;
    text-align: left;
    font-size: 11px;
    max-width: 200px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

/* Controls */
.controls {
    margin-bottom: 20px;
    padding: 15px;
    background-color: #f8f9fa;
    border-radius: 6px;
    border: 1px solid #e9ecef;
}

.controls label {
    margin-right: 20px;
    cursor: pointer;
    font-weight: 500;
}

.controls input[type="radio"] {
    margin-right: 5px;
}

/* Legend */
.matrix-legend {
    background-color: #f8f9fa;
    padding: 15px;
    border-radius: 6px;
    border-left: 4px solid #3498db;
}

.matrix-legend p {
    font-weight: 600;
    margin-bottom: 10px;
    color: #2c3e50;
}

.matrix-legend ul {
    list-style-type: none;
    padding-left: 0;
}

.matrix-legend li {
    margin-bottom: 5px;
    padding-left: 20px;
    position: relative;
}

.matrix-legend li:before {
    content: "•";
    color: #3498db;
    font-weight: bold;
    position: absolute;
    left: 0;
}

/* Responsive design */
@media (max-width: 768px) {
    .container {
        padding: 10px;
    }
    
    header h1 {
        font-size: 2em;
    }
    
    .tab-button {
        padding: 10px 16px;
        font-size: 14px;
    }
    
    .sample-header {
        writing-mode: horizontal-tb;
        text-orientation: initial;
        min-height: auto;
        padding: 8px;
    }
    
    .matrix-container {
        font-size: 10px;
    }
}'
  
  writeLines(styles_css, file.path(js_dir, "styles.css"))
  
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
  
  # Prepare detailed family data
  cat("Preparing detailed family data...\n")
  detailed_data <- prepare_detailed_family_data(parsed_data)
  
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
                      detailed_data, clustering_result, opt$output)
  
  cat("Report generation completed successfully!\n")
  cat("Open", opt$output, "in your web browser to view the interactive report.\n")
}

# Parse command line arguments
option_list <- list(
  make_option(
    c("-i", "--input"), type="character", default=NULL,
    help="Input TSV file from tc_comparative_analysis.R script (required)"),
  make_option(
    c("-o", "--output"), type="character", default="trc_comparative_report.html",
    help="Output HTML file path [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is mandatory. Use -i or --input to specify it.")
}

if (!file.exists(opt$input)) {
  stop("Input file does not exist: ", opt$input)
}

# Run main analysis
main(opt)