#!/user/bin/env Rscript

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(hwriter))
suppressPackageStartupMessages(library(R2HTML))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))

kmers2graph <- function(kmers, mode = "strong", prop = NULL) {
  kmerLength <- nchar(kmers[1, 1])
  if (ncol(kmers) == 2) {
    kmers$size <- kmers[, 2]/sum(kmers[, 2])
  }
  colnames(kmers) <- c("name", "count", "size")
  if (!is.null(prop)) {  # tohle se nepouziva(prop je null), a je to asi spatne - filtuje se to pred tridenim!!
    p <- cumsum(kmers$size)
    kmers <- kmers[p < prop, ]
  }
  kmers <- kmers[order(kmers$size), ]
  ## convert kmers to fasta file
  kms <- data.frame(kmer = substring(kmers$name, 1, kmerLength - 1), ids = 1:nrow(kmers), stringsAsFactors = FALSE)
  kme <- data.frame(kmer = substring(kmers$name, 2), ide = 1:nrow(kmers), stringsAsFactors = FALSE)

  ## df = merge(kms,kme, by = 'kmer',all=FALSE)[,c(2,3,1)]
  df <- inner_join(kme, kms, by = 'kmer', multiple='all')[, c(2, 3)]
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
  }
  
  whg <- apply(cbind(kmers[df[, 1], 2], V2 = kmers[df[, 2], 2]), 1, gm_mean)
  G <- graph.data.frame(data.frame(V1 = kmers$name[df[, 1]], V2 = kmers$name[df[,
                                                                               2]], weight = whg), vertices = kmers[, 1:3])
                                        # separate to connected components:
  ccs <- clusters(G, mode = mode)$membership
  sel_cls <- which(tabulate(ccs) > 1)
  Gs <- list()
  for (i in seq_along(sel_cls)) {
    Gs[[i]] <- induced.subgraph(G, vids = which(ccs %in% sel_cls[i]))
  }
  ## reorder!!!
  Gs <- Gs[order(sapply(Gs, vcount), decreasing = TRUE)]
  return(Gs)
}


OGDFlayout <- function(G, ncol = NULL, alg = "fmmm", OGDF = getOption("OGDF")) {
  ## is ogdf binary available?
  if (is.null(OGDF)) {
    OGDF <- Sys.getenv("OGDF")
    if ("" == OGDF) {
      options(warn = -1)
      OGDF <- system("which runOGDFlayout", intern = TRUE)
      options(warn = 0)
      if (length(OGDF) == 0) {
        cat("path to runOGDFlayout not found\n")
        return(NULL)
      }
      
    }
  }
  if (is.null(ncol)) {
    if (is.null(E(G)$weight)) {
      el <- data.frame(get.edgelist(G, names = TRUE), rep(1, ecount(G)))
    } else {
      el <- data.frame(get.edgelist(G, names = TRUE), E(G)$weight)
    }
    ncol <- paste(tempfile(pattern = as.character(Sys.getpid())), ".layout", sep = "")
    write.table(el, file = ncol, row.names = FALSE, col.names = FALSE, sep = "\t", 
                quote = FALSE)
  } else {
                                        # copy ncol:
    ncol_tmp <- paste(tempfile(pattern = as.character(Sys.getpid())), ".layout",
                      sep = "")
    file.copy(ncol, ncol_tmp)
    ncol <- ncol_tmp
  }
  algopt <- c("fmmm", "sm", "fme")
  if (!(alg %in% c("fmmm", "sm", "fme") && TRUE)) {
    stop("alg must by :", algopt, "\n")
    
  }
  
                                        # output file:
  Louts <- list()
  layout_file <- tempfile(pattern = as.character(Sys.getpid()))
  for (i in alg) {
    cmd <- paste(OGDF, "-alg", i, "-iff layout -off layout -out", layout_file,
                 ncol)
    system(cmd, intern = TRUE)
    L <- read.table(layout_file, skip = ecount(G))
    L <- L[match(V(G)$name, L$V2), ]
    Lout <- as.matrix(L[, -(1:2)])
    unlink(layout_file)
    Louts[[i]] <- Lout
    
  }
                                        # clean up
  unlink(ncol)
  return(Louts)
}

xcolor_code <- c(A = "#FF0000", C = "#00FF00", T = "#0000FF", G = "#FFFF00")

kmers2color <- function(s, position = NULL) {
  if (is.null(position)) {
    position <- round(nchar(s[1])/2, 0)
    ## position = 1
    position <- nchar(s[1])
  }
  color_code <- c(A = "#FF0000", C = "#00FF00", T = "#0000FF", G = "#FFFF00")
  color_base <- substring(s, position, position)
  colors <- color_code[color_base]
  names(colors) <- color_base
  return(colors)
}
get_sequence <- function(g, v, position = NULL) {
  s <- V(g)$name
  if (is.null(position)) {
    position <- round(nchar(s[1])/2, 0)
    ## position = 1
    position <- nchar(s[1])
  }
  nt <- paste(substring(s[v], position, position), collapse = "")
  return(nt)
}

get_mimimal_cc <- function(km, thr = 20, min_coverage = 0.45, step = 2, start = NULL) {
  if (is.null(start)) {
    i <- sum(cumsum(km$freq) < 0.5)
  } else {
    i <- sum(cumsum(km$freq) < start)
  }
  continue <- TRUE
  while (continue) {
    if (i > nrow(km)) {
      i <- nrow(km)
      continue <- FALSE
      step <- 1
    }
    GG <- kmers2graph(km[1:i, ])
    if (length(GG) > 0) {
      if (vcount(GG[[1]]) > thr) {
        if (sum(V(GG[[1]])$size) >= min_coverage) {
          GG[[1]]$input_coverage <- sum(km$freq[1:i])
          GG[[1]]$L <- OGDFlayout(GG[[1]])[[1]]
          return(GG[[1]])
        }
      }
    }
    i <- round(i * step)
    
  }
  if (length(GG) == 0 | is.null(GG)) {
    return(NULL)
  }
  
  GG[[1]]$input_coverage <- sum(km$freq[1:i])
  GG[[1]]$L <- OGDFlayout(GG[[1]])[[1]]
  return(GG[[1]])
}

paths2string <- function(paths) {
  pathstring <- sapply(lapply(lapply(paths, as_ids), substring, 1, 1), paste, collapse = "")
  return(pathstring)
}

align_paths <- function(paths, G) {
  shift <- rep(NA, length(paths))
  thr <- 0  # minimal length
  tr_paths <- list()
  Seqs <- list()
  centre_node <- as.numeric(names(sort(table(unlist(paths)), decreasing = TRUE)))[[1]]
  
  for (i in seq_along(paths)) {
    if (centre_node %in% paths[[i]]) {
      S <- which(paths[[i]] %in% centre_node)
      shift[i] <- S
      if (S == 1) {
        tr_paths[[i]] <- paths[[i]]
      } else {
        tr_paths[[i]] <- c(paths[[i]][S:length(paths[[i]])], paths[[i]][1:(S -
                                                                          1)])
      }
      Seqs[[i]] <- get_sequence(G, tr_paths[[i]])
    } else {
      shift[i] <- NA
    }
  }
  paths_n <- lapply(paths, as.numeric)
  tr_paths_n <- do.call(cbind, lapply(tr_paths, as.numeric))
  new_shift <- shift
  for (i in which(is.na(shift))) {
    score <- numeric(length(paths_n))
    for (S in seq_along(paths_n[[i]])) {
      if (S == 1) {
        path_tmp_n <- paths_n[[i]]
      } else {
        path_tmp_n <- c(paths_n[[i]][S:length(paths_n[[i]])], paths_n[[i]][1:(S -
                                                                             1)])
      }
      score[S] <- sum(tr_paths_n == path_tmp_n)
    }
    if (sum(score) != 0) {
      S <- which.max(score)
      new_shift[i] <- S
      if (S == 1) {
        tr_paths[[i]] <- paths[[i]]
      } else {
        tr_paths[[i]] <- c(paths[[i]][S:length(paths[[i]])], paths[[i]][1:(S -
                                                                          1)])
      }
      Seqs[[i]] <- get_sequence(G, tr_paths[[i]])
    }
  }
  shift <- new_shift
                                        # try to shift based on the sequences itself
  return(list(Seqs = Seqs[!is.na(shift)], tr_paths = tr_paths[!is.na(shift)], shift = shift[!is.na(shift)]))
}

make_consensus <- function(paths_info, G) {
  include <- !is.na(paths_info$shift)
  ## make alignments using mafft
  aln <- mafft(unlist(paths_info$Seqs[include]))
  CM <- calculate_consensus_matrix(aln = aln, tr_paths = paths_info$tr_paths[include],
                                   G = G)
  gaps <- get_gaps_from_alignment(aln)
  CMnorm <- CM/rowSums(CM)
  bases <- colnames(CM)
  consensus <- sapply(apply(CMnorm, 1, function(x) bases[which(x > 0.2)][order(x[x >
                                                                                0.2], decreasing = TRUE)]), paste, collapse = "")
  consensus2 <- gsub("-", "", paste0(substring(consensus, 1, 1), collapse = ""),
                     fixed = TRUE)
  number_of_SNP <- sum(rowSums(CM > 0) > 1)
  SNP_positions <- which(rowSums(CM > 0) > 1)
  number_of_position_with_indels <- sum(colSums(do.call(rbind, strsplit(as.character(aln),
                                                                        "")) == "-") > 0)
  indel_positions <- which(colSums(do.call(rbind, strsplit(as.character(aln), "")) ==
                                  "-") > 0)
  if (length(SNP_positions) > 0) {
    variable_sites <- unlist(c(c(mapply(paste, strsplit((sapply(apply(CMnorm,
                                                                      1, function(x) bases[which(x > 0.2)]), paste, collapse = ""))[SNP_positions],
                                                        ""), SNP_positions, sep = "_")), paste("-", indel_positions, sep = "_")))
  } else {
    variable_sites <- NULL
  }
  variable_positions <- unique(SNP_positions, indel_positions)
  return(list(aln = aln, CM = CM, CMnorm = CMnorm, consensus = consensus, consensus2 = consensus2, 
              number_of_SNP = number_of_SNP, SNP_positions = SNP_positions, number_of_position_with_indels = number_of_position_with_indels, 
              indel_positions = indel_positions, variable_positions = variable_positions, 
              variable_sites = variable_sites, gaps = gaps))
}

estimate_monomer <- function(G, weights = NULL, limit = NULL) {
  if (is.null(G)) {
    return(NULL)
  }
  ## estimate monomer from kmer based graph
  V(G)$id <- 1:vcount(G)
  GS <- induced_subgraph(G, vids = which(degree(G) == 2))  ## keep only vertices without branching
  cls <- clusters(GS)$membership
  
  
  ids <- mapply(FUN = function(x, y) x[which.max(y)], split(V(GS)$id, cls), split(V(GS)$size,
                                                                                  cls))  ## from each branch use only one vertex with larges size!
  
  
  ids <- ids[order(V(G)$size[ids], decreasing = TRUE)]
  ids_size <- V(G)$size[ids]
  N50 <- sum(cumsum(ids_size)/sum(ids_size) < 0.5)
  if (length(ids) > 10000) {
    ids <- ids[1:N50]
  }
  ## use only large vertices in search!  how many?
  el <- get.edgelist(G, names = FALSE)
  node_use <- numeric(vcount(G))
  LL <- numeric()
  ## W=numeric()
  i <- 0
  paths <- list()
  if (is.null(weights)) {
    weights <- (max(E(G)$weight) - (E(G)$weight) + 1)
    weights <- E(G)$weight^(-3)
  }
  included <- rep(FALSE, vcount(G))
  W_total <- sum(V(G)$size)
  
  coverage <- numeric()
  t0 <- c()
  i <- 0
  j <- 0
  for (n in ids) {
    j <- j + 1
    t0[j] <- Sys.time()
    if (included[n]) {
      next
    }
    m <- which(el[, 1] %in% n)
    i <- i + 1
    s <- get.shortest.paths(G, el[m, 2], el[m, 1], weights = weights, output = "vpath")
    included[as.numeric(s$vpath[[1]])] <- TRUE
    paths[[i]] <- s$vpath[[1]]
    LL[i] <- (length(s$vpath[[1]]))
  }
  
  ## evaluate if paths should be divided to variants - by length and path weight
  paths_clusters <- split(paths, LL)
  paths_clusters_tr <- mclapply(paths_clusters,
                                FUN = align_paths, G = G, mc.cores = getOption("CPU"),
                                mc.preschedule = FALSE
  )
  
  ## paths_clusters_tr = lapply(paths_clusters, FUN = align_paths, G = G) consensus
  paths_consensus <- mclapply(paths_clusters_tr, make_consensus,
                              G = G, mc.cores = getOption("CPU"),
                              mc.preschedule = FALSE
  )
  
  ## evaluate weight for individual paths:
  for (v in seq_along(paths_consensus)) {
    p <- paths_clusters_tr[[v]]$tr_paths
    ## clean
    p <- p[!sapply(p, function(x) anyNA(x) | is.null(x))]
    L <- sapply(p, length)
    p_groups <- split(p, L)
    w_groups <- sapply(p_groups,
                       function(x) sum(V(G)$size[unique(c(sapply(x,as.numeric)))]))
    total_score <- sum(V(G)$size[unique(c(unlist(sapply(p, as.numeric))))])
    LW <- data.frame(`Length estimate` = unique(L), weight = w_groups, stringsAsFactors = FALSE)
    LW <- LW[order(LW$weight, decreasing = TRUE), ]
    rownames(LW) <- NULL
    paths_consensus[[v]]$total_score <- total_score
    paths_consensus[[v]]$length_variant_score <- LW
    
  }
  return(list(estimates = paths_consensus, paths = paths_clusters_tr))
}

subject_coverage <- function(blastdf) {
  ## calculate coverage for all blast subjects
  coverage_profiles <- by(blastdf, INDICES = blastdf$sseqid, FUN = function(x) {
    as.numeric(coverage(IRanges(start = x$sstart, end = x$send)))
  })
  return(coverage_profiles)
  
}

## in alingment calculate significance from kmers
calculate_consensus_matrix <- function(aln, tr_paths, G) {
  bases <- c("A", "C", "G", "T", "-")
  positions <- lapply(strsplit(as.character(aln), split = ""), function(x) which(!x %in%
                                                                                "-"))
  base_matrix <- do.call(rbind, strsplit(as.character(aln), split = ""))
  weights <- matrix(0, nrow = length(aln), ncol = nchar(aln))
  kmer_rel_proportions <- (V(G)$size/table(factor(unlist(tr_paths), levels = 1:vcount(G))))
  for (i in seq_along(positions)) {
    weights[i, positions[[i]]] <- kmer_rel_proportions[tr_paths[[i]]]
  }
  names(kmer_rel_proportions) <- V(G)$name
  ## get weights for gaps by approximation
  fitgaps <- function(y) {
    if (sum(y == 0) == 0) {
      return(y)
    } else {
      y0 <- rep(y[y != 0], 3)
      x0 <- which(rep(y, 3) != 0)
      fitted <- approx(x0, y0, xout = seq_along(rep(y, 3)), rule = 1)
    }
    return(fitted$y[-seq_along(y)][seq_along(y)])
  }
  weights_with_gaps <- t(apply(weights, 1, fitgaps))

  ## get consensus matrix
  CM <- sapply(1:nchar(aln[[1]]), function(i) {
    sapply(split(weights_with_gaps[, i], factor(base_matrix[, i], levels = bases)), 
           sum)
  })
  return(t(CM)[, 1:5])
}


get_gaps_from_alignment <- function(aln) {
  as.character(aln)
  gaps_positions <- unique(do.call(rbind, str_locate_all(as.character(aln), "-+")))
  return(gaps_positions)
}

plot_kmer_graph <- function(G, L = NULL, vertex.size = NULL, ord = NULL, upto = NULL,
                            highlight = NULL) {
  if (!is.null(G$L) & is.null(L)) {
    L <- G$L
  }
  if (is.null(L)) {
    ## L=layout.kamada.kawai(G)
    L <- OGDFlayout(G)[[1]]
  }
  clr <- kmers2color(V(G)$name)
  if (!is.null(highlight)) {
    clr[highlight] <- "#00000080"
  }
  if (!is.null(ord)) {
    clr[ord[1:upto]] <- paste(clr[ord[1:upto]], "30", sep = "")
  }
  if (is.null(vertex.size)) {
    vertex.size <- rescale(V(G)$size, to = c(0.5, 6))
    
  }
  plot(G, layout = L, vertex.label = "", vertex.size = vertex.size, edge.curved = FALSE, 
       vertex.color = clr, vertex.frame.color = "#00000020", edge.width = 2, edge.arrow.mode = 1, 
       edge.arrow.size = 0.2)
  
}

mafft <- function(seqs, params = "--auto --thread 1 ") {
  if (length(seqs) < 2) {
    return(seqs)
  }
  infile <- tempfile()
  if (class(seqs) == "character") {
    seqs <- DNAStringSet(seqs)
  }
  writeXStringSet(seqs, file = infile)
  outfile <- tempfile()
  cmd <- paste("mafft --quiet --nuc ", params, infile, "2> /dev/null > ", outfile)
  system(cmd, intern = TRUE, ignore.stderr = FALSE)
  aln <- readDNAStringSet(outfile)
  unlink(c(outfile, infile))
  return(aln)
}

estimate_sample_size <- function(NV, NE, maxv, maxe) {
  ## density
  d <- (2 * NE)/(NV * (NV - 1))
  eEst <- (maxv * (maxv - 1) * d)/2
  nEst <- (d + sqrt(d^2 + 8 * d * maxe))/(2 * d)
  if (eEst >= maxe) {
    N <- round(nEst)
    E <- round((N * (N - 1) * d)/2)
    
  }
  if (nEst >= maxv) {
    N <- maxv
    E <- round((N * (N - 1) * d)/2)
    
  }
  return(N)
}


list2dictionary <- function(l){
    dict <- "{"
    q <- '"'
    for (i in 1:length(l)){
        if (class(l[[i]])=="character" | is.null(l[[i]])){
            q2 <- "'''"
        }else{
            q2 <- ''
        }
        dict <- paste0(
            dict,
            q,names(l)[i],q,":",q2, l[[i]], q2,", "
        )
    }
    dict <- paste0(dict, "}")
    return(dict)
}

wrap <- Vectorize(function(s, width = 80) {
  i1 <- seq(1, nchar(s), width)
  i2 <- seq(width, by = width, length.out = length(i1))
  return(paste(substring(s, i1, i2), collapse = "\n"))
})


tarean <- function(input_sequences, output_dir, min_kmer_length = 11, max_kmer_length = 27,
                   CPU = 2, sample_size = 10000,
                   include_layout = TRUE, paired = TRUE, lock_file=NULL) {
  options(nwarnings = 100000, warn = 1)
  options(CPU = CPU)
  time0 <- Sys.time()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  input_sequences_copy <- paste(output_dir, "/", basename(input_sequences), sep = "")
  
  if (!file.copy(input_sequences, input_sequences_copy, overwrite = TRUE)) {
    cat(paste("cannot copy", input_sequences, " to", output_dir), "\n")
    stop()
  }
  
  # lock_file <- waitForRAM(lock_file = lock_file)
  pair_completeness <- NULL
  ## sampling
  if (sample_size != 0) {
    s <- readDNAStringSet(input_sequences_copy)
    N <- length(s)
    ## pair completness must be calculated before sampling!
    if (N > sample_size) {
      set.seed(N)
      writeXStringSet(sample(s, sample_size), filepath = input_sequences_copy)
      if (paired) {
        pair_counts <- tabulate(table(gsub(".$", "", names(s))))
        pair_completeness <- 1 - pair_counts[1]/sum(pair_counts)
      }
    }
    rm(s)
  }
  input_sequences_oriented <- input_sequences_copy
  graph_info <- NULL

  


  kmer_counts <- list()
  kmer_selected <- c(3,4,5,6,7,11,13,17,19,23,27,31,35, 39)

  kmer_lengths <- seq(min_kmer_length, max_kmer_length, by = 4)
  kmer_lengths <- sort(intersect(kmer_selected, kmer_lengths))
  for (i in kmer_lengths) {
    ## use pythonis.null(l[[i]])) function - faster
    cmd <- paste(script.dir, "/kmer_counting.py ", input_sequences_oriented, " ",
                 i, sep = "")
    f <- system(cmd, intern = TRUE)
    x <- read.table(f, as.is = TRUE, sep = "\t")
    kmer_counts[[as.character(i)]] <- data.frame(x, freq = x$V2/sum(x$V2))
  }
  
  ## number of kmers:
  N50 <- sapply(kmer_counts, function(x) {
    sum(cumsum(x$freq) < 0.5)
  })
  
  N70 <- sapply(kmer_counts, function(x) {
    sum(cumsum(x$freq) < 0.7)
  })
  
  
  time1 <- Sys.time()
  ggmin <- mclapply(kmer_counts, get_mimimal_cc, start = 0.5, mc.cores = CPU)
  time2 <- Sys.time()
  names(ggmin) <- names(kmer_counts)
  
  
  ## estimate monomer
  monomers <- mclapply(ggmin, estimate_monomer, mc.cores = CPU)
  
  names(monomers) <- names(kmer_counts)
  monomers <- monomers[!sapply(monomers, is.null)]
  ## error handling:
  error <- sapply(monomers, class) == "try-error"
  if (any(error)) {
    cat("\nError in monomers estimation: ")
    cat("calculation failed for monomers length ", names(monomers)[error], "\n")
    print(monomers[error])
    if (any(!error)) {
      monomers <- monomers[!error]
    } else {
      stop("monomer estimation failed")
    }
  }
  
  ## summary - make function!!
  total_score <- list()
  k <- 0
  length_best <- numeric()
  score_bn <- numeric()
  consensus <- character()
  
  for (i in names(monomers)) {
    for (v in seq_along(monomers[[i]]$estimates)) {
      k <- k + 1
      total_score[[k]] <- c(kmer = as.numeric(i), variant = v, total_score = monomers[[i]]$estimates[[v]]$total_score)
      score_bn[k] <- min(rowSums(monomers[[i]]$estimates[[v]]$CM))
      length_best[k] <- monomers[[i]]$estimates[[v]]$length_variant_score[1,
                                                                          1]
      consensus[[k]] <- monomers[[i]]$estimates[[v]]$consensus2
    }
  }

  summary_table <- as.data.frame(do.call(rbind, total_score))
  summary_table$monomer_length <- length_best
  summary_table$consensus_length <- nchar(consensus)
  summary_table$score_bn <- score_bn
  summary_table$consensus <- paste("<pre>", wrap(consensus, width = 80), "<pre>",
                                   sep = "")
  consensus_DS <- DNAStringSet(consensus)
  names(consensus_DS) <- with(summary_table, paste0(kmer, "_", variant, "_sc_",
                                                    signif(total_score), "_l_", monomer_length))
  
  ## filter by score - keep
  
  ## reorder by score
  consensus_DS <- consensus_DS[order(summary_table$total_score, decreasing = TRUE)]
  summary_table <- summary_table[order(summary_table$total_score, decreasing = TRUE),
                                ]
  rownames(summary_table) <- NULL
  N <- nrow(summary_table)
  ## concatenate concensus(ie. make dimer head 2 tail) for pbs detection and other
  ## make something like 'pseudo contig' multimer for mapping - min length 200 bp

  ## searches
  consensus_DS_dimer <- DNAStringSet(paste0(consensus_DS, consensus_DS))
  tarean_contigs <- DNAStringSet(sapply(consensus_DS, function(x)
    ifelse(nchar(x)<200,
           paste(rep(as.character(x),round(300/nchar(as.character(x))+1)),collapse=''),
           as.character(x)))
    )

  names(consensus_DS_dimer) <- names(consensus_DS)
                                        # save sequences:
  consensus_DS_dimer_file <- paste0(output_dir, "/consensus_dimer.fasta")
  consensus_DS_file <- paste0(output_dir, "/consensus.fasta")
  tarean_contig_file <- paste0(output_dir, "/tarean_contigs.fasta")
  writeXStringSet(consensus_DS, consensus_DS_file)
  writeXStringSet(tarean_contigs, tarean_contig_file)
  writeXStringSet(consensus_DS_dimer, consensus_DS_dimer_file)

  dir.create(paste0(output_dir, "/img"), recursive = TRUE, showWarnings = FALSE)
  summary_table$monomer_length_graph <- numeric(N)
  summary_table$graph_image <- character(N)
  summary_table$monomer_length_logo <- numeric(N)
  summary_table$logo_image <- character(N)
  summary_table$graph_link <- character(N)
  summary_table$logo_link <- character(N)
  summary_table$n_gap50 <- numeric(N)
  summary_table$n_gap80 <- numeric(N)
  summary_table$n_gap90 <- numeric(N)

  ## export graph nd consensus estimate to cluster directory type of output may
  ## change in future
  save(ggmin, file = paste0(output_dir, "/ggmin.RData"))
  save(monomers, file = paste0(output_dir, "/monomers.RData"))
  

  for (i in 1:N) {
    kmer <- as.character(summary_table$kmer[i])
    variant <- summary_table$variant[i]
    ## export graph
    fout_link <- paste0("img/graph_", kmer, "mer_", variant, ".png")
    fout <- paste0(output_dir, "/", fout_link)
    ## summary_table$monomer_length_graph[i] = summary_table$monomer_length[i]
    ## summary_table$monomer_length_logo[[i]] = nrow(monomers[[kmer]]$estimates[[variant]]$CM)
    summary_table$monomer_length[[i]] <- length(monomers[[kmer]]$estimates[[variant]]$consensus)

    # get number of gaps in consensus fomr CMnorm
    summary_table$n_gap50[i] <-  sum(monomers[[kmer]]$estimates[[variant]]$CMnorm[,5] > 0.5)
    summary_table$n_gap80[i] <-  sum(monomers[[kmer]]$estimates[[variant]]$CMnorm[,5] > 0.8)
    summary_table$n_gap90[i] <-  sum(monomers[[kmer]]$estimates[[variant]]$CMnorm[,5] > 0.9)

    if (i <= 10) {
      png(fout, width = 800, height = 800)
      plot_kmer_graph(ggmin[[kmer]], highlight = unlist(monomers[[kmer]]$paths[[variant]]$tr_paths))
      dev.off()
      summary_table$graph_image[i] <- hwriteImage(fout_link, link = fout_link,
                                                  table = FALSE, width = 100, height = 100)
      ## export logo
      png_link <- paste0("img/logo_", kmer, "mer_", variant, ".png")
      fout <- paste0(output_dir, "/", png_link)
      png(fout, width = 1200, height = round(summary_table$monomer_length[i] * 
                                             1) + 550)
      try(plot_multiline_logo(monomers[[kmer]]$estimates[[variant]]$CM, W = 100))
      dev.off()
      ## export corresponding position probability matrices
      ppm_file <- paste0(output_dir, '/ppm_', kmer, "mer_", variant, ".csv")
      ppm_link <- paste0('ppm_', kmer, "mer_", variant, ".csv")
      write.table(monomers[[kmer]]$estimates[[variant]]$CM,
                  file = ppm_file,
                  col.names = TRUE, quote = FALSE,
                  row.names = FALSE, sep="\t")
      summary_table$logo_image[i] <- hwriteImage(png_link, link = ppm_link,
                                                 table = FALSE, width = 200, height = 100)
      summary_table$graph_link[i] <- fout_link
      summary_table$logo_link[i] <-  png_link
    }
    
  }
  
  ## html_report = HTMLReport()
  
  htmlfile <- paste0(output_dir, "/report.html")
  cat(htmlheader, file = htmlfile)
  included_columns <- c('kmer', 'variant', 'total_score', 'consensus_length', 'consensus', 'graph_image', 'logo_image')
  summary_table_clean <- summary_table[, included_columns]
  colnames(summary_table_clean) <- c('k-mer length',
                                     'Variant index',
                                     'k-mer coverage score',
                                     'Consensus length',
                                     'Consensus sequence',
                                     'k-mer based graph',
                                     'Sequence logo')
  HTML(summary_table_clean, file = htmlfile, sortableDF = TRUE)
  HTMLEndFile(file = htmlfile)
  # export summary table as csv
  summary_table_file <- paste0(output_dir, "/summary_table.csv")
  write.csv(summary_table, summary_table_file, row.names = FALSE, quote = TRUE)
  time4 <- Sys.time()

  if (!is.null(lock_file)){
    print("------removing-lock--------")
    removelock(lock_file)
  }
  return(list2dictionary(list(htmlfile = htmlfile, TR_score = summary_table$total_score[1],
              TR_monomer_length = as.numeric(summary_table$consensus_length[1]),
              TR_consensus = summary_table$consensus[1], graph_info = graph_info,
              tarean_contig_file = tarean_contig_file)))
}


## graph loop index stability
loop_index_instability <- function(G) {
  N <- 50
  s <- seq(vcount(G), vcount(G)/10, length.out = N)
  p <- seq(1, 0.1, length.out = N)
  li <- numeric()
  for (i in seq_along(s)) {

    gs <- induced_subgraph(G, sample(1:vcount(G), s[i]))
    li[i] <- max(clusters(gs, "strong")$csize)/vcount(gs)
  }
  instability <- lm(li ~ p)$coefficient[2]
  return(instability)
}

isSatellite <- function(x, y, model) {
  p <- get_prob(x, y, model)
  if (p > model$cutoff) {
    return("Putative Satellite")
  } else {
    return("")
  }
}

get_prob <- function(x, y, model) {
  pm <- model$prob_matrix
  N <- ncol(pm)
  i <- round(x * (N - 1)) + 1
  j <- round(y * (N - 1)) + 1
  p <- pm[i, j]
  return(p)
}


detectMemUsage <- function() {
  con <- textConnection(gsub(" +", " ", readLines("/proc/meminfo")))
  memInfo <- read.table(con, fill = TRUE, row.names = 1)
  close(con)
  memUsage <- 1 - (memInfo["MemFree", 1] + memInfo["Cached", 1])/memInfo["MemTotal",
                                                                         1]
  return(memUsage)
}


makelock<-function(lockfile,lockmsg,CreateDirectories=TRUE){
    lockdir <- dirname(lockfile)
    if(!file.exists(lockdir)){
        if(CreateDirectories) dir.create(lockdir,recursive=TRUE, showWarnings=FALSE)
        else stop("Lock Directory for lockfile ",lockfile," does not exist")
    } 
    if(missing(lockmsg)) lockmsg <- paste(system('hostname', intern=TRUE), Sys.getenv("R_SESSION_TMPDIR"))
    if (file.exists(lockfile)) return (FALSE)
                                        # note the use of paste makes the message writing atomic
    cat(paste(lockmsg,"\n",sep=""),file=lockfile,append=TRUE,sep="")
    firstline <- readLines(lockfile, n=1)
    if(firstline!=lockmsg){
                                        # somebody else got there first
        return(FALSE)
    } else return(TRUE)
}


removelock<-function(lockfile){
  if(unlink(lockfile)!=0) {
    warning("Unable to remove ",lockfile)
    return (FALSE)
  }
  return (TRUE)
}


waitForRAM <- function(p = 0.5, lock_file=NULL) {
  if (detectMemUsage() < p) {
    return(NULL)
    ## check lock file:
  } else {
    cat("waiting for RAM \n")
    free_count <- 0
    while (TRUE) {
        if (makelock(lock_file)){
            print("---------locking--------")
            return(lock_file)
        }
      if (detectMemUsage() < p) {
        cat("RAM freed \n")
        return(NULL)
      }
      Sys.sleep(5)
      if (evaluate_user_cpu_usage() == 'free'){
        free_count <- free_count + 1
      }else{
        free_count <- 0
      }
      if (detectMemUsage() < 0.8 & free_count > 100){
        cat("RAM not free but nothing else is running \n")
        return(NULL)
      }
    }
  }
}

lsmem <- function() {
  g <- globalenv()
  out_all <- envs <- list()
  envs <- append(envs, g)
  total_size <- numeric()
  while (environmentName(g) != "R_EmptyEnv") {
    g <- parent.env(g)
    envs <- append(envs, g)
  }
  for (e in envs) {
    
    obj <- ls(envir = e)
    if (length(obj) == 0) {
      break
    }
    obj.size <- list()
    for (i in obj) {
      obj.size[[i]] <- object.size(get(i, envir = e))
    }
    out <- data.frame(object = obj, size = unlist(obj.size), stringsAsFactors = FALSE)
    out <- out[order(out$size, decreasing = TRUE), ]
    out_all <- append(out_all, out)
    total_size <- append(total_size, sum(out$size))
  }
  return(list(objects = out_all, total_size = total_size))
} 

evaluate_user_cpu_usage <- function(){
  user <- Sys.info()["user"]
  a <- sum(as.numeric (system(paste ("ps -e -o %cpu -u", user), intern = TRUE)[-1]))
  s <- substring (system(paste ("ps -e -o stat -u", user), intern = TRUE)[-1], 1, 1)
  if (a<5 & sum(s %in% 'D')==0 & sum(s%in% 'R')<2){
    status <- 'free'
  }else{
    status <- 'full'
  }
  return(status)
}
