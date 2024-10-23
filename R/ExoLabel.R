ExoLabel <- function(edgelistfiles, outfile=tempfile(),
                          mode=c("undirected", "directed"),
                          add_self_loops=FALSE,
                          ignore_weights=FALSE,
                          normalize_weights=FALSE,
                          iterations=0L, inflation=1.05,
                          return_table=FALSE,
                          consensus_cluster=FALSE,
                          verbose=interactive(),
                          sep='\t',
                          tempfiledir=tempdir(),
                          cleanup_files=TRUE){
  on.exit(.C("cleanup_ondisklp_global_values"))
  if(!is.numeric(iterations)){
    stop("'iterations' must be an integer or numeric.")
  } else {
    iterations <- as.integer(iterations)
  }
  if(!is.numeric(inflation)){
    stop("'inflation' must be a numeric value")
  } else {
    inflation <- as.numeric(inflation)
  }
  if(is.na(iterations) || is.null(iterations) || is.infinite(iterations) || iterations < 0){
    warning("Invalid value of 'iterations', will determine automatically from node degree.")
    iterations <- 0L
  }
  if(is.infinite(inflation) || is.na(inflation) || inflation < 0){
    warning("Invalid value of 'inflation', defaulting to 1.0")
    inflation <- 1.0
  }
  if(!is.numeric(add_self_loops) && !is.logical(add_self_loops)){
    stop("value of 'add_self_loops' should be numeric or logical")
  }
  if(add_self_loops < 0){
    warning("self loops weight supplied is negative, setting to zero.")
    add_self_loops <- 0
  } else if(is.logical(add_self_loops)){
    add_self_loops <- ifelse(add_self_loops, 1, 0)
  }
  if(!is.logical(ignore_weights)){
    stop("'ignore_weights' must be logical")
  } else if(is.na(ignore_weights) || is.null(ignore_weights)){
    stop("invalid value for 'ignore_weights' (should be TRUE or FALSE)")
  }
  if(!is.logical(normalize_weights)){
    stop("'normalize_weights' must be logical")
  } else if(is.na(normalize_weights) || is.null(normalize_weights)){
    stop("invalid value for 'normalize_weights' (should be TRUE or FALSE)")
  }
  if(ignore_weights && normalize_weights){
    warning("Cannot both ignore weights and normalize them")
  }
  # verify that the first few lines of each file are correct
  if(!all(file.exists(edgelistfiles))) stop("edgelist file does not exist")
  edgelistfiles <- normalizePath(edgelistfiles, mustWork=TRUE)
  for(f in edgelistfiles){
    v <- readLines(f, n=10L)
    v <- strsplit(v, sep)
    lv <- lengths(v)
    if(any(lv != lv[1]) || lv[1] < 2) stop("file ", f, " is misformatted")
    lv <- lv[1] # now we know they're all the same
    if(!ignore_weights && lv == 2) stop("file ", f, " is missing weights!")
    if(!ignore_weights && any(vapply(v, \(x) is.na(as.numeric(x[3])), logical(1L))))
      stop("file ", f, " has malformed weights")
  }

  if(is.logical(consensus_cluster)){
    if(consensus_cluster){
      consensus_cluster <- c(0,0.2,0.4,0.6,0.8,1,1.33,1.67,2)
    } else {
      consensus_cluster <- numeric(0L)
    }
  } else {
    if(!is.numeric(consensus_cluster) || any(is.na(consensus_cluster) | is.null(consensus_cluster)))
      stop("'consensus_cluster' must be a logical or numeric vector")
    if(any(consensus_cluster < 0))
      stop("'consensus_cluster' cannot contain negative values")
  }
  tempfiledir <- normalizePath(tempfiledir, mustWork=TRUE)
  tempfiledir <- file.path(tempfiledir, "ExoLabelTemp")
  if(dir.exists(tempfiledir)){
    for(f in list.files(tempfiledir, full.names=TRUE))
      file.remove(f)
  } else {
    dir.create(tempfiledir)
  }
  mode <- match.arg(mode)
  is_undirected <- mode == "undirected"
  outfile <- file.path(normalizePath(dirname(outfile), mustWork=TRUE), basename(outfile))

  if(verbose) cat("Temporary files stored at ", tempfiledir, "\n")

  seps <- paste(sep, "\n", sep='')
  ctr <- 1
  # R_hashedgelist(tsv, csr, clusters, queues, hashdir, seps, 1, iter, verbose)
  .Call("R_LPOOM_cluster", edgelistfiles, length(edgelistfiles),
        tempfiledir, outfile, seps, ctr, iterations,
        verbose, is_undirected, add_self_loops, ignore_weights, normalize_weights,
        consensus_cluster, inflation)

  # R_write_output_clusters(clusters, hashes, length(hashes), out_tsvpath, seps)
  #.Call("R_LP_write_output", counter_cluster_binary, hashdir,
  #      outfile, seps, verbose)
  if(cleanup_files){
    for(f in list.files(tempfiledir, full.names=TRUE))
      if(file.exists(f)) file.remove(f)
    file.remove(tempfiledir)
  }
  if(return_table){
    tab <- read.table(outfile, sep=sep)
    colnames(tab) <- c("Vertex", "Cluster")
    if(file.exists(outfile)) file.remove(outfile)
    return(tab)
  } else {
    return(outfile)
  }
}

EstimateExoLabel <- function(num_v, avg_degree=2,
                          num_edges=num_v*avg_degree, node_name_length=10L){
  if(!missing(avg_degree) && !missing(num_edges)){
    warning("Only one of 'avg_degree' and 'num_edges' are needed, ignoring num_edges")
  } else if (missing(avg_degree)){
    avg_degree = num_edges / num_v
  }
  lv <- num_v*node_name_length
  # assuming file is v1 v2 %.3f, which is 2*node_name_len + 3 + 5
  exp_size_file <- (2*node_name_length+8)*num_edges
  exp_size_internal <- 41*num_v+12*num_edges
  exp_size_final <- (2+node_name_length+log10(num_v))*num_v
  exp_size_ram_lower <- (24 + 16)*num_v + 104857600 # 1e8 is roughly the cache size
  exp_size_ram_upper <- (24*node_name_length + 16)*num_v + 104857600
  exp_ratio <- exp_size_internal / exp_size_file
  v <- c(exp_size_ram_lower, exp_size_ram_upper, exp_size_file, exp_size_internal, exp_size_final, exp_ratio)
  names(v) <- c("Minimum RAM Usage", "Maximum RAM Usage", "Expected Input File Size", "Expected Internal File Size",
    "Expected Final File Size", "Disk Usage Ratio")

  max_nchar <- max(nchar(names(v)[-length(v)]))
  unitsizes <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB")
  for(i in seq_along(v)){
    if(i == length(v)){
      if(exp_ratio < 0.001){
        cat("\nExoLabel total disk consumption is <0.001x that of the original files\n")
      } else {
        cat("\nExoLabel total disk consumption is about ", round(exp_ratio, 2), "x that of the initial files.\n", sep='')
      }
      next # skip other stats for disk usage ratio
    }
    unit <- ""
    for(j in seq_along(unitsizes)){
      unit <- unitsizes[j]
      p <- 1024^(j-1)
      if((v[i] / p) < 1024)
        break
    }
    n <- names(v)[i]
    padding_required <- max_nchar - nchar(n)
    pad <- paste(rep(' ', padding_required), collapse='')
    cat(pad, names(v)[i], ": ", sprintf("%5.1f ", v[i]/p), unit, '\n', sep='')
  }
  invisible(v)
}
