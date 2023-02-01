SuperTree <- function(myDendList, NAMEFUN=NULL, Verbose=TRUE, Processors=1){
  # error checking
  if (!is(myDendList, 'list') || 
      !all(vapply(myDendList, \(x) is(x, 'dendrogram') || is.null(x), FUN.VALUE=TRUE))){
    stop("SuperTree requires input to be a list of dendrograms")
  }
  
  stopifnot("'Processors' must be a positive integer"=
              (Processors >=1 && Processors == as.integer(Processors)))
  if (!is(Processors, 'integer')) Processors <- as.integer(Processors)
  
  useNFUN <- is(NAMEFUN, 'function')
  
  if (useNFUN){
    l1 <- flatdendrapply(myDendList[[1]], LEAFFUN=attr, which='label')
    l2 <- NAMEFUN(l1)
    if(length(l1)!=length(l2) || !all(is.atomic(l2))){
      stop("NAMEFUN should operate on a character vector and return a character
           vector of the same size as input.")
    }
  }
  
  # Get list of species
  allspecies <- character(0)
  if (Verbose){ 
    cat("  Compiling species...\n")
    pb = txtProgressBar(max=length(myDendList), style=3)
    start <- Sys.time()
  }
  ctr <- 0
  for ( i in seq_along(myDendList) ){
    if (Verbose) setTxtProgressBar(pb, i)
    lst <- myDendList[[i]]
    if (is.null(lst)) next
    specs <- flatdendrapply(lst, LEAFFUN=\(x) attr(x, 'label'))
    
    if (useNFUN)
      specs <- NAMEFUN(specs)
    
    allspecies <- unique(c(allspecies, specs))
    ctr <- ctr + 1
  }
  
  # Initialize distance matrix and counts matrix
  dmat <- countmat <- matrix(0, nrow=length(allspecies), ncol=length(allspecies))
  rownames(dmat) <- colnames(dmat) <- allspecies
  rownames(countmat) <- colnames(countmat) <- allspecies
  if (Verbose){ 
    cat("\n  Done.\n\n  Calculating distance matrices...\n")
    pb <- txtProgressBar(max = ctr, style=3)
  }
  ctr <- 0
  for ( dend in myDendList ){
    if (is.null(dend)) next
    # Calculate Cophenetic
    d <- as.matrix(Cophenetic(dend))
    
    # Get rownames
    if (useNFUN)
      rownames(d) <- colnames(d) <- NAMEFUN(rownames(d))
    rn <- rownames(d)
    cn <- colnames(d)
    lu <- length(unique(rn))
    lr <- length(rn)
    # Paralogs (currently averaging all rows corresponding to same genome)
    while (lu != lr){
      firsttocom <- names(which.max(table(rownames(d)) > 1))
      pos <- which(rownames(d) == firsttocom)
      newrow <- colSums(d[pos,]) / length(pos)
      d[pos[1], ] <- d[,pos[1]] <- newrow
      d <- d[-pos[2], -pos[2]]
      lr <- lr - 1
      rn <- rn[-pos[2]]
    }
    cn <- rn
    dmat[rn, cn] <- dmat[rn, cn] + d
    countmat[rn, cn] <- countmat[rn, cn] + 1
    ctr <- ctr + 1
    if (Verbose) setTxtProgressBar(pb, ctr)
  }
  
  # Average result
  dmat <- as.dist(dmat / countmat)
  
  # Build species tree with NJ
  if (Verbose) cat("\n  Building species tree...\n")
  newTree <- TreeLine(myDistMatrix=dmat, method="NJ", 
                      verbose=Verbose, processors = Processors)
  
  if (Verbose){
    dt <- difftime(start, Sys.time())
    cat("Done.\n  Time difference of", round(abs(dt), 2), attr(dt, "units"), '\n', sep=' ')
  }
  return(newTree)
}
