###### Residue Methods for ProtWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Residue Mutual Information
#  - Natural Vector + Di/Trinucleotide Freq.
##########################

#### S3 Generic Definitions ####
ResidueMI <- function(pw, ...) UseMethod('ResidueMI')
NVDT <- function(pw, ...) UseMethod('NVDT')
################################

ResidueMI.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                 precalcSubset=NULL, ...){
  useResidue <- attr(pw, 'useResidue')
  useMT <- attr(pw, 'useMT')
  useColoc <- attr(pw, 'useColoc')
  
  stopifnot('ProtWeaver object must be initialized with dendrograms to run Residue methods'=
              useMT)
  stopifnot('ProtWeaver dendrograms must have ancestral states to run Residue methods'=
              useResidue)
  
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  
  pairscores <- rep(NA_real_, l*(l-1)/2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    tree1 <- pw[[uval1]]
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        tree2 <- pw[[uval2]]
        pairscores[ctr+1] <- ResidueMIDend(tree1, tree2, useColoc=useColoc, ...)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  n <- n[uvals]
  pairscores <- as.simMat(pairscores, NAMES=n, DIAG=FALSE)
  Diag(pairscores) <- 1
  if (Verbose) cat('\n')
  
  return(pairscores)
}

NVDT.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                 precalcSubset=NULL, extended=TRUE, ...){
  #source('/Users/aidan/Nextcloud/RStudioSync/comps/NVDT/calcNVDT.R')
  useResidue <- attr(pw, 'useResidue')
  useMT <- attr(pw, 'useMT')
  
  stopifnot('ProtWeaver object must be initialized with dendrograms to run Residue methods'=
              useMT)
  stopifnot('ProtWeaver dendrograms must have ancestral states to run Residue methods'=
              useResidue)
  
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  outl <- ifelse(extended, 92L, 12L)
  vecs <- matrix(NA_real_, nrow=l, ncol=outl)
  if(Verbose){
    cat('Calculating sequence vectors...\n')
    pb <- txtProgressBar(max=l, style=3)
  }
  for (i in seq_len(l)){
    uv <- uvals[i]
    tree <- pw[[uv]]
    #seqs <- DNAStringSet(flatdendrapply(tree, NODEFUN=\(x) attr(x, 'state')))
    seqs <- flatdendrapply(tree, NODEFUN=\(x) attr(x, 'state'))
    # Remove Gaps
    #seqs <- gsub('[-.+]', '', seqs)
    #s1 <- gsub('[-.+]', '', seqs)
    #seqs <- gsub('.', '', seqs, fixed=TRUE, useBytes=TRUE)
    #if(!all(s1==s2)){
    #  which(s1!=s2)
    #  print("here")
    #}
    #seqs <- gsub('.', '', seqs, fixed=TRUE, useBytes=TRUE)
    outv <- numeric(outl)
    for (j in seq_along(seqs)){
      outv <- outv + .Call("StringToNVDC", 
                           charToRaw(seqs[j]),
                           FALSE,
                           extended,
                           PACKAGE="SynExtend")
    }
    outv <- outv / sqrt(sum(outv**2))
    vecs[i,] <- outv
    if(Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\nDone.\n')
  pairscores <- rep(NA_real_, l*(l-1)/2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    v1 <- vecs[i,]
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        v2 <- vecs[j,]
        pairscores[ctr+1] <- cor(v1,v2)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  n <- n[uvals]
  pairscores <- as.simMat(pairscores, NAMES=n, DIAG=FALSE)
  Diag(pairscores) <- 1
  if (Verbose) cat('\n')
  
  return(pairscores)
  
}