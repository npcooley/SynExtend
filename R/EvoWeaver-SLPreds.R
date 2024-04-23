###### Residue Methods for EvoWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Residue Mutual Information
#  - Natural Vector + Di/Trinucleotide Freq.
##########################

#### S3 Generic Definitions ####
ResidueMI <- function(ew, ...) UseMethod('ResidueMI')
NVDT <- function(ew, ...) UseMethod('NVDT')
Ancestral <- function(ew, ...) UseMethod('Ancestral')
################################

ResidueMI.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                 precalcSubset=NULL, gapCutoff=0.5,
                                 useDNA=FALSE, Processors=1L, useWeights=TRUE, ...){
  useResidue <- attr(ew, 'useResidue')
  useMT <- attr(ew, 'useMT')
  useColoc <- attr(ew, 'useColoc')
  Processors <- NormArgProcessors(Processors)

  stopifnot('EvoWeaver object must be initialized with dendrograms to run Residue methods'=
              useMT)
  stopifnot('EvoWeaver dendrograms must have ancestral states to run Residue methods'=
              useResidue)

  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)

  uvals <- subs$uvals
  evalmap <- subs$evalmap
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  if(useDNA)
    lookup <- 'ATGC'
  else
    lookup <- "ARNDCQEGHILKMFPSTWYV"
  lookup <- strsplit(lookup, '')[[1]]
  lookupMap <- seq_along(lookup)
  names(lookupMap) <- lookup

  ResidueSets <- vector('list', length(uvals))
  if(Verbose){
    cat('  Preproccessing sequence sets\n')
    pb <- txtProgressBar(max=length(uvals), style=3)
  }
  for(i in seq_along(uvals)){
    tree <- ew[[uvals[i]]]
    seqs <- rapply(tree, attr, which='state')
    ncharSeq <- nchar(seqs)
    stopifnot('Sequences are not aligned!'=all(ncharSeq==ncharSeq[1]))
    seqs <- toupper(seqs)
    if(useDNA){
      maskPos <- which(!MaskAlignment(DNAStringSet(seqs),
                                      type='values', includeTerminalGaps = TRUE)[,3])
    } else {
      maskPos <- which(!MaskAlignment(AAStringSet(seqs), threshold=0.3,
                                      type='values', includeTerminalGaps = TRUE)[,3])
    }

    if(length(maskPos)==0){
      seqs <- matrix(NA_integer_, nrow=length(seqs), ncol=0L)
    } else {
      seqs <- vapply(seqs,
                     \(s) lookupMap[(strsplit(s, '')[[1]])[maskPos]],
                                integer(length(maskPos)), USE.NAMES = FALSE)
      seqs[is.na(seqs)] <- 0L
      seqs <- t(seqs)
    }
    labs <- labels(tree)
    if(useColoc){
      labs <- gsub('([^_]*)_.*', '\\1', labs)
    }
    rownames(seqs) <- labs
    ResidueSets[[i]] <- seqs
    if(Verbose) setTxtProgressBar(pb, i)
  }

  pairscores <- rep(NA_real_, l*(l-1)/2)
  ctr <- 0
  if (Verbose){
    cat('\nDone.\n')
    pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  }
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    #tree1 <- ew[[uval1]]
    seqs1 <- ResidueSets[[i]]
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        #tree2 <- ew[[uval2]]
        seqs2 <- ResidueSets[[j]]
        #pairscores[ctr+1] <- ResidueMIDend(tree1, tree2, useColoc=useColoc, ...)
        score <- ResidueMISeqs(seqs1, seqs2, lookup=lookup,
                                useColoc=useColoc, Processors=Processors, ...)
        pairscores[ctr+1] <- ifelse(is.na(score), 0L, score)
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

NVDT.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                            precalcSubset=NULL, extended=TRUE,
                            DNAseqs=TRUE, centerObservations=FALSE,
                            sqrtCorrelation=TRUE, ...){
  #source('/Users/aidan/Nextcloud/RStudioSync/comps/NVDT/calcNVDT.R')
  useResidue <- attr(ew, 'useResidue')
  useMT <- attr(ew, 'useMT')

  stopifnot('EvoWeaver object must be initialized with dendrograms to run Residue methods'=
              useMT)
  stopifnot('EvoWeaver dendrograms must have ancestral states to run Residue methods'=
              useResidue)

  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)

  uvals <- subs$uvals
  evalmap <- subs$evalmap
  l <- length(uvals)
  alllabs <- lapply(uvals, \(x) labels(ew[[x]]))
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  outl <- ifelse(DNAseqs, 12L, 60L)
  if(extended && DNAseqs)
    outl <- 92L
  vecs <- matrix(NA_real_, nrow=l, ncol=outl)
  if(Verbose){
    cat('Calculating sequence vectors...\n')
    pb <- txtProgressBar(max=l, style=3)
  }
  for (i in seq_len(l)){
    uv <- uvals[i]
    tree <- ew[[uv]]
    #seqs <- DNAStringSet(flatdendrapply(tree, NODEFUN=\(x) attr(x, 'state')))
    seqs <- rapply(tree, attr, which='state')
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
      nvdtvec <- .Call("StringToNVDT",
                       charToRaw(seqs[j]),
                       FALSE,    # RemoveGaps
                       extended, # Use di/tri freqs
                       DNAseqs,     # Use DNA base pairs
                       PACKAGE="SynExtend")
      seqlen <- nchar(seqs[j])
      if(DNAseqs){
        divval <- c(rep(seqlen, 8L), rep(1L, 4L))
        if(extended){
          #correcting 2-mers and 3-mers
          divval <- c(divval, rep(seqlen-1L, 16L), rep(seqlen-2L, 64L))
        }
      } else {
        divval <- rep(seqlen, 60L)
      }
      nvdtvec <- nvdtvec / divval
      outv <- outv + nvdtvec
    }
    if(centerObservations && !DNAseqs){
      s <- seq_len(20L)
      outv[s] <- (outv[s] - mean(outv[s])) / sd(outv[s])
      # Center should always be 0.5 for mean
      s <- s + 20L
      outv[s] <- (outv[s] - 0.5) / sd(outv[s])
      s <- s + 20L
      outv[s] <- (outv[s] - mean(outv[s])) / sd(outv[s])
    } else {
      outv <- outv / sqrt(sum(outv**2))
    }
    vecs[i,] <- outv
    if(Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\nDone.\n')
  pairscores <- rep(NA_real_, l*(l-1)/2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    l1 <- alllabs[[i]]
    if(sqrtCorrelation)
      v1 <- sqrt(abs(vecs[i,])) * sign(vecs[i,])
    else
      v1 <- vecs[i,]
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      l2 <- alllabs[[j]]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        if(sqrtCorrelation)
          v2 <- sqrt(abs(vecs[j,])) * sign(vecs[j,])
        else
          v2 <- vecs[j,]
        if(length(intersect(l1,l2))==0){
          pairscores[ctr+1] <- 0
        } else {
          # first entry is score, second is t-value
          # cval <- .Call('fastPearsonC', v1, v2)
          # pval <- pt(cval[2], cval[3]-2)
          # pval <- ifelse(pval < 0.5, 2*pval, 2*(pval-0.5))

          cval <- suppressWarnings(cor(v1, v2,
                                      use='pairwise.complete.obs',
                                      method='spearman'))
          # should be absolute value, negative correlation is same effect
          cval <- abs(cval[1])

          num_obs <- length(v1)
          pval <- 1 - exp(pt(cval, num_obs-2, lower.tail=FALSE, log.p=TRUE))
          #pairscores[ctr+1] <- cor(v1,v2)
          pairscores[ctr+1] <- pval*cval
        }
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

Ancestral.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                 precalcSubset=NULL, ...){
  useResidue <- attr(ew, 'useResidue')
  useMT <- attr(ew, 'useMT')
  useColoc <- attr(ew, 'useColoc')

  stopifnot('EvoWeaver object must be initialized with dendrograms to run Residue methods'=
              useMT)
  stopifnot('EvoWeaver dendrograms must have ancestral states to run Residue methods'=
              useResidue)

  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)

  uvals <- subs$uvals
  evalmap <- subs$evalmap
  l <- length(uvals)
  n <- names(ew)

  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  pmlst <- vector('list', length=l)
  if (Verbose){
    cat("Pre-processing dendrograms...\n")
    pb <- txtProgressBar(max=l, style=3)
  }
  for (i in seq_len(l)){
    uv <- uvals[i]
    pmlst[[i]] <- find_dists_pos(ew[[uv]], useColoc)
    if(Verbose) setTxtProgressBar(pb, i)
  }

  if(Verbose) cat('\nDone.\n')
  pairscores <- rep(NA_real_, l*(l-1)/2)
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  ctr <- 0
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    v1 <- pmlst[[i]]
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        v2 <- pmlst[[j]]
        res <- pair_residues(v1,v2)
        if(is.na(res[1]) || is.nan(res$R))
          pairscores[ctr+1] <- NA
        else
          pairscores[ctr+1] <- abs(res$R) * (1-res$P)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
      #if (Verbose) cat("\r", ctr, " / ", 2963, sep='')
    }
  }
  n <- n[uvals]
  pairscores <- as.simMat(pairscores, NAMES=n, DIAG=FALSE)
  Diag(pairscores) <- 1
  if (Verbose) cat('\n')

  return(pairscores)
}
