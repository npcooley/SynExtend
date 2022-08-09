###### Residue Methods for ProtWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Naive Coloc
##########################

#### S3 Generic Definitions ####
ResidueMI <- function(pw, ...) UseMethod('ResidueMI')
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
  
  pairscores <- matrix(NA_real_, nrow=l, ncol=l)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l) ){
    uval1 <- uvals[i]
    tree1 <- pw[[uval1]]
    for ( j in i:l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (i!=j && (is.null(evalmap) || entry %in% evalmap[[accessor]])){
        tree2 <- pw[[uval2]]
        pairscores[i,j] <- pairscores[j,i] <- ResidueMIDend(tree1, tree2, 
                                                            useColoc=useColoc, ...)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  diag(pairscores) <- 1
  if (Verbose) cat('\n')
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  
  return(pairscores)
}