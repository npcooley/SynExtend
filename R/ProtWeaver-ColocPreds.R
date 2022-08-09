###### Co-localization Methods for ProtWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Naive Coloc
##########################

#### S3 Generic Definitions ####
Coloc <- function(pw, ...) UseMethod('Coloc')
################################

Coloc.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                             precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  stopifnot('Colocalization is disabled.'=attr(pw,'useColoc'))
  if (attr(pw, 'useMT')){
    labvecs <- lapply(pw[uvals], labels)
  } else {
    labvecs <- pw[uvals]
  }
  l <- length(labvecs)
  n <- names(pw)[uvals]
  pairscores <- matrix(NA, nrow=l, ncol=l)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l)){
    uval1 <- uvals[i]
    lab1 <- labvecs[[i]]
    l1sp <- gsub('(.*)_.*$', '\\1', lab1) # species + chromosome number
    for ( j in i:l){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        lab2 <- labvecs[[j]]
        l2sp <- gsub('(.*)_.*$', '\\1', lab2)
        
        shared <- intersect(l1sp, l2sp)
        score <- 0
        mult <- ifelse(length(shared)==0, NA, 1/length(shared))
        for ( k in seq_along(shared) ){
          spk <- shared[k]
          v1 <- lab1[match(spk, l1sp)]
          v2 <- lab2[match(spk, l2sp)]
          idx1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', v1, useBytes=TRUE))
          idx2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', v2, useBytes=TRUE))
          # Double loop is SIGNIFICANTLY faster than expand.grid
          # ...even though it looks so gross :(
          diffs <- 0
          for (v1i in idx1)
            for (v2i in idx2)
              diffs <- diffs + exp(1 - abs(v1i - v2i))
          score <- score + sum(diffs) / (length(v1) * length(v2))
        }
        score <- score * mult
        pairscores[i,j] <- pairscores[j,i] <- score
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  pairscores <- pairscores / m
  diag(pairscores) <- 1
  rownames(pairscores) <- colnames(pairscores) <- n
  return(pairscores)
}
