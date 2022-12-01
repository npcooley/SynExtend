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
  
  allsp <- lapply(labvecs, \(x) gsub('(.*)_.*$', '\\1', x))
  n <- names(pw)[uvals]
  
  ARGS <- list(allsp=allsp)
  FXN <- function(lab1, lab2, ARGS, ii, jj){
    l1sp <- gsub('(.*)_.*$', '\\1', lab1)
    l2sp <- gsub('(.*)_.*$', '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    score <- 0
    mult <- ifelse(length(shared)==0, NA, 1/length(shared))
    for ( k in seq_along(shared) ){
      spk <- shared[k]
      vec1 <- lab1[match(spk, l1sp)]
      vec2 <- lab2[match(spk, l2sp)]
      idx1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec1, useBytes=TRUE))
      idx2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec2, useBytes=TRUE))
      # Double loop is SIGNIFICANTLY faster than expand.grid
      # ...even though it looks so gross :(
      diffs <- 0
      for (v1i in idx1)
        for (v2i in idx2)
          diffs <- diffs + exp(1 - abs(v1i - v2i))
      score <- score + sum(diffs) / (length(vec1) * length(vec2))
    }
    return(score * mult)
  }
  pairscores <- BuildSimMatInternal(labvecs, uvals, evalmap, l, n, 
                                    FXN, ARGS, Verbose,
                                    InputIsList=TRUE)
  
  m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  pairscores <- pairscores / m
  Diag(pairscores) <- 1
  return(pairscores)
}
