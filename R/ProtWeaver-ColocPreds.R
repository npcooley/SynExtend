###### Co-localization Methods for ProtWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Naive Coloc
##########################

#### S3 Generic Definitions ####
Coloc <- function(pw, ...) UseMethod('Coloc')
ColocMoran <- function(pw, ...) UseMethod('ColocMoran')
TranscripMI <- function(pw, ...) UseMethod('TranscripMI')
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

ColocMoran.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                  MySpeciesTree=NULL,
                                  precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  stopifnot('Colocalization is disabled.'=attr(pw,'useColoc'))
  if ( is.null(MySpeciesTree) || !is(MySpeciesTree, 'dendrogram')){
    MySpeciesTree <- findSpeciesTree(pw, Verbose)
  }
  stopifnot('Missing MySpeciesTree'=!is.null(MySpeciesTree))
  stopifnot('MySpeciesTree must be a dendrogram'=is(MySpeciesTree, 'dendrogram'))
  if (attr(pw, 'useMT')){
    labvecs <- lapply(pw[uvals], labels)
  } else {
    labvecs <- pw[uvals]
  }
  l <- length(labvecs)
  specCoph <- as.matrix(Cophenetic(MySpeciesTree))
  
  allsp <- lapply(labvecs, \(x) gsub('(.*)_.*$', '\\1', x))
  n <- names(pw)[uvals]
  
  ARGS <- list(allsp=allsp, cMat=specCoph)
  FXN <- function(lab1, lab2, ARGS, ii, jj){
    l1sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab1)
    l2sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    if (length(shared) <= 3){
      return(0)
    }
    score <- 0
    mult <- 1/length(shared)
    vals <- rep(NA_real_, length(shared))
    for ( k in seq_along(shared) ){
      spk <- shared[k]
      vec1 <- lab1[match(spk, l1sp)]
      vec2 <- lab2[match(spk, l2sp)]
      idx1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec1, useBytes=TRUE))
      idx2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec2, useBytes=TRUE))
      diffs <- 0
      vals[k] <- min(abs(c(outer(as.integer(idx1), as.integer(idx2), '-'))))
      vals[k] <- exp(-vals[k])
      # exp(-sqrt(vals[k])) is a little better but feels overfit-y
    }
    if(all(vals==vals[1])){
      # if there's 100% consistency it'll break MoransI,
      # but honestly at that point isn't that great
      # evidence of colocalization?
      return(1)
    }
    shared <- gsub('([^_]*).*', '\\1', shared)
    w <- ARGS$cMat
    w <- w[shared, shared]
    w <- as.dist(exp(-w))
    # two.sided or less ?
    # 1-greater is better
    res <- MoransI(vals, w, alternative = 'greater')
    score <- mean(vals)
    pval <- res$p.value
    return(score*pval)
  }
  
  pairscores <- BuildSimMatInternal(labvecs, uvals, evalmap, l, n, 
                                    FXN, ARGS, Verbose,
                                    InputIsList=TRUE)
  
  #m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  #pairscores <- pairscores / m
  Diag(pairscores) <- 1
  return(pairscores)
}

TranscripMI.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                          precalcProfs=NULL, precalcSubset=NULL, ...){
  stopifnot('Some labels are missing strand identifiers!'=attr(pw, 'useStrand'))
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
    l1sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab1)
    l2sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    score <- 0
    mult <- ifelse(length(shared)==0, NA, 1/length(shared))
    vals <- rep(NA_real_, length(shared))
    conttable <- matrix(rep(1,4), nrow=2)
    for ( k in seq_along(shared) ){
      spk <- shared[k]
      vec1 <- lab1[match(spk, l1sp)]
      vec2 <- lab2[match(spk, l2sp)]
      idx1 <- as.integer(gsub('.*_([01])_[0-9]*$', '\\1', vec1, useBytes=TRUE))
      idx2 <- as.integer(gsub('.*_([01])_[0-9]*$', '\\1', vec2, useBytes=TRUE))
      gene1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec1, useBytes=TRUE))
      gene2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec2, useBytes=TRUE))
      for (x in seq_along(idx1)){
        if(gene1[x] > gene2[x]){
          v1 <- idx2[x]+1
          v2 <- idx1[x]+1
        } else {
          v1 <- idx1[x]+1
          v2 <- idx2[x]+1
        }
        conttable[v1,v2] <- conttable[v1,v2]+1
      }
    }
    totalc <- sum(conttable)
    px1 <- sum(conttable[1,]) / totalc
    py1 <- sum(conttable[1,]) / totalc
    pmf <- conttable / totalc
    mutinf <- pmf[1,1]*log2(pmf[1,1] / (px1*py1)) +
      -1*pmf[1,2]*log2(pmf[1,2] / (px1*(1-py1))) +
      -1*pmf[2,1]*log2(pmf[2,1] / ((1-px1)*py1)) +
      pmf[2,2]*log2(pmf[2,2] / ((1-px1)*(1-py1)))
    
    jointentropy <- -1*sum(pmf*log2(pmf))
    mutinf <- mutinf / jointentropy
  
    # mutinf <- ifelse(pmf[1,1]==0, 0, pmf[1,1]*log(pmf[1,1] / (px1*py1)))
    # mutinf <- mutinf - ifelse(pmf[1,2]==0, 0, pmf[1,2]*log(pmf[1,2] / (px1*(1-py1))))
    # mutinf <- mutinf - ifelse(pmf[2,1]==0, 0, pmf[2,1]*log(pmf[2,1] / ((1-px1)*py1)))
    # mutinf <- mutinf + ifelse(pmf[2,2]==0, 0, pmf[2,2]*log(pmf[2,2] / ((1-px1)*(1-py1))))
    pval <- 1-fisher.test(conttable)$p.value
    return(mutinf*pval)
  }
  
  pairscores <- BuildSimMatInternal(labvecs, uvals, evalmap, l, n, 
                                    FXN, ARGS, Verbose,
                                    InputIsList=TRUE)
  
  #m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  #pairscores <- pairscores / m
  Diag(pairscores) <- 1
  return(pairscores)
}