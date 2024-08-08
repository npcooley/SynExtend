###### Co-localization Methods for EvoWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - Naive GeneDistance
##########################

#### S3 Generic Definitions ####
GeneDistance <- function(ew, ...) UseMethod('GeneDistance')
MoransI <- function(ew, ...) UseMethod('MoransI')
OrientationMI <- function(ew, ...) UseMethod('OrientationMI')
################################

GeneDistance.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                             precalcProfs=NULL, precalcSubset=NULL,
                             minimumGenomeSize=2500, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  stopifnot('Colocalization is disabled.'=attr(ew,'useColoc'))
  if (attr(ew, 'useMT')){
    labvecs <- lapply(ew[uvals], labels)
  } else {
    labvecs <- ew[uvals]
  }
  l <- length(labvecs)

  allsp <- lapply(labvecs, \(x) gsub('([^_]+)_.*$', '\\1', x))
  n <- names(ew)

  ARGS <- list(allsp=allsp)
  FXN <- function(lab1, lab2, ARGS, ii, jj){
    l1sp <- gsub('([^_]+)_.*$', '\\1', lab1)
    l2sp <- gsub('([^_]+)_.*$', '\\1', lab2)
    l1chrom <- gsub("([^_]+_[^_]+)_.*", '\\1', lab1)
    l2chrom <- gsub("([^_]+_[^_]+)_.*", '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    if(length(shared) == 0) return(0)
    score <- numeric(length(shared))
    for ( k in seq_along(shared) ){
      spk <- shared[k]
      match1 <- match(spk, l1sp)
      match2 <- match(spk, l2sp)
      vec1 <- lab1[match1]
      chrom1 <- l1chrom[match1]
      vec2 <- lab2[match2]
      chrom2 <- l2chrom[match2]
      idx1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec1, useBytes=TRUE))
      idx2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', vec2, useBytes=TRUE))
      # Double loop is SIGNIFICANTLY faster than expand.grid
      # ...even though it looks so gross :(
      diffs <- Inf
      for(i in seq_along(idx1)){
        for(j in seq_along(idx2)){
          if(chrom1[i] == chrom2[j]){
            diffs <- min(diffs, abs(idx1[i] - idx2[j]))
          }
        }
      }
      if(diffs==0) diffs <- 1
      score[k] <- ifelse(is.infinite(diffs), 0, exp(1-diffs))
    }
    pval <- 0
    if(length(score) > 0){
      rawscore <- -1*log(score)
      rawscore[is.infinite(rawscore)] <- minimumGenomeSize
      pvals <- (1-(rawscore / minimumGenomeSize))**2
      pvals <- pmax(pvals, 0)
      score <- mean(score * pvals)
    } else {
      score <- 0
    }
    return(score)
  }

  pairscores <- BuildSimMatInternal(labvecs, uvals, evalmap, l, n,
                                    FXN, ARGS, Verbose,
                                    InputIsList=TRUE)

  m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  pairscores <- pairscores / m
  Diag(pairscores) <- 1
  return(pairscores)
}

MoransI.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                  MySpeciesTree=NULL,
                                  precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  stopifnot('Colocalization is disabled.'=attr(ew,'useColoc'))
  if ( is.null(MySpeciesTree) || !is(MySpeciesTree, 'dendrogram')){
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }
  stopifnot('Missing MySpeciesTree'=!is.null(MySpeciesTree))
  stopifnot('MySpeciesTree must be a dendrogram'=is(MySpeciesTree, 'dendrogram'))
  if (attr(ew, 'useMT')){
    labvecs <- lapply(ew[uvals], labels)
  } else {
    labvecs <- ew[uvals]
  }
  l <- length(labvecs)
  specCoph <- as.matrix(fastCoph(MySpeciesTree))

  allsp <- lapply(labvecs, \(x) gsub('(.*)_.*$', '\\1', x))
  n <- names(ew)

  ARGS <- list(allsp=allsp, cMat=specCoph)
  FXN <- function(lab1, lab2, ARGS, ii, jj){
    l1sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab1)
    l2sp <- gsub('([^_]*_[0-9]*)_.*$', '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    if (length(shared) <= 3){
      return(0)
    }
    score <- 0
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
    res <- MoranI(vals, w, alternative = 'greater')
    #score <- mean(vals)
    score <- res$observed - res$expected
    pval <- 1-res$p.value
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

OrientationMI.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                          precalcProfs=NULL, precalcSubset=NULL, ...){
  stopifnot('Some labels are missing strand identifiers!'=attr(ew, 'useStrand'))
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  stopifnot('Colocalization is disabled.'=attr(ew,'useColoc'))
  if (attr(ew, 'useMT')){
    labvecs <- lapply(ew[uvals], labels)
  } else {
    labvecs <- ew[uvals]
  }
  l <- length(labvecs)

  allsp <- lapply(labvecs, \(x) gsub('([^_]+)_.*$', '\\1', x))
  n <- names(ew)

  ARGS <- list(allsp=allsp)
  FXN <- function(lab1, lab2, ARGS, ii, jj){
    l1sp <- gsub('([^_]*_[^_]*)_.*$', '\\1', lab1)
    l2sp <- gsub('([^_]*_[^_]*)_.*$', '\\1', lab2)
    shared <- intersect(l1sp, l2sp)
    if(length(shared) == 0) return(0)
    score <- 0
    vals <- rep(NA_real_, length(shared))

    ## a pseudocount of 1 is added
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
    mutinf <- 0
    vx1 <- rep(c(px1, 1-px1), times=2)
    vy1 <- rep(c(py1, 1-py1), 2)
    jointent <- 0
    for(i in seq_len(4L)){
      if(pmf[i] != 0){
        mutinf <- mutinf +
          ifelse(i%in%2:3,-1,1)*pmf[i]*log2(pmf[i] / (vx1[i]*vy1[i]))
        jointent <- jointent + pmf[i]*log2(pmf[i])
      }
    }
    # mutinf <- pmf[1,1]*log2(pmf[1,1] / (px1*py1)) +
    #   -1*pmf[1,2]*log2(pmf[1,2] / (px1*(1-py1))) +
    #   -1*pmf[2,1]*log2(pmf[2,1] / ((1-px1)*py1)) +
    #   pmf[2,2]*log2(pmf[2,2] / ((1-px1)*(1-py1)))

    #jointentropy <- -1*sum(pmf*log2(pmf))
    mutinf <- ifelse(jointent==0, mutinf, mutinf / jointent)

    pval <- 1-fisher.test(conttable)$p.value
    return(abs(mutinf)*pval)
  }

  pairscores <- BuildSimMatInternal(labvecs, uvals, evalmap, l, n,
                                    FXN, ARGS, Verbose,
                                    InputIsList=TRUE)

  #m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  #pairscores <- pairscores / m
  Diag(pairscores) <- 1
  return(pairscores)
}
