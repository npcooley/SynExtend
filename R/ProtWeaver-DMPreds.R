###### Distance Matrix Methods for ProtWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - MirrorTree
#  - ContextTree
##########################

#### S3 Generic Definitions ####
MirrorTree <- function(pw, ...) UseMethod('MirrorTree')
ContextTree <- function(pw, ...) UseMethod('ContextTree')
################################

MirrorTree.ProtWeaver <- function(pw, MTCorrection=c(),
                                  Subset=NULL, Verbose=TRUE,
                                  MySpeciesTree=NULL, 
                                  precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  pl <- length(uvals)
  
  if (is.null(precalcProfs)){
    if (Verbose) cat('Pre-processing distance matrices...\n')
    spl <- NULL
    if (!is.null(MySpeciesTree)) spl <- labels(MySpeciesTree)
    CPs <- CophProfiles(pw, uvals, Verbose=Verbose, speciesList=spl)
  } else {
    CPs <- precalcProfs
  }
  
  l <- ncol(CPs)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- uvals
    return(mat)
  }
  
  MTCorrection <- tolower(MTCorrection)
  if ('speciestree' %in% MTCorrection){
    if (Verbose) cat('Correcting with species tree...\n')
    stopifnot('Missing MySpeciesTree'=!is.null(MySpeciesTree))
    stopifnot('MySpeciesTree must be a dendrogram'=is(MySpeciesTree, 'dendrogram'))
    corrvec <- as.matrix(Cophenetic(MySpeciesTree))
    
    corrvec <- corrvec[upper.tri(corrvec)]
    stopifnot('MySpeciesTree has incorrect number of leaf nodes'=
                length(corrvec)==nrow(CPs))
    CPs <- CPs - corrvec
  }
  if ('normalize' %in% MTCorrection){
    if (Verbose) cat('Normalizing profiles...\n')
    means <- colMeans(CPs, na.rm=TRUE)
    vars <- apply(CPs, MARGIN=2, var, na.rm=TRUE)
    for ( i in seq_len(ncol(CPs)) ){
      CPs[,i] <- (CPs[,i] - means[i]) / (ifelse(vars[i]!=0, sqrt(vars), 1))
    }
  }
  if ('satoaverage' %in% MTCorrection){
    means <- rowMeans(CPs, na.rm = TRUE)
    if (Verbose) cat('Calculating Sato projection vectors...\n')
    
    # Big profiles lead to space issues that crash R
    if (nrow(CPs)**2 < (2**28)){
      if (Verbose) pb <- txtProgressBar(max=ncol(CPs), style=3)
      proj_op <- diag(nrow=nrow(CPs)) - (means %*% t(means))
      for ( i in seq_len(ncol(CPs)) ){
        CPs[,i] <- c(CPs[,i] %*% proj_op)
        if ( Verbose ) setTxtProgressBar(pb, i)
      }
    } else {
      if (Verbose) pb <- txtProgressBar(max=(ncol(CPs)*nrow(CPs)), style=3)
      for ( i in seq_len(ncol(CPs)) ){
        v <- projv <- CPs[,i]
        multv <- means * v
        for ( j in seq_len(nrow(CPs)) ){
          if ( !is.na(v[j]) )
            projv[j] <- sum(multv * means[j])
          if ( Verbose ) setTxtProgressBar(pb, (i-1)*nrow(CPs) + j)
        }
        CPs[,i] <- v - projv
      }
    }
    if (Verbose) cat('\n')
  }
  
  #pairscores <- matrix(NA, nrow=pl, ncol=pl)
  #pairscores <- sim(nelem=pl)
  pairscores <- rep(NA_real_, pl*(pl-1) / 2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(pl*(pl-1) / 2), style=3)
  for ( i in seq_len(pl-1) ){
    uval1 <- uvals[i]
    v1 <- CPs[,i]
    sd1 <- sd(v1, na.rm=TRUE)
    # Should only be NA if there's only one entry
    if (is.na(sd1) || sd1 == 0){
      pairscores[(ctr+1):(ctr+pl)] <- 0
      ctr <- ctr + length((i+1):pl)
      if (Verbose) setTxtProgressBar(pb, ctr)
    } else {
      for ( j in (i+1):pl ){
        uval2 <- uvals[j]
        accessor <- as.character(min(uval1, uval2))
        entry <- max(uval1, uval2)
        if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
          v2 <- CPs[,j]
          sd2 <- sd(v2, na.rm=TRUE)
          if (is.na(sd2) || sd2 == 0)
            val <- 0
          else
            val <- suppressWarnings(cor(v1, v2, 
                                        use='pairwise.complete.obs', 
                                        method='pearson'))
          pairscores[ctr+1] <- ifelse(is.na(val), 0, val)
        }
        ctr <- ctr + 1
        if (Verbose) setTxtProgressBar(pb, ctr)
      }
    }
  }
  if (Verbose) cat('\n')
  pairscores <- as.sim(pairscores, NAMES=names(pw)[uvals], DIAG=FALSE)
  diag(pairscores) <- 1
  if ('partialcorrelation' %in% MTCorrection){
    flag <- TRUE
    pairscores <- as.matrix(pairscores)
    if (!is.null(Subset)){
      opsm <- pairscores
      pairscores <- pairscores[uvals, uvals]
      if (any(is.na(pairscores))){
        pairscores <- opsm
        warning('Partial correlation requires a square matrix. Skipping.')
        flag <- FALSE
      }
    }
    if (flag){
      d <- det(pairscores)
      if (d == 0){
        warning('Matrix is exactly singular, cannot use partial correlation correction.')
      } else {
        inv <- solve(pairscores)
        cols <- matrix(diag(inv), nrow=nrow(inv), ncol=ncol(inv))
        rows <- matrix(diag(inv), nrow=nrow(inv), ncol=ncol(inv), byrow=TRUE)
        divisor <- sqrt(cols * rows)
        pairscores <- (-1 * inv) / divisor
      }
      if ( !is.null(Subset) ){
        opsm[uvals,uvals] <- pairscores
        pairscores <- opsm
      }
    }
    pairscores <- as.sim(pairscores)
  }
  diag(pairscores) <- 1
  
  return(abs(pairscores))
}

ContextTree.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, precalcProfs=NULL, 
                                   MySpeciesTree=NULL, ...){
  
  if ( !is.null(MySpeciesTree) && is(MySpeciesTree, 'dendrogram')){
    MTCorrection <- c('speciestree', 'normalize', 'partialcorrelation')
  } else { 
    MTCorrection <- c('normalize', 
                      'partialcorrelation')
  }
  
  return(MirrorTree(pw, MTCorrection=MTCorrection,
                    Verbose=Verbose, 
                    precalcCProfs=precalcProfs,
                    MySpeciesTree=MySpeciesTree))
}