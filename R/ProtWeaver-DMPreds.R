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
TreeDistance <- function(pw, ...) UseMethod('TreeDistance')
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
    if (Verbose) pb <- txtProgressBar(max=ncol(CPs), style=3)
    for ( i in seq_len(ncol(CPs)) ){
      CPs[,i] <- (CPs[,i] - means[i]) / (ifelse(vars[i]!=0, sqrt(vars), 1))
      if (Verbose) setTxtProgressBar(pb, i)
    }
    if(Verbose) cat('\n')
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
  #pairscores <- simMat(nelem=pl)
  pairscores <- rep(NA_real_, pl*(pl-1) / 2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(pl*(pl-1) / 2), style=3)
  for ( i in seq_len(pl-1) ){
    uval1 <- uvals[i]
    v1 <- CPs[,i]
    sd1 <- sd(v1, na.rm=TRUE)
    # Should only be NA if there's only one entry
    if (is.na(sd1) || sd1 == 0){
      pairscores[(ctr+1):(ctr+pl-i-1)] <- 0
      ctr <- ctr + pl - i
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
  pairscores <- as.simMat(pairscores, NAMES=names(pw)[uvals], DIAG=FALSE)
  Diag(pairscores) <- 1
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
    pairscores <- as.simMat(pairscores)
    Diag(pairscores) <- 1
  }
  
  return(abs(pairscores))
}

ContextTree.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, precalcProfs=NULL, 
                                   MySpeciesTree=NULL, ...){
  
  if ( is.null(MySpeciesTree) || !is(MySpeciesTree, 'dendrogram')){
    MySpeciesTree <- findSpeciesTree(pw, Verbose)
  }
  
  #MTCorrection <- c('speciestree', 'normalize', 'partialcorrelation')
  MTCorrection <- c('speciestree', 'normalize')
  
  return(MirrorTree(pw, MTCorrection=MTCorrection,
                    Verbose=Verbose, 
                    precalcCProfs=precalcProfs,
                    MySpeciesTree=MySpeciesTree))
}

TreeDistance.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE,
                                      precalcSubset=NULL, 
                                      TreeMethods="GRF", JRFk=2, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  
  useColoc <- attr(pw, "useColoc")
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  
  
  bmn <- c("GRF", "RF", "JRF", "Nye", "KF")
  if ('all' %in% TreeMethods){
    bitmask <- rep(TRUE, length(bmn))
  } else {
    bitmask <- rep(FALSE, length(bmn))
    bitmask <- vapply(bmn, \(x) x %in% TreeMethods, logical(1))
  }
  
  bmn <- bmn[bitmask]
  pairscoresList <- vector('list', length=length(bmn))
  for(i in seq_along(bmn))
    pairscoresList[[i]] <- rep(NA_real_, l*(l-1)/2)
  names(pairscoresList) <- bmn
  
  ctr <- 0
  pArray <- vector('list', length=l)
  labelsArray <- vector('list', length=l)
  for ( i in seq_len(l) ){
    tree <- pw[[uvals[i]]]
    if (useColoc){
      tree <- rapply(tree, \(x){
        if (!is.null(attr(x, 'leaf'))){
          attr(x, 'label') <- gsub("(.*)_.*_.*", '\\1', attr(x, 'label'))
        }
        return(x)
      }, how='replace')
    }
    labs <- labels(tree)
    ptr <- .Call("initCDend", tree)
    pArray[[i]] <- ptr
    labelsArray[[i]] <- labs
  }
  on.exit(vapply(seq_len(l), \(x){
    ptr <- pArray[[x]]
    rm(ptr)
    return(0L)
  }, integer(1)))
  
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    p1 <- pArray[[i]]
    uval1 <- uvals[i]
    for ( j in (i+1):l ){
      p2 <- pArray[[j]]
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        interlabs <- intersect(labelsArray[[i]], labelsArray[[j]])
        # GRF
        if (bitmask[1]){
          # GRF
          s <- .Call("GRFInfo", p1, p2, interlabs, FALSE, 0)
          normval <- 0.5*(s[2] + s[3])
          if (is.na(normval) || normval == 0)
            pairscoresList$GRF[ctr+1] <- NA
            #pairscoresList$GRF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$GRF[ctr+1] <- 1 - (normval - s[1]) / normval
        }
        # RF
        if (bitmask[2]){
          s <- .Call("RFDist", p1, p2, interlabs)
          normval <- s[2] + s[3]
          if (is.na(normval) || normval == 0)
            pairscoresList$RF[ctr+1] <- NA
            #pairscoresList$RF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else{
            s <- s[1] / normval
            pairscoresList$RF[ctr+1] <- 1 - s
            # p-value calculation
            #s <- ppois(length(interlabs) - 3 - 0.5*s[1],
            #                                  lambda=1/8, log.p=TRUE, lower.tail=FALSE)
            #pairscoresList$RF[ctr+1] <- 1 - exp(s)
          }
        }
        # JRF
        if (bitmask[3]){
          s <- .Call("GRFInfo", p1, p2, interlabs, TRUE, JRFk)
          normval <- (s[2] + s[3])
          if (is.na(normval) || normval == 0)
            pairscoresList$JRF[ctr+1] <- NA
            #pairscoresList$JRF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$JRF[ctr+1] <- 1 - (s[1] / normval)
        }
        # Nye
        if (bitmask[4]){
          s <- .Call("GRFInfo", p1, p2, interlabs, TRUE, 1)
          normval <- (s[2] + s[3])
          if (is.na(normval) || normval == 0)
            pairscoresList$Nye[ctr+1] <- NA
            #pairscoresList$Nye[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$Nye[ctr+1] <- 1 - (s[1] / normval)
        }
        # KF
        if (bitmask[5]){
          s <- .Call("KFDist", p1, p2, interlabs)
          normval <- s[2]
          if (is.na(normval) || normval == 0)
            pairscoresList$KF[ctr+1] <- NA
            #pairscoresList$KF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$KF[ctr+1] <- s[1] / normval
        }
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  
  if (bitmask[5]){
    pairscoresList$KF <- pairscoresList$KF - min(pairscoresList$KF, na.rm=TRUE)
    pairscoresList$KF <- pairscoresList$KF / max(pairscoresList$KF, na.rm=TRUE)
    pairscoresList$KF <- 1 - pairscoresList$KF
  }
  n <- n[uvals]
  for (i in seq_along(pairscoresList)){
    pairscoresList[[i]] <- as.simMat(pairscoresList[[i]], NAMES=n, DIAG=FALSE)
  }

  return(pairscoresList) 
}