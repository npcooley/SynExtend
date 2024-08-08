###### Presence/Absence Methods for EvoWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu
# TODO: clean up the naming here

#### Implemented Methods: ####
#  - Jaccard Distance at tips
#  - Hamming Distance
#  - Mutual Information at tips
#  - Mutual Information of gain/loss events
#  - Jaccard Distance, centered w/ conserved clades collapsed
#  - Inverse Potts Model
#  - Behdenna 2016 gain/loss
#  - Pearson correlation of gain/loss (w/ p-value)
#  - Distance between gain/loss events
#  - Overlap of ancestral state coverage on tree
##########################

#### S3 Generic Definitions ####
ExtantJaccard <- function(ew, ...) UseMethod('ExtantJaccard')
Hamming <- function(ew, ...) UseMethod('Hamming')
CorrGL <- function(ew, ...) UseMethod('CorrGL')
GLMI <- function(ew, ...) UseMethod('GLMI')
ProfileDCA <- function(ew, ...) UseMethod('ProfileDCA')
Behdenna <- function(ew, ...) UseMethod('Behdenna')
GLDistance <- function(ew, ...) UseMethod('GLDistance')
PAJaccard <- function(ew, ...) UseMethod("PAJaccard")
PAOverlap <- function(ew, ...) UseMethod("PAOverlap")
PAPV <- function(ew, ...) UseMethod('PAPV')
################################

ExtantJaccard.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                               precalcProfs=NULL, precalcSubset=NULL, ...){

  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  nr <- nrow(pap)
  pap[] <- as.integer(pap)
  ARGS <- list(nr=nr)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    return(.Call("calcScoreJaccard", v1, v2, ARGS$nr, PACKAGE="SynExtend"))
  }
  pairscores <- BuildSimMatInternal(pap, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}

Hamming.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- simMat(1, nelem=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  nr <- nrow(pap)
  pairscores <- rep(NA_real_, l*(l-1)/2)
  #nc <- ncol(pap)
  pap[] <- as.integer(pap)
  ARGS <- list(nr=nr)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    return(.Call("calcScoreHamming", v1, v2, ARGS$nr, 1, PACKAGE="SynExtend"))
  }
  pairscores <- BuildSimMatInternal(pap, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}

PAJaccard.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                    MySpeciesTree=NULL,
                                    precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if (is.null(MySpeciesTree)){
    stopifnot("Method 'GLDistance' requires a species tree"=attr(ew, 'useMT'))
    if (Verbose) cat('Calculating Species Tree...\n')
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree), ...)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- simMat(1, nelem=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  # Initialize Dendrogram in C
  y <- .Call("initCDend", MySpeciesTree, PACKAGE="SynExtend")
  on.exit(rm(y))
  numnodes <- .Call("getTreeNodesCount", y, PACKAGE="SynExtend")
  glvs <- matrix(NA_integer_, nrow=numnodes, ncol=l)
  rn <- rownames(pap)

  # Calculate Gain/Loss Vectors
  if (Verbose) cat('  Calculating ancestral states:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for (i in seq_len(l)){
    v <- rn[pap[,i]]
    if (length(v) == 0){
      glv <- rep(0L, numnodes)
    } else {
      glv <- .Call("calcGainLoss", y, v, FALSE, PACKAGE="SynExtend")
    }
    glvs[,i] <- glv
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat("\n")

  ARGS <- list(numnodes=numnodes, y=y)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    if(all(v1[1]==v1) || all(v2[1]==v2)){
      # catch case where standard deviation is 0
      # (all elements of vector the same)
      return(0)
    }
    y <- ARGS$y
    # 00 = 1, 01 = 2, 10 = 3, 11 = 4
    # 0 reserved for "ignore me"
    v1 <- 2L*v1 + v2 + 1L
    v1 <- .Call("cladeCollapsePA", y, v1)
    jpd <- tabulate(v1, nbins=4L)
    ## now jpd is: 00 01 10 11
    #pval <- 1-fisher.test(jpd, alternative='two.sided')$p.value
    #jpd <- jpd / sum(jpd)

    # now do we use Jaccard, or MI?
    # Jaccard does have a nice p-value we could use
    # see Chung et al. 2019
    # note indices are reversed relative to paper
    m <- sum(jpd)
    jaccard_index <- ifelse(sum(jpd[2:4]) == 0, 0, jpd[4] / sum(jpd[2:4]))
    p_i <- sum(jpd[3:4]) / m
    p_j <- sum(jpd[c(2,4)]) / m
    q1 <- p_i * p_j
    q2 <- p_i + p_j - 2*p_i*p_j
    expected_jaccard <- ifelse(q2+q1 == 0, 0, q1 / (q2+q1))
    centered_jaccard <- jaccard_index - expected_jaccard

    # Bootstrap p-value calculation
    # note that this is better than asymptotics,
    # which are anti-conservative estimators
    NUM_BOOTSTRAP <- min(m*5L, 1000L)
    NUM_BOOTSTRAP <- max(NUM_BOOTSTRAP, 5L)
    bootstrap_replicates <- numeric(NUM_BOOTSTRAP)
    for(iter in seq_along(bootstrap_replicates)){
      s <- (rbinom(m, 1, p_i) * 2) + (rbinom(m, 1, p_j) + 1L)
      s <- tabulate(s, nbins=4)
      jv <- s[4] / sum(s[2:4])
      if(is.nan(jv) || is.na(jv)) jv <- 0
      bootstrap_replicates[iter] <- jv - expected_jaccard
    }

    p_value_left <- sum(bootstrap_replicates >= centered_jaccard) / NUM_BOOTSTRAP
    p_value_right <- sum(bootstrap_replicates <= centered_jaccard) / NUM_BOOTSTRAP
    p_value <- 2*min(p_value_left, p_value_right)

    weighting <- 1-sqrt(p_value)
    # the centered jaccard is bounded on [-1/3,2/3] due to how it's calculated
    if(centered_jaccard < 0) centered_jaccard <- centered_jaccard * 3
    else centered_jaccard <- centered_jaccard * 3/2

    # clip to [-1,1] just in case
    if(centered_jaccard > 1) centered_jaccard <- 1
    if(centered_jaccard < -1) centered_jaccard <- -1

    return(weighting * centered_jaccard)
  }
  pairscores <- BuildSimMatInternal(glvs, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}


PAOverlap.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                MySpeciesTree=NULL,
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  N_PERM <- 200L
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if (is.null(MySpeciesTree)){
    stopifnot("Method 'PAOverlap' requires a species tree"=attr(ew, 'useMT'))
    if (Verbose) cat('Calculating Species Tree...\n')
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree), ...)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)

  if ( l == 1 ){
    mat <- simMat(1, nelem=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  # Initialize Dendrogram in C
  y <- .Call("initCDend", MySpeciesTree, PACKAGE="SynExtend")
  on.exit(rm(y))
  numnodes <- .Call("getTreeNodesCount", y, PACKAGE="SynExtend")
  rn <- rownames(pap)
  glvs <- a_states <- matrix(NA_integer_, nrow=numnodes, ncol=l)

  # Calculate Gain/Loss Vectors
  if (Verbose) cat('  Calculating gain/loss vectors:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for (i in seq_len(l)){
    v <- rn[pap[,i]]
    if (length(v) == 0){
      glv <- a_state <- rep(0L, numnodes)
    } else {
      glv <- .Call("calcGainLoss", y, v, TRUE, PACKAGE="SynExtend")
      a_state <- .Call("calcGainLoss", y, v, FALSE, PACKAGE='SynExtend')
    }
    glvs[,i] <- glv
    a_states[,i] <- a_state
    if (Verbose) setTxtProgressBar(pb, i)
  }
  glvs <- rbind(glvs, a_states)
  if(Verbose) cat("\n")

  all_lengths <- .Call("calcAllTreeLengths", y, PACKAGE='SynExtend')

  ARGS <- list(numnodes=numnodes, all_lens=all_lengths,
               total_len=sum(all_lengths), n_perm=N_PERM)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    nn <- ARGS$numnodes
    pa1 <- v1[seq_len(nn)+nn]
    pa2 <- v2[seq_len(nn)+nn]
    if(sum(pa1) * sum(pa2) == 0) return(0)
    glv1 <- v1[seq_len(nn)]
    glv2 <- v2[seq_len(nn)]

    total_len <- ARGS$total_len
    lens <- ARGS$all_lens
    n_perm <- ARGS$n_perm

    # ignore length 0 edges unless there's an event on one
    pos_real <- which(lens > 0 | (lens==0 & abs(glv1)+abs(glv2) != 0))
    pa1 <- pa1[pos_real]
    pa2 <- pa2[pos_real]
    glv1 <- glv1[pos_real]
    glv2 <- glv2[pos_real]
    lens <- lens[pos_real]
    numnodes <- length(pa1)

    # scoring
    time_len1 <- sum((pa1 - (0.5*abs(glv1))) * lens)
    time_len2 <- sum((pa2 - (0.5*abs(glv2))) * lens)
    denom <- time_len1 + time_len2
    if(denom == 0) denom <- Inf

    #weighting <- c(time_len1 / total_len, time_len2 / total_len)
    #weighting <- ifelse(any(weighting==0), 0, 2 / sum(1/sqrt(weighting)))
    pos_inter <- which(pa1*pa2==1 & glv1*glv2 != -1)
    # subsetting here means we won't have simultaneous G/L or L/G, only G/G and L/L
    inter_lens <- lens[pos_inter]
    p <- abs(glv1[pos_inter]) + abs(glv2[pos_inter])!=0
    inter_lens[p] <- inter_lens[p] * 0.5
    inter_lens <- 2*sum(inter_lens) / denom
    score <- inter_lens - (time_len1 * time_len2) / (total_len*total_len)
    ## commented lines change it from permutation testing to bootstrapping
    #probs_pa <- c(sum(pa1) / nn, sum(pa2)/nn)
    #probs_glv1 <- tabulate(glv1+2L, nbins=3) / nn
    #probs_glv2 <- tabulate(glv2+2L, nbins=3) / nn
    replicates <- numeric(n_perm)
    for(iter in seq_len(n_perm)){
      pa1 <- .C("shuffleRInt", pa1, numnodes, PACKAGE="SynExtend")[[1]]
      pa2 <- .C("shuffleRInt", pa2, numnodes, PACKAGE="SynExtend")[[1]]
      glv1 <- .C("shuffleRInt", glv1, numnodes, PACKAGE="SynExtend")[[1]]
      glv2 <- .C("shuffleRInt", glv2, numnodes, PACKAGE="SynExtend")[[1]]
      # pa1 <- rbinom(nn, 1, probs_pa[1])
      # pa2 <- rbinom(nn, 1, probs_pa[2])
      # glv1 <- sample(c(-1:1), nn, replace=TRUE, probs_glv1)
      # glv2 <- sample(c(-1:1), nn, replace=TRUE, probs_glv2)
      time_len1 <- sum((pa1 - (0.5*abs(glv1))) * lens)
      time_len2 <- sum((pa2 - (0.5*abs(glv2))) * lens)
      denom <- time_len1 + time_len2
      if(denom == 0) denom <- Inf
      pos_inter <- which((pa1-0.5*glv1) == (pa2-0.5*glv2) & pa1*pa2==1)
      inter_lens <- lens[pos_inter]
      p <- abs(glv1[pos_inter]) + abs(glv2[pos_inter])!=0
      inter_lens[p] <- inter_lens[p] * 0.5
      inter_lens <- 2*sum(inter_lens) / denom
      replicates[iter] <- abs(inter_lens - (time_len1 * time_len2) / (total_len*total_len))
    }
    p_lv <- sum(replicates <= score) / n_perm
    p_rv <- sum(replicates >= score) / n_perm
    p_v <- 2*min(c(p_lv, p_rv))
    if(p_v > 1) p_v <- 1

    return(sqrt(1-p_v)*score)
    #return(sqrt(1-p_v)*score*weighting)
  }
  pairscores <- BuildSimMatInternal(glvs, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}

CorrGL.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                MySpeciesTree=NULL,
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if (is.null(MySpeciesTree)){
    stopifnot("Method 'GLDistance' requires a species tree"=attr(ew, 'useMT'))
    if (Verbose) cat('Calculating Species Tree...\n')
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree), ...)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)

  if ( l == 1 ){
    mat <- simMat(1, nelem=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  # Initialize Dendrogram in C
  y <- .Call("initCDend", MySpeciesTree, PACKAGE="SynExtend")
  on.exit(rm(y))
  numnodes <- .Call("getTreeNodesCount", y, PACKAGE="SynExtend")
  rn <- rownames(pap)
  glvs <- matrix(NA_integer_, nrow=numnodes, ncol=l)

  # Calculate Gain/Loss Vectors
  if (Verbose) cat('  Calculating gain/loss vectors:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for (i in seq_len(l)){
    v <- rn[pap[,i]]
    if (length(v) == 0){
      glv <- rep(0L, numnodes)
    } else {
      glv <- .Call("calcGainLoss", y, v, TRUE, PACKAGE="SynExtend")
    }
    glvs[,i] <- glv
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat("\n")

  ARGS <- list(numnodes=numnodes)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    if(all(v1[1]==v1) || all(v2[1]==v2)){
      # catch case where standard deviation is 0
      # (all elements of vector the same)
      return(0)
    }
    val <- cor(v1, v2)
    #pval <- 1 - pt(val, ARGS$numnodes - 2, lower.tail=FALSE)
    pval <- 1-fisher.test(v1, v2, simulate.p.value=TRUE)$p.value
    return(pval*val)
  }
  pairscores <- BuildSimMatInternal(glvs, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}


MutualInformationPA.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                         precalcProfs=NULL, precalcSubset=NULL,
                                        switchval=0L, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  FXN <- function(v1, v2, ARGS, ii, jj){
    score <- 0
    v1l <- length(v1)
    if(v1l == 0) return(0)
    tt <- sum(v1 & v2) / v1l
    tf <- sum(v1 & !v2) / v1l
    ft <- sum(!v1 & v2) / v1l
    ff <- sum(!v1 & !v2) / v1l
    jpd <- c(tt, tf, ft, ff)

    tv1 <- sum(v1) / v1l
    tv2 <- sum(v2) / v1l
    fv1 <- 1 - tv1
    fv2 <- 1 - tv2
    mpdv1 <- c(tv1, tv1, fv1, fv1)
    mpdv2 <- c(tv2, fv2, tv2, fv2)

    mult <- c(1,-1,-1,1)

    for ( k in seq_along(jpd) ){
      val <- jpd[k] * log(jpd[k] / (mpdv1[k] * mpdv2[k]), base=2) * mult[k]
      score <- score + ifelse(is.nan(val), 0, val)
    }
    pv <- ifelse(all(v1==v1[1]) || all(v2==v2[1]), 1, fisher.test(v1,v2)$p.value)
    return(score*(1-pv))
  }

  # APC Correction, removing for now
  # CORRECTION <- function(ps){
  #   apccorr <- mean(ps, na.rm=TRUE)
  #   ps <- ps - apccorr
  #   ps <- abs(ps)
  #   # Normalize
  #   denom <- max(ps, na.rm=TRUE)
  #   ps <- ps / ifelse(denom==0, 1, denom)
  # }
  pairscores <- BuildSimMatInternal(pap, uvals, evalmap, l, n,
                                    FXN, NULL, Verbose)
                                    #CORRECTION=CORRECTION)

  Diag(pairscores) <- 1
  #pairscores <- pairscores #because distance
  return(pairscores)
}

GLMI.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                         precalcProfs=NULL, precalcSubset=NULL,
                                        MySpeciesTree=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if (is.null(MySpeciesTree)){
    stopifnot("Method 'GLMI' requires a species tree"=attr(ew, 'useMT'))
    if (Verbose) cat('Calculating Species Tree...\n')
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  y <- .Call("initCDend", MySpeciesTree, PACKAGE="SynExtend")
  on.exit(rm(y))
  numnodes <- .Call("getTreeNodesCount", y, PACKAGE="SynExtend")
  rn <- rownames(pap)
  glvs <- matrix(NA_integer_, nrow=numnodes, ncol=l)

  # Calculate Gain/Loss Vectors
  if (Verbose) cat('  Calculating gain/loss vectors:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for (i in seq_len(l)){
    v <- rn[pap[,i]]
    if (length(v) == 0){
      glv <- rep(0L, numnodes)
    } else {
      glv <- .Call("calcGainLoss", y, v, TRUE, PACKAGE="SynExtend")
    }
    glvs[,i] <- glv
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat("\n")

  FXN <- function(v1, v2, ARGS, ii, jj){
    # subset just to where events actually occur
    p <- abs(v1) + abs(v2) != 0
    v1 <- v1[p]
    v2 <- v2[p]
    v1p <- v2p <- rep(0, length(v1))
    v1p[v1==v2] <- 1
    v2p[v1==v2] <- 1
    v1p[v1==1 & v2==-1] <- 1
    v2p[v2==1 & v1==-1] <- 1
    v1 <- v1p
    v2 <- v2p
    if(all(v1==v1[1]) || all(v2==v2[1])) return(0)

    score <- 0
    v1l <- length(v1)
    if(v1l == 0) return(0)
    tt <- sum(v1 & v2) / v1l
    tf <- sum(v1 & !v2) / v1l
    ft <- sum(!v1 & v2) / v1l
    ff <- sum(!v1 & !v2) / v1l
    jpd <- c(tt, tf, ft, ff)

    tv1 <- sum(v1) / v1l
    tv2 <- sum(v2) / v1l
    fv1 <- 1 - tv1
    fv2 <- 1 - tv2
    mpdv1 <- c(tv1, tv1, fv1, fv1)
    mpdv2 <- c(tv2, fv2, tv2, fv2)

    # +1 for same event, -1 for opposite, 0 for event + no change
    mult <- c(1,-1,-1,0)

    for ( k in seq_along(jpd) ){
      val <- jpd[k] * log(jpd[k] / (mpdv1[k] * mpdv2[k]), base=2) * mult[k]
      score <- score + ifelse(is.nan(val), 0, val)
    }
    pv <- fisher.test(v1,v2)$p.value
    return(score*(1-pv))
  }

  # APC Correction, removing for now
  # CORRECTION <- function(ps){
  #   apccorr <- mean(ps, na.rm=TRUE)
  #   ps <- ps - apccorr
  #   ps <- abs(ps)
  #   # Normalize
  #   denom <- max(ps, na.rm=TRUE)
  #   ps <- ps / ifelse(denom==0, 1, denom)
  # }
  pairscores <- BuildSimMatInternal(glvs, uvals, evalmap, l, n,
                                    FXN, NULL, Verbose)
                                    #CORRECTION=CORRECTION)

  Diag(pairscores) <- 1
  #pairscores <- pairscores #because distance
  return(pairscores)
}

ProfileDCA.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE, Processors=1L,
                                  precalcProfs=NULL, precalcSubset=NULL, useAbs=TRUE, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- ncol(pap)
  n <- colnames(pap)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }

  pairscores <- DCA_logRISE(pap, Verbose=Verbose, Processors=Processors, ...)
  rownames(pairscores) <- colnames(pairscores) <- n
  pairscores <- as.simMat(pairscores)
  if (useAbs) pairscores <- abs(pairscores)
  if (max(pairscores) != 0)
    pairscores <- pairscores / max(abs(pairscores))

  return(pairscores)
}


Behdenna.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                MySpeciesTree=NULL, useSubtree=FALSE,
                                useACCTRAN=TRUE, rawZScores=FALSE,
                                precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap

  if (is.null(MySpeciesTree)){
      stopifnot("Method 'Behdenna' requires a species tree"=attr(ew, 'useMT'))
      if (Verbose) cat('Calculating Species Tree...\n')
      MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree))
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  stopifnot(nrow(pap) > 1)
  fd <- FastDend(MySpeciesTree)
  v1 <- abs(generateGainLossVec(fd, pap[,1], moveEventsUpward=useACCTRAN))
  glmat <- matrix(0, nrow=length(v1), ncol=ncol(pap))
  glmat[,1] <- v1
  if (Verbose) cat('  Calculating gain/loss vectors:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for ( i in 2:ncol(pap) ){
    glmat[,i] <- abs(generateGainLossVec(fd, pap[,i], moveEventsUpward=useACCTRAN))
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if (Verbose) cat('\n')
  vals <- calc_SId_mat(fd, IdOnly=!useSubtree)
  if (useSubtree)
    M <- vals$S
  else
    M <- vals$Id
  Cmat <- vals$Cmat
  bl <- vals$blengths

  glmat <- abs(glmat)

  ARGS <- list(M=M, Cmat=Cmat, bl=bl)
  FXN <- function(v1, v2, ARGS, ii, jj){
    n1 <- sum(v1)
    n2 <- sum(v2)
    score <- 0
    if ( n1*n2 != 0 ){
      score <- t(v1) %*% ARGS$M %*% v2
      exp_mean <- n2 * (t(v1) %*% ARGS$M %*% ARGS$bl)
      exp_var <- n2*t(v1) %*% ARGS$M %*% ARGS$Cmat %*% t(M) %*% v1
      score <- (score - exp_mean) / sqrt(abs(exp_var))
    }

    return(score)
  }
  CORRECTION <- NULL
  if (!rawZScores){
    CORRECTION <- function(ps){
      ps <- abs(ps)
      ps <- ps / ifelse(max(ps,na.rm=TRUE) != 0,
                        max(ps, na.rm=TRUE), 1)
      return(ps)
    }
  }
  pairscores <- BuildSimMatInternal(glmat, uvals, evalmap, l, names(ew),
                                    FXN, ARGS, Verbose,
                                    CORRECTION=CORRECTION)

  Diag(pairscores) <- ifelse(rawZScores, 1, 0)

  return(pairscores)
}

GLDistance.EvoWeaver <- function(ew, Subset=NULL,
                     Verbose=TRUE, MySpeciesTree=NULL,
                     precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  BOOTSTRAP_NUM <- 50L

  if (is.null(MySpeciesTree)){
    stopifnot("Method 'GLDistance' requires a species tree"=attr(ew, 'useMT'))
    if (Verbose) cat('Calculating Species Tree...\n')
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree), ...)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  stopifnot(nrow(pap) > 1)
  pa_counts <- colSums(pap)

  # Initialize Dendrogram in C
  y <- .Call("initCDend", MySpeciesTree, PACKAGE="SynExtend")
  on.exit(rm(y))
  numnodes <- .Call("getTreeNodesCount", y, PACKAGE="SynExtend")
  rn <- rownames(pap)
  glvs <- matrix(NA_integer_, nrow=numnodes, ncol=l)
  allnonzero <- logical(l)
  # Calculate Gain/Loss Vectors
  if (Verbose) cat('  Calculating gain/loss vectors:\n')
  if (Verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for (i in seq_len(l)){
    v <- rn[pap[,i]]
    if (length(v) == 0){
      glv <- rep(0L, numnodes)
      allnonzero[i] <- TRUE
    } else {
      glv <- .Call("calcGainLoss", y, v, TRUE, PACKAGE="SynExtend")
    }
    glvs[,i] <- glv
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat("\n")
  lenglv <- nrow(glvs)
  ARGS <- list(allnonzero=allnonzero, y=y, bsn=BOOTSTRAP_NUM, l=lenglv)
  FXN <- function(v1, v2, ARGS, ii, jj){
    allnonzero <- ARGS$allnonzero
    #labs <- ARGS$labs
    #n1 <- ARGS$counts[ii]
    #n2 <- ARGS$counts[jj]
    if (allnonzero[ii] || allnonzero[jj]){
      return(0)
    } else {
      #res1 <- .Call("calcScoreGL", ARGS$y, v1, v2, PACKAGE="SynExtend")
      #res2 <- .Call("calcScoreGL", ARGS$y, v2, v1, PACKAGE="SynExtend")
      res1 <- .Call("calcScoreGL", ARGS$y, v1, v2, PACKAGE="SynExtend") / sum(abs(v1))
      res2 <- .Call("calcScoreGL", ARGS$y, v2, v1, PACKAGE="SynExtend") / sum(abs(v2))
      res <- res1 + res2
      if (is.na(res) || is.infinite(res) || res==0) return(0)
      num_bs <- ARGS$bsn
      if(num_bs > 0){
        v1 <- as.integer(v1)
        v2 <- as.integer(v2)
        l <- as.integer(ARGS$l)
        replicates <- rep(NA_real_, num_bs)
        #ss <- seq_along(labs)
        #ssl <- length(ss)
        for(i in seq_len(num_bs)){
          v1tmp <- .C("shuffleRInt", v1, l, PACKAGE="SynExtend")[[1]]
          #pa1 <- labs[v1tmp]
          v2tmp <- .C("shuffleRInt", v2, l, PACKAGE="SynExtend")[[1]]
          #pa2 <- labs[v2tmp]
          #v1tmp <- .Call("calcGainLoss", y, pa1, TRUE, PACKAGE="SynExtend")
          #v2tmp <- .Call("calcGainLoss", y, pa2, TRUE, PACKAGE="SynExtend")
          replicates[i] <- .Call("calcScoreGL", ARGS$y, v1tmp, v2tmp, PACKAGE="SynExtend") / sum(abs(v1tmp)) +
                            .Call("calcScoreGL", ARGS$y, v2tmp, v1tmp, PACKAGE="SynExtend") / sum(abs(v2tmp))
        }

        # pair up gain/loss and loss/gain
        #replicates <- replicates / normer
        #res <- res / normer

        # doing abs instead so we get a 1-sided pv
        # also fixes cases where all values are equal
        pv <- sum(abs(res) > abs(replicates)) / num_bs
        # left.pv <- sum(res <= replicates) / num_bs
        # pv <- 2 * min(right.pv, left.pv)
        # pv <- sum(replicates < res) / num_bs
        # if(pv > 0.5){
        #   pv <- 2*(1-pv)
        # } else {
        #   pv <- 2*pv
        # }
        res <- (pv)*(res/2)
      } else{
        normer <- mean(sum(abs(v1)), sum(abs(v2)))
        normer <- sum(abs(v1), abs(v2))
        normer <- ifelse(normer==0, 1, normer)
        res <- 2*res / normer
      }
      return(res)
    }
  }
  pairscores <- BuildSimMatInternal(glvs, uvals, evalmap, l, names(ew),
                                    FXN, ARGS, Verbose)

  return(pairscores)
}

PAPV.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                               precalcProfs=NULL, precalcSubset=NULL, ...){

  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('  Calculating PA Profiles...\n')
    pap <- PAProfiles(ew, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  nr <- nrow(pap)
  pap[] <- as.integer(pap)
  ARGS <- list(nr=nr)
  FXN <- function(v1, v2, ARGS, ii, jj) {
    if(all(v1==v1[1]) || all(v2==v2[1])) return(0)
    return(1-fisher.test(v1, v2, simulate.p.value=TRUE)$p.value)
  }
  pairscores <- BuildSimMatInternal(pap, uvals, evalmap, l, n, FXN, ARGS, Verbose)

  return(pairscores)
}
