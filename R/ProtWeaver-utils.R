#### Helper Functions for ProtWeaver class ####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### S3 Generic Definitions ####
PAProfiles <- function(pw, ...) UseMethod('PAProfiles')
CophProfiles <- function(pw, ...) UseMethod('CophProfiles')
RandCophProfiles <- function(pw, ...) UseMethod('RandCophProfiles')
################################

BuildSimMatInternal <- function(vecs, uvals, evalmap, l, n, FXN, ARGS, Verbose,
                                CORRECTION=NULL, InputIsList=FALSE){
  pairscores <- rep(NA_real_, l*(l-1)/2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    if (InputIsList){
      v1 <- vecs[[i]]
    } else {
      v1 <- vecs[,i]
    }
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        if (InputIsList){
          v2 <- vecs[[j]]
        } else {
          v2 <- vecs[,j]
        }
        pairscores[ctr+1] <- FXN(v1, v2, ARGS, i, j)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  
  n <- n[uvals]
  if (!is.null(CORRECTION)){
    pairscores <- CORRECTION(pairscores)
  }
  pairscores <- as.simMat(pairscores, NAMES=n, DIAG=FALSE)
  return(pairscores)
}

PAProfiles.ProtWeaver <- function(pw, toEval=NULL, Verbose=TRUE, 
                                  speciesList=NULL, ...){
  cols <- names(pw)
  ao <- attr(pw, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(pw, 'useColoc')
  useMT <- attr(pw, 'useMT')
  if (useMT)
    pw <- lapply(pw, labels)
  if (useColoc)
    pw <- lapply(pw, gsub, pattern='(.+)_.+_[0-9]+', replacement='\\1')
  
  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  profiles <- matrix(FALSE, nrow=length(allOrgs), ncol=length(pw))
  rownames(profiles) <- allOrgs
  colnames(profiles) <- cols
  if (Verbose) pb <- txtProgressBar(max=length(pw), style=3)
  for ( i in seq_len(length(pw)) ){
    if( !skip || i %in% locs)
      profiles[pw[[i]],i] <- TRUE
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if (Verbose) cat('\n')
  if (!is.null(toEval))
    profiles <- profiles[,locs]
  return(profiles)
}

CophProfiles.ProtWeaver <- function(pw, toEval=NULL, Verbose=TRUE, 
                                    speciesList=NULL, ...){
  ## TODO: Some way to handle paralogs
  cols <- names(pw)
  ao <- attr(pw, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(pw, 'useColoc')
  useMT <- attr(pw, 'useMT')
  
  stopifnot('ProtWeaver object must be initialized with dendrograms to run MirrorTree methods'=
              useMT)
  
  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  l <- length(allOrgs)
  num_entries <- (l * (l-1)) / 2
  outmat <- matrix(0, nrow=num_entries, ncol=length(pw))
  dummycoph <- matrix(NA, nrow=l, ncol=l)
  ut <- upper.tri(dummycoph)
  rownames(dummycoph) <- colnames(dummycoph) <- allOrgs
  if (Verbose) pb <- txtProgressBar(max=length(pw), style=3)
  for ( i in seq_along(pw) ){
    if ( !skip || i %in% locs ){
      dummycoph[] <- NA
      cop <- NA
      # This is occasionally throwing errors that don't affect output for some reason
      cop <- as.matrix(Cophenetic(pw[[i]]))
      copOrgNames <- rownames(cop)
      if (useColoc){
        copOrgNames <- vapply(copOrgNames, gsub, pattern='(.+)_.+_[0-9]+', 
                              replacement='\\1', FUN.VALUE=character(1))
        rownames(cop) <- colnames(cop) <- copOrgNames
      }
      dummycoph[copOrgNames, copOrgNames] <- cop
      outmat[,i] <- dummycoph[ut]
    }
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\n')
  colnames(outmat) <- cols
  if (!is.null(toEval)){
    outmat <- outmat[,locs]
    #ltr <- vapply(seq_len(nrow(outmat)), function(x) all(is.na(outmat[x,])),
    #              FUN.VALUE=logical(1))
    #outmat <- outmat[!ltr,]
  }
  return(outmat)
}

RandCophProfiles.ProtWeaver <- function(pw, toEval=NULL, Verbose=TRUE, 
                                      speciesList=NULL, outdim=-1, 
                                      speciesCorrect=FALSE, mySpeciesTree=NULL, ...){
  ## TODO: Some way to handle paralogs
  cols <- names(pw)
  ao <- attr(pw, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(pw, 'useColoc')
  useMT <- attr(pw, 'useMT')
  
  stopifnot('ProtWeaver object must be initialized with dendrograms to run MirrorTree methods'=
              useMT)
  
  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  l <- length(allOrgs)
  num_entries <- as.integer((l * (l-1)) / 2)
  outdim <- ifelse(outdim < 1, l, outdim)
  outdim <- as.integer(outdim)
  outmat <- matrix(0, nrow=outdim, ncol=length(pw))
  dummycoph <- matrix(NA, nrow=l, ncol=l)
  ut <- upper.tri(dummycoph, diag=FALSE)
  
  if (speciesCorrect && !is.null(mySpeciesTree)){
    specvec <- c(Cophenetic(mySpeciesTree))
    spv2 <- specvec
    spv2[spv2==0] <- 1
    #nonzeros <- which(specvec != 0)
    #specvec <- .Call("randomProjection", specvec, nonzeros, length(nonzeros), outdim)  
  }
  
  rownames(dummycoph) <- colnames(dummycoph) <- allOrgs
  if (Verbose) pb <- txtProgressBar(max=length(pw), style=3)
  for ( i in seq_along(pw) ){
    if ( !skip || i %in% locs ){
      dummycoph[] <- 0
      cop <- 0
      # This is occasionally throwing errors that don't affect output for some reason
      cop <- as.matrix(Cophenetic(pw[[i]]))
      copOrgNames <- rownames(cop)
      if (useColoc){
        copOrgNames <- vapply(copOrgNames, gsub, pattern='(.+)_.+_[0-9]+', 
                              replacement='\\1', FUN.VALUE=character(1))
        rownames(cop) <- colnames(cop) <- copOrgNames
      }
      dummycoph[copOrgNames, copOrgNames] <- cop
      copvec <- dummycoph[ut]
      pos <- which(copvec != 0)
      if (speciesCorrect){
       copvec[pos] <- (copvec[pos] - specvec[pos]) / spv2[pos]
      }
      copvec <- .Call("randomProjection", copvec, pos, length(pos), outdim)  
      copvec[copvec==0] <- NA
      outmat[,i] <- copvec
    }
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\n')
  colnames(outmat) <- cols
  if (!is.null(toEval)){
    outmat <- outmat[,locs]
    #ltr <- vapply(seq_len(nrow(outmat)), function(x) all(is.na(outmat[x,])),
    #              FUN.VALUE=logical(1))
    #outmat <- outmat[!ltr,]
  }
  return(outmat)
}

ProcessSubset <- function(pw, Subset=NULL){
  pl <- length(pw)
  evalmap <- NULL
  uvals <- seq_len(pl)
  if (!is.null(Subset)){
    if (is(Subset, 'data.frame')){
      Subset <- as.matrix(Subset)
    }
    n <- names(pw)
    stopifnot("'Subset' must be either character, numeric, or matrix"=
                (is(Subset, 'character') || is(Subset, 'numeric') || is(Subset, 'matrix')))
    if (is(Subset, 'matrix')){
      if( ncol(Subset) != 2)
        stop('If Subset is a matrix, it must have 2 columns')
      if( is(Subset[1], 'character') ){
        Subset <- matrix(vapply(c(Subset), function(x) {
          val <- which(x==n)
          val <- ifelse(length(val) == 0, -1, val[1])
          }, 0), ncol=2)
        excise <- (Subset[,1] < 0) | (Subset[,2] < 0)
        if (sum(excise) > 0) 
          Subset <- Subset[!excise,]
      }
      for ( i in seq_len(nrow(Subset))){
        pos <- Subset[i,]
        i1 <- as.character(min(pos))
        i2 <- max(pos)
        evalmap[[i1]] <- c(evalmap[[i1]], i2)
      }
      uvals <- unique(c(Subset))
    } else {
      if (is(Subset, 'character'))
        Subset <- which(vapply(uvals, function(x) x == n, FUN.VALUE=logical(1)))
      if (length(Subset)==1){
        entry <- Subset[1]
        if (entry > 1){
          evalmap <- lapply(seq_len(entry-1), \(x) entry)
        } else {
          evalmap <- list()
        }
        
        if (entry < length(n)){
          evalmap[[entry]] <- seq(entry+1, length(n))
        }
        names(evalmap) <- as.character(seq(1, entry))
        uvals <- seq_along(n)
      } else {
        uvals <- unique(Subset)
        evalmap <- lapply(uvals, function(x) uvals)
        names(evalmap) <- as.character(uvals)
      }
    }
  }
  
  return(list(evalmap=evalmap, uvals=uvals))
}

flatdendrapply <- function(dend, NODEFUN=NULL, LEAFFUN=NODEFUN, 
                           INCLUDEROOT=TRUE, ...){
  stopifnot("flatdendrapply only works on objects of class 'dendrogram'"=
              is(dend, 'dendrogram'))
  if (!is(NODEFUN, 'function') && !is(LEAFFUN, 'function'))
    stop("At least one of NODEFUN and LEAFFUN must be a function!")
  
  val <- lapply(dend, 
                \(x){
                  if (is.null(attr(x, 'leaf'))){
                    if (!is(NODEFUN, 'function'))
                      v <- list()
                    else
                      v <- list(NODEFUN(x, ...))
                    for ( child in x ) v <- c(v, Recall(child))
                    return(v)
                  } 
                  else if (!is(LEAFFUN, 'function'))
                    return(list())
                  else 
                    return(list(LEAFFUN(x, ...)))
                }
  )
  retval <- unlist(val, recursive=FALSE)
  if (!INCLUDEROOT)
    retval[[1]] <- NULL
  
  lens <- vapply(retval, length, FUN.VALUE=0L)
  atom <- vapply(retval, is.atomic, FUN.VALUE=TRUE)
  
  if(all(lens == 1L) && all(atom))
    retval <- unlist(retval)
  
  return(retval)
}

AdjMatToDf <- function(preds, Verbose=TRUE){
  stopifnot(length(preds) > 0)
  n <- names(preds)
  prednames <- names(preds[[1]])
  lp <- length(prednames)
  v1 <- v2 <- character(lp*(lp+1)/2)
  ctr <- 1
  # can't use expand.grid because of duplicates but this is still fast
  for ( i in seq_len(lp) ){
    for (j in i:lp){
      v1[ctr] <- prednames[i]
      v2[ctr] <- prednames[j]
      ctr <- ctr + 1
    }
  }
  AdjDf <- data.frame(Gene1=v1, Gene2=v2)
  if (Verbose) pb <- txtProgressBar(max=length(n), style=3)
  for ( i in seq_len(length(n))){
    AdjDf[,n[i]] <- unclass(preds[[i]])
    if (Verbose) setTxtProgressBar(pb, i)
  }
  cat('\n')
  
  nc <- ncol(AdjDf)
  rtk <- vapply(seq_len(nrow(AdjDf)), function(x) sum(is.na(AdjDf[x,3:nc])) < (nc/2 - 1),
                FUN.VALUE=TRUE)
  AdjDf <- AdjDf[rtk,]
  rownames(AdjDf) <- NULL
  return(AdjDf)
}

PAStats <- function(predDf, paps){
  l <- nrow(paps)
  ocv1 <- vapply(predDf[,1], function(x) sum(paps[,x]) / l, 0)
  ocv2 <- vapply(predDf[,2], function(x) sum(paps[,x]) / l, 0)
  
  d <- abs(ocv1 - ocv2)
  av <- (ocv1 + ocv2) / 2
  return(list(avg=av, diff=d))
}

CheckBifurcating <- function(dend){
  # Checks if a dendrogram is bifurcating
  helperfunc <- function(node){
    if (length(node) == 1) return(TRUE)
    if (length(node) > 2) return(FALSE)
    return(helperfunc(node[[1]]) & helperfunc(node[[2]]))
  }
  return(helperfunc(dend))
}

DCA_minimize_fxn <- function(params, R, spins, i){
  sigi <- spins[,i]
  sigj <- spins[,-i]
  Hi <- params[i]
  Jij <- params[-i]
  
  sigprods <- sigi * sigj
  firstterm <- colSums(Jij * t(sigprods))
  secondterm <- Hi * sigi
  Si <- mean(exp( -1 * (firstterm + secondterm)))
  
  regularizer <- R * sum(abs(Jij))
  retval <- log(Si + regularizer)
  if (is.infinite(retval)){
    retval <- -1 * .Machine$double.xmax
  }
  return(retval)
}

DCA_gradient_minimize_fxn <- function(params, R, spins, i){
  grad <- numeric(length(params))
  
  sigi <- spins[,i]
  sigj <- spins[,-i]
  Hi <- params[i]
  Jij <- params[-i]
  
  #Calculate gradient
  sigprods <- sigi * sigj
  firstterm <- colSums(Jij * t(sigprods))
  secondterm <- Hi * sigi
  interior <- exp(-1 * firstterm + secondterm)
  entries <- interior * sigi * sigj
  Si <- mean(interior)
  if (is.null(dim(entries))){
    Jpartials <- (mean(entries) / Si) - R
  }
  else
    Jpartials <- (colMeans(entries) / Si) - R
  Hpartial <- mean(interior * sigi) / Si
  grad[-i] <- Jpartials
  grad[i] <- Hpartial
  return(grad)
}

DCA_logrise_run <- function(spins, links, regterm, printProgress=FALSE, NumCores=1){
  if (NumCores != 1){
    availableCores <- max(detectCores() - 1, 1) #leave an extra core just in case
    if (is.na(availableCores)) 
      NumCores <- 1
    else
      NumCores <- ifelse(NumCores < 0, availableCores, min(availableCores, NumCores))
  }
  nnodes <- ncol(spins)
  if (printProgress){
    cat('  Finding topology...\n')
    charsperline <- 73
    charsvec <- diff(floor(seq(0, charsperline, by=charsperline/nnodes))) == 1
  }
  
  if (printProgress) cat('  |')
  links <- simplify2array(mclapply(seq_len(nnodes), 
                                   function(i){
                                     probs <- links[,i]
                                     val <- optim(probs, DCA_minimize_fxn, 
                                                  gr=DCA_gradient_minimize_fxn,
                                                  method='BFGS',
                                                  control=list(reltol=1e-6),
                                                  R=regterm, spins=spins, i=i)$par
                                     if (printProgress && charsvec[i]) 
                                       system2('printf', '=')
                                     return(val)
                                   }, mc.cores=NumCores, mc.preschedule = TRUE
  )
  )
  if (printProgress) cat('| ')
  for ( i in seq_len(nrow(links)-1) ){
    for ( j in (i+1):ncol(links)){
      links[i,j] <- links[j,i] <- mean(links[i,j], links[j,i])
    }
  }
  
  if (printProgress) cat('\n  Done.\n  Eliminating edges close to zero (may take a moment)...\n')
  # Shrink close to zero to zero
  vals <- links[upper.tri(links)]
  h <- hist(vals, breaks=40, plot=FALSE)
  bin0 <- which(h$breaks==0)
  countsNorm <- h$counts == 0
  lb <- length(countsNorm)
  if (length(bin0) != 0 && lb > bin0){
    up <- countsNorm[bin0:lb]
    down <- countsNorm[bin0:1]
    double_up <- which(up[-length(up)] & up[-1])
    double_down <- which(down[-length(down)] & down[-1])
    safeguards <- quantile(vals, c(0.15,0.85))
    
    if ( length(double_up) == 0 ){
      tmp <- which(up)
      double_up <- ifelse(length(tmp) != 0, which(up)[1], NA)
    } else {
      double_up <- double_up[1]
    }
    uthresh <- ifelse(is.na(double_up), safeguards[2], min(h$breaks[bin0+double_up+1], safeguards[2]))
    
    if ( length(double_down) == 0 ){
      tmp <- which(down)
      double_down <- ifelse(length(tmp) != 0, which(down)[1], NA)
    } else {
      double_down <- double_down[1]
    }
    lthresh <- ifelse(is.na(double_down), safeguards[1], max(h$breaks[bin0-double_down-1], safeguards[1]))
    
    d <- diag(links)
    links[(links < uthresh & links > 0) | (links > lthresh & links < 0)] <- 0
    diag(links) <- d
  } 
  
  max_degree <- max(rowSums(links != 0))
  max_val <- max(links) * 2 #conservative estimate on link strength
  bounding <- exp(max_degree * max_val)
  if (printProgress) cat('  Refining edge weights...\n')
  
  if (printProgress) cat('  |')
  links <- simplify2array(
    mclapply(seq_len(nnodes), 
             function(i) {
               probs <- links[,i]
               mask <- seq_len(nnodes) == i
               nonzeros <- (probs != 0) | mask #have to include the self val
               adjustment <- ifelse(i==1, 0, sum(!nonzeros[seq_len(i-1)]))
               
               if ( sum(nonzeros) > 1 ){
                 pspins <- spins[,nonzeros]
                 pprobs <- probs[nonzeros]
                 val <- optim(pprobs, DCA_minimize_fxn, 
                              gr=DCA_gradient_minimize_fxn,
                              method='L-BFGS-B', 
                              lower=-bounding, upper=bounding,
                              R=0, spins=pspins, i=(i-adjustment))$par
                 probs[nonzeros] <- val
               }
               if (printProgress && charsvec[i]) system2('printf', '=')
               return(probs)
             }, 
             mc.cores=NumCores, mc.preschedule=TRUE))  
  if (printProgress) cat('| ')
  
  for ( i in seq_len(nrow(links)-1) ){
    for ( j in (i+1):ncol(links)){
      links[i,j] <- links[j,i] <- mean(links[i,j], links[j,i])
    }
  }
  
  if (max(abs(links)) != 0)
    links <- links / max(abs(links))
  if (printProgress) cat('\n  Done.\n')
  return(links)
}


DCA_logRISE <- function(PAProfiles, niter=1, reg_const=1, 
                        NumCores=1, zero_cutoff=0, Verbose=TRUE, ...){
  
  mult <- 1/niter
  intPA <- PAProfiles + 0
  intPA[intPA==0] <- -1
  nc <- ncol(intPA)
  pp <- FALSE
  
  if (Verbose & niter==1) pp <- TRUE  
  else if (Verbose){
    cat('Running DCA with', niter, 'iterations:\n')
    pb <- txtProgressBar(max=niter, style=3)
  } 
  
  truelinks <- countsmat <- matrix(0, nrow=nc, ncol=nc)
  for ( i in seq_len(niter) ){
    #initlinks <- matrix(0, nrow=nc, ncol=nc)
    initlinks <- matrix(rnorm(nc**2), nrow=nc)
    iterlink <- DCA_logrise_run(intPA, initlinks, reg_const, 
                                printProgress=pp, NumCores=NumCores)
    countsmat <- countsmat + (iterlink != 0)
    truelinks <- truelinks + mult * iterlink
    if(Verbose & !pp) setTxtProgressBar(pb, i)
  }
  if(Verbose & !pp) cat('\n')
  
  if (zero_cutoff > 0) {
    cutoff <- ifelse(zero_cutoff > 1, zero_cutoff, zero_cutoff * niter)
    countsmat[countsmat < cutoff] <- 0
    countsmat[countsmat > 0] <- 1
    truelinks <- truelinks * countsmat
  } 
  
  return(truelinks)
}

MICalc_C <- function(v1, v2, uv, pseudocount=1L){
  stopifnot("'pseudocount' must be an integer"=is(pseudocount, 'integer'))
  pseudocount <- min(0, pseudocount)
  # Psuedocount=0.01 taken as default by Gerardos et al. (2022) in PLoS
  a <- .Call('calcMIcVec', v1, v2, uv, pseudocount)
  return(a)
}

CorrComp_C <- function(fm, fsp, ssp, nv, nr){
  on.exit(.C("cleanupFxn"))
  stopifnot(nv == as.integer(nv))
  stopifnot(nr == as.integer(nr))
  stopifnot(all(fsp == as.integer(fsp)))
  stopifnot(all(ssp == as.integer(ssp)))
  
  a <- .Call('trimCovar', fm, as.integer(fsp), as.integer(ssp), as.integer(nv), as.integer(nr))
  return(a)
}

ResidueMIDend <- function(dend1, dend2, cutoff=0.9, comppct=0.25, useColoc, ...){
  if (useColoc){
    l2 <- gsub('([0-9]*)_.*', '\\1', labels(dend1))
    l1 <- gsub('([0-9]*)_.*', '\\1', labels(dend2))
  } else {
    l1 <- labels(dend1)
    l2 <- labels(dend2)
  }
  completeSet <- intersect(l1, l2)
  if (length(completeSet) == 0){
    return(0)
  }
  
  edges1 <- flatdendrapply(dend1, 
                           \(x) list(vals=as.character(unlist(x)), 
                                     state=attr(x, 'state')), 
                           NULL)
  edges2 <- flatdendrapply(dend2, 
                           \(x) list(vals=as.character(unlist(x)), 
                                     state=attr(x, 'state')), 
                           NULL)
  
  
  jsscore <- matrix(Inf, nrow=length(edges1), ncol=length(edges2))
  for ( i in seq_along(edges1) ){
    v1 <- intersect(edges1[[i]]$vals, completeSet)
    for ( j in seq_along(edges2) ){
      v2 <- intersect(edges2[[j]]$vals, completeSet)
      s <- 1 - length(intersect(v1, v2)) / length(union(v1, v2))
      jsscore[i,j] <- ifelse(is.nan(s), 1, s)
    }
  }
  
  nr <- nrow(jsscore)
  nc <- ncol(jsscore)
  if (nr < nc){
    tm <- edges1
    edges1 <- edges2
    edges2 <- tm
    jsscore <- t(jsscore)
    tm <- nc
    nc <- nr
    nr <- tm
  }
  #now guaranteed to have the larger dimension be nrow
  
  rownames(jsscore) <- as.character(seq_len(nr))
  colnames(jsscore) <- as.character(seq_len(nc))
  
  # I'm just using a greedy matching here, couldn't figure out Hungarian
  # and this also scales much better
  pairings <- rep(NA, nc)
  allvals <- rownames(jsscore)
  for ( i in seq_len(nc) ){
    ordered <- allvals[order(jsscore[,i])]
    pos <- which.min(ordered %in% pairings)
    if (jsscore[ordered[pos], i] < cutoff)
      pairings[i] <- ordered[pos]
  }
  # need to catch this here, in case we don't find enough elements just pick
  # random ones
  checksum <- sum(is.na(pairings))
  if (checksum > 0){
    possible <- allvals[!(allvals %in% pairings)]
    pairings[is.na(pairings)] <- sample(possible, checksum)
  }
  names(pairings) <- colnames(jsscore) 
  
  seqset1 <- seqset2 <- NULL
  n <- names(pairings)
  for ( i in seq_along(pairings) ){
    a1 <- as.integer(pairings[i])
    a2 <- as.integer(n[i])
    if (i == 1){
      seqset1 <- BStringSet(edges1[[a1]]$state)
      seqset2 <- BStringSet(edges2[[a2]]$state)
    } else {
      seqset1 <- append(seqset1, edges1[[a1]]$state)
      seqset2 <- append(seqset2, edges2[[a2]]$state)
    }
  }
  
  names(seqset1) <- names(seqset2) <- seq_len(length(pairings))
  res <- MISeqLevel(seqset1, seqset2, compressionpct=comppct)
  return(res)
}

MISeqLevel <- function(seqSet1, seqSet2, compressionpct=0.25){
  stopifnot('seqSets must be XStringSets'=is(seqSet1, 'XStringSet') && is(seqSet2, 'XStringSet'))
  stopifnot('seqSetq sequences have differing lengths. Ensure you are using an aligned sequence set.'=
              all(width(seqSet1) == width(seqSet1[1])))
  stopifnot('seqSet2 sequences have differing lengths. Ensure you are using an aligned sequence set.'=
              all(width(seqSet2) == width(seqSet2[1])))
  stopifnot('compressionpct must be between 0 and 1'=
              compressionpct < 1 && compressionpct > 0)
  stopifnot('seqSets must be named'=!is.null(names(seqSet1)) && !is.null(names(seqSet2)))
  stopifnot('seqSets must be named'=all(!is.na(c(names(seqSet1), names(seqSet2)))))
  #stopifnot('Both inputs must be DNAStringSets'=
  #            is(seqSet1, 'DNAStringSet') && is(seqSet2, 'DNAStringSet'))
  start2 <- width(seqSet1)[1] + 1
  cali <- ConcatSeqs(seqSet1, seqSet2)
  if (length(cali) == 0){
    #warning('No sequences shared. Check seqSet names!')
    return(0)
  }
  
  v <- CorrCompressSeqs(cali, start2, mvalpct=compressionpct)
  if (!is.null(v$warn)){
    #warning('Sequences identical.')
    return(1)
  }
  compali <- v$xstrset
  pos <- v$pos
  newstart2 <- which.max(pos >= start2)
  miscore <- CalcMIReduced(compali, newstart2)
  
  # APC correction
  nr <- nrow(miscore)
  nc <- ncol(miscore)
  if (nr == 0 || nc == 0){
    return(0)
  }
  APC_corr <- matrix(colMeans(miscore), nr, nc, byrow = TRUE) * 
    matrix(rowMeans(miscore), nr, nc, byrow = FALSE) / mean(miscore)
  miscore <- miscore - APC_corr
  
  # scoring
  miscore <- apply(abs(miscore), 2, max)
  retval <- mean(miscore)
  if (is.nan(retval)) 
    retval <- 0
  return(retval)
}

ConcatSeqs <- function(seqSet1, seqSet2){
  unames <- intersect(names(seqSet1), names(seqSet2))
  concatAli <- xscat(seqSet1[unames], seqSet2[unames])
  names(concatAli) <- unames
  return(concatAli)
}

CorrCompressSeqs <- function(myStringSet, start2, pseudocount=2, mvalpct=0.5, 
                             gapLetters=c('-', '.'),
                             uncertainty_cutoff=0.158, MAF_cutoff=0.15){
  freqMat <- consensusMatrix(myStringSet, as.prob=FALSE)
  freqMat <- freqMat[rowSums(freqMat) != 0,]
  freqMat <- freqMat + pseudocount
  
  freqMat <- t(t(freqMat) / colSums(freqMat))
  
  
  to_keep <- rep(FALSE, ncol(freqMat))
  nongaploc <- !(rownames(freqMat) %in% gapLetters)
  for ( i in seq_along(to_keep) ){
    pos <- freqMat[,i]
    pos_no_gap <- pos[nongaploc]
    
    missing_prob <- sum(pos[gapLetters], na.rm=TRUE)
    vals <- sort(pos_no_gap, decreasing=TRUE)
    
    MAF <- ifelse(vals[1] == 0, 0, vals[2] / (vals[1] + vals[2]))
    
    to_keep[i] <- missing_prob < uncertainty_cutoff && MAF > MAF_cutoff
  }
  # Need a guard case here
  if (!any(to_keep)){
    to_keep <- !to_keep
  }
  trimmedFreqMat <- freqMat[,to_keep]
  colnames(trimmedFreqMat) <- which(to_keep)
  fm <- unique(trimmedFreqMat, MARGIN=2)
  
  nc <- ncol(fm)
  num_vals <- ceiling(mvalpct * nc) 
  
  isInSecondSeq <- as.integer(colnames(fm)) >= start2
  s2 <- which.max(isInSecondSeq)
  firstSeqPos <- seq_len(s2-1)
  if (length(firstSeqPos) == 0){
    return(list(warn=TRUE, 1))
  }
  secondSeqPos <- seq_len(length(isInSecondSeq) - s2 + 1L) + s2 - 1L
  # keeping the entire covariance matrix in memory is really hard
  # this is more code but significantly more memory efficient
  # the cost is slightly more runtime
  corrs <- CorrComp_C(fm, firstSeqPos, secondSeqPos, num_vals, nrow(fm))
  corrs <- unique(corrs)
  if ( length(corrs) < num_vals ){
    fsp <- firstSeqPos[!(firstSeqPos %in% corrs)]
    ssp <- secondSeqPos[!(secondSeqPos %in% corrs)]
    d <- ceiling((num_vals - length(corrs)) / 2)
    l1 <- min(d, length(fsp))
    l2 <- min(d, length(ssp))
    corrs <- c(corrs, sample(fsp, l1), sample(ssp, l2))
  }
  upositions <- sort(as.integer(colnames(fm)[unique(corrs)]))
  subsetXStr <- extractAt(myStringSet, IRanges(upositions, width=1))
  return(list(xstrset=unstrsplit(subsetXStr), pos=upositions))
}

CalcMIReduced <- function(trimmedXStringSet, start2,
                          secondgroupstart=-1){
  matxss <- as.matrix(trimmedXStringSet)
  group1 <- seq_len(start2-1)
  group2 <- seq_len(width(trimmedXStringSet[1]) - start2) + start2
  
  u <- unique(c(matxss))
  converter <- seq_len(length(u)) - 1L
  names(converter) <- u
  umat <- matrix(converter[matxss], ncol=ncol(matxss))
  scores <- matrix(NA, nrow=length(group1), ncol=length(group2))
  gapnum <- ifelse('-' %in% u, which(u=='-'), -1)
  numunique <- length(u) - (gapnum!=-1)
  ctr <- 0
  for ( i in seq_along(group1) ){
    p1 <- umat[,group1[i]]
    subsetloc <- p1 != gapnum
    for ( j in seq_along(group2) ){
      p2 <- umat[,group2[j]]
      
      fullsub <- subsetloc & (p2 != gapnum)
      p1p <- p1[fullsub]
      p2p <- p2[fullsub]
      MI <- MICalc_C(p1p, p2p, numunique)
      scores[i,j] <- MI
      ctr <- ctr + 1
    }
  }
  return(scores)
}

predictWithBuiltins <- function(preds){
  # Key: (val is binary + 1)
  # 000 => 1: Jaccard, Hamming, MI, ProfileDCA (base)
  # 001 => 2: base and MT
  # 010 => 3: base and Behdenna
  # 011 => 4: base and Behdenna, MT
  # 100 => 5: base and Coloc
  # 101 => 6: base and Coloc, MT
  # 110 => 7: base and Behdenna, Coloc
  # 111 => 8: base and Behdenna, Coloc, MT
  modelsToUse <- rep(1, nrow(preds))
  relevant_cnames <- c('MirrorTree', 'Behdenna', 'Coloc')
  pred_cnames <- colnames(preds)
  for (i in seq_along(relevant_cnames)){
    rcn <- relevant_cnames[i]
    if (rcn %in% pred_cnames){
      idxs <- !is.na(preds[,rcn])
      modelsToUse[idxs] <- modelsToUse[idxs] + (2**(i-1)) 
    }
  }
  builtins <- get(data('BuiltInEnsembles', envir=environment()))
  if (all(modelsToUse == modelsToUse[1])){
    return(predict(builtins[[modelsToUse[1]]], preds, type='response'))
  } else {
    builtInPredictions <- rep(NA, nrow(preds))
    for (i in seq_along(builtInPredictions)){
      model <- builtins[[modelsToUse[i]]]
      builtInPredictions[i] <- predict(model, preds[i,], type='response')
    }
    return(builtInPredictions)
  }
}

findSpeciesTree <- function(pw, Verbose=TRUE, NameFun=NULL){
  stopifnot("ProtWeaver object must contain dendrograms"=attr(pw, "useMT"))
  if (attr(pw, "useColoc") && is.null(NameFun)){
    NameFun <- function(x) gsub('([^_])_.*', '\\1', x)
  }
  
  SpecTree <- SuperTree(unclass(pw), NAMEFUN=NameFun, Verbose=Verbose)
  
  return(SpecTree)
}
########