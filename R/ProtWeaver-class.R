##### -- ProtWeaver Class to find clusters of functionally linked genes -------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### DEFINITION ####
# Class expects as input one of the following:
#   1. A list of dendrograms
#      a) Dendrograms with labels [ORG]_[INDEX]_[POS] will use co-loc code
#      b) Dendrograms without will not use co-loc
#   2. A list of character vectors defining COGs (will not use MirrorTree)
#   3. A list of character vectors defining COGs with labels like 1(a)
####################


#### S3 Generic Definitions ####
PAProfiles <- function(pw, ...) UseMethod('PAProfiles')
CophProfiles <- function(pw, ...) UseMethod('CophProfiles')
Jaccard <- function(pw, ...) UseMethod('Jaccard')
Hamming <- function(pw, ...) UseMethod('Hamming')
MutualInformation <- function(pw, ...) UseMethod('MutualInformation')
MirrorTree <- function(pw, ...) UseMethod('MirrorTree')
ContextTree <- function(pw, ...) UseMethod('ContextTree')
ProfileDCA <- function(pw, ...) UseMethod('ProfileDCA')
Coloc <- function(pw, ...) UseMethod('Coloc')
Behdenna <- function(pw, ...) UseMethod('Behdenna')
ResidueMI <- function(pw, ...) UseMethod('ResidueMI')
Ensemble <- function(pw, ...) UseMethod('Ensemble')
########


#### Class Constructor ####
new_ProtWeaver <- function(validatedInput){
  structure(validatedInput$ipt,
            allOrgs=validatedInput$allgenomes,
            useMT=validatedInput$flags$usemirrortree,
            useColoc=validatedInput$flags$usecoloc,
            useResidue=validatedInput$flags$useresidue,
            class='ProtWeaver')
}

ProtWeaver <- function(ListOfData, NoWarn=FALSE){
  vRes <- validate_ProtWeaver(ListOfData, noWarn=NoWarn)
  new_ProtWeaver(vRes)
}

validate_ProtWeaver <- function(ipt, noWarn=FALSE){
  bitflags <- list(usecoloc=FALSE, usemirrortree=FALSE, useresidue=FALSE)
  stopifnot('ProtWeaver expects a list of dendrograms or character vectors as input.'=
              is(ipt, 'list'))
  
  stopifnot('Input has no groups!'=length(ipt)>0)
  checkdend <- vapply(ipt, is, class2='dendrogram', FUN.VALUE = TRUE)
  checkchar <- vapply(ipt, is, class2='character', FUN.VALUE = TRUE)
  stopifnot('Input list must be all vectors or all dendrograms'=
              !any(checkchar) || all(checkchar))
  stopifnot('Input list must be all vectors or all dendrograms'=
              !any(checkdend) || all(checkdend))
  
  stopifnot('Input list elements must be dendrograms or vectors of character'=
              all(checkdend) || all(checkchar))
  
  # Now we know that the input is either of type 'character' or 'dendrogram'
  if (all(checkchar)){
    if (!noWarn) message('Disabling Residue and MirrorTree-based algorithms. ',
                  'Input list must include inputs of type "dendrogram" ',
                  'for MT algorithms. Consult the documentation for more info.\n')
    bitflags[['usemirrortree']] <- FALSE
    allentries <- unique(unlist(ipt))
  } else {
    bitflags[['usemirrortree']] <- TRUE
    allentries <- unique(unlist(lapply(ipt, labels)))
  }
  
  if (bitflags[['usemirrortree']]){
    useResidueMI <- TRUE
    for ( tree in ipt){
      if (any(unlist(flatdendrapply(tree, \(x) is.null(attr(x, 'state')))))){
        useResidueMI <- FALSE
        if (!noWarn) message('Disabling Residue methods. Input dendrograms must',
                             ' include ancenstral state reconstruction for residue',
                             ' methods. Consult the documentation for more info.\n')
        break
      }
    }
    bitflags[['useresidue']] <- useResidueMI
  }
  
  checkforcoloc <- grepl('.+_.+_[0-9]+', allentries)
  if ( all(checkforcoloc) ){
    bitflags[['usecoloc']] <- TRUE
    allentries <- unique(gsub('(.+)_.+_[0-9]+', '\\1', allentries))
  }
  else{
    bitflags[['usecoloc']] <- FALSE
    if (!noWarn) message('Co-localization disabled. Labels must be in the format ',
              '[GENOME]_[INDEX]_[ORDER] to use co-localization ',
            '(where ORDER is a numeric). Consult the documentation for more info.\n')
  }
  
  if (is.null(names(ipt))){
    if (!noWarn) message('Adding character labels to input data.\n')
    names(ipt) <- as.character(seq_len(length(ipt)))
  }
  if (any(names(ipt) == '')){
    if (!noWarn) message('Adding labels to unnamed groups.\n')
    safe <- paste0('Unnamed_Grp_', as.character(seq_len(length(ipt))))
    n <- names(ipt)
    names(ipt)[n==''] <- safe[n=='']
  }
  
  if (all(grepl('^[0-9]+$', allentries)))
    allentries <- allentries[order(as.integer(allentries))]
  else 
    allentries <- sort(allentries)
  
  return(list(ipt=ipt, allgenomes=allentries, flags=bitflags))
}
########

#### User-Exposed S3 Methods ####
show.ProtWeaver <- function(x, ...){
  if (length(x) == 1){
    cat(paste('a ProtWeaver object with', length(x),
              'group and', length(attr(x,'allOrgs')), 'genomes.\n'))
  } else {
    cat(paste('a ProtWeaver object with', length(x),
              'groups and', length(attr(x,'allOrgs')), 'genomes.\n'))
  }
}

print.ProtWeaver <- function(x, ...){
  if (length(x) == 1){
    cat(paste('a ProtWeaver object with', length(x),
              'group and', length(attr(x,'allOrgs')), 'genomes.\n'))
  } else {
    cat(paste('a ProtWeaver object with', length(x),
              'groups and', length(attr(x,'allOrgs')), 'genomes.\n'))
  }
}

`[.ProtWeaver` <- function(x, i){
  y <- unclass(x)
  newv <- validate_ProtWeaver(y[i], noWarn=TRUE)
  new_ProtWeaver(newv)
}

predict.ProtWeaver <- function(object, Method='Ensemble', Subset=NULL, NumCores=1,
                               MySpeciesTree=NULL, PretrainedModel=NULL,
                               RawZScores=FALSE, NoPrediction=FALSE, 
                               ReturnRawData=FALSE, Verbose=TRUE, ...){
  pw <- object
  func <- getS3method(Method, 'ProtWeaver')
  if(Verbose && !ReturnRawData) starttime <- Sys.time()
  
  preds <- func(pw, Subset=Subset, Verbose=Verbose, 
                MySpeciesTree=MySpeciesTree, NumCores=NumCores,
                PretrainedModel=PretrainedModel, RawZScores=RawZScores, 
                NoPrediction=NoPrediction, ...)
  
  if (Verbose && !ReturnRawData){
    cat('Done.\n\nTime difference of', 
        round(difftime(Sys.time(), starttime, units = 'secs'), 2),
        'seconds.\n')
  } 
  if (is(preds, 'list') && !is.null(preds$noPostFormatting))
    return(invisible(preds$res))
  if (ReturnRawData)
    return(invisible(preds))
  
  pc <- ProcessSubset(pw, Subset)
  n <- names(pw)[pc$uvals]
  rownames(preds) <- colnames(preds) <- n
  rs <- structure(preds,
                  method=Method,
                  class='ProtWeb')
  
  return(invisible(rs))
}

########


#### Helper Functions ####

ProcessSubset <- function(pw, Subset=NULL){
  pl <- length(pw)
  evalmap <- NULL
  uvals <- seq_len(pl)
  if (!is.null(Subset)){
    n <- names(pw)
    stopifnot("'Subset' must be either character, numeric, or matrix"=
                (is(Subset, 'character') || is(Subset, 'numeric') || is(Subset, 'matrix')))
    if (is(Subset, 'matrix')){
      if( ncol(Subset) != 2)
        stop('If Subset is a matrix, it must have 2 columns')
      if( is(Subset[1], 'character') ){
        Subset <- matrix(vapply(c(Subset), function(x) which(x==n), 0), ncol=2)
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

flatdendrapply <- function(dend, NODEFUN, LEAFFUN=NODEFUN, INCLUDEROOT=TRUE, ...){
  val <- lapply(dend, 
                \(x){
                  if (is.null(attr(x, 'leaf'))){
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
  
  return(retval)
}

AdjMatToDf <- function(preds){
  stopifnot(length(preds) > 0)
  expred <- preds[[1]]
  pair_locs <- upper.tri(expred)
  pairnames <- which(pair_locs, arr.ind=TRUE)
  pairentry1 <- rownames(expred)[pairnames[,'row']]
  pairentry2 <- colnames(expred)[pairnames[,'col']]
  AdjDf <- data.frame(Gene1=pairentry1, Gene2=pairentry2)
  n <- names(preds)
  for ( i in seq_along(n))
    AdjDf[,n[i]] <- (preds[[i]])[pair_locs]
  
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
  
  return(log(Si + regularizer))
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
    charsperline <- 82
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
  
  if (printProgress) cat('Done.\n  Eliminating edges close to zero (may take a moment)...\n')
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
  if (printProgress) cat('Done.\n')
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
  return(mean(miscore))
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
########


#### Internal S3 Methods ####

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
  
  pairscores <- matrix(NA, nrow=pl, ncol=pl)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(pl*(pl-1) / 2), style=3)
  for ( i in seq_len(pl-1) ){
    uval1 <- uvals[i]
    v1 <- CPs[,i]
    for ( j in (i+1):pl ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        v2 <- CPs[,j]
        val <- cor(v1, v2, use='na.or.complete', method='pearson')
        pairscores[i,j] <- pairscores[j,i] <- ifelse(is.na(val), 0, val)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  diag(pairscores) <- 1
  if ('partialcorrelation' %in% MTCorrection){
    flag <- TRUE
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
  }
  diag(pairscores) <- 1
  rownames(pairscores) <- colnames(pairscores) <- names(pw)[precalcSubset$uvals]
  
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

Jaccard.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE,
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap

  if ( is.null(precalcProfs) ){
    if (Verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  pairscores <- matrix(NA, nrow=l, ncol=l)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l) ){
    uval1 <- uvals[i]
    p1 <- pap[,i]
    for ( j in i:l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        p2 <- pap[,j]
        m10 <- sum(p1 & !p2)
        m01 <- sum(!p1 & p2)
        m00 <- sum(!p1 & !p2)
        dJ <- (m10 + m01) / (m10 + m01 + m00)
        dJ <- ifelse(is.nan(dJ), 0, dJ) 
        pairscores[i,j] <- pairscores[j,i] <- dJ
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  diag(pairscores) <- 0
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  pairscores <- 1 - pairscores # because distance
  return(pairscores)
}

Hamming.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  pairscores <- matrix(NA, nrow=l, ncol=l)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l) ){
    uval1 <- uvals[i]
    p1 <- pap[,i]
    for ( j in i:l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        p2 <- pap[,j]
        pairscores[i,j] <- pairscores[j,i] <- sum(xor(p1,p2)) / ncol(pap)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  diag(pairscores) <- 0
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  mp <- max(pairscores, na.rm=TRUE)
  pairscores <- (mp - pairscores) / mp #because distance
  return(pairscores)
}

MutualInformation.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                         precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, Verbose=Verbose)
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  n <- names(pw)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  pairscores <- matrix(NA, nrow=l, ncol=l)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l) ){
    uval1 <- uvals[i]
    v1 <- pap[,i]
    for ( j in i:l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        v2 <- pap[,j]
        score <- 0
        v1l <- length(v1)
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
        pairscores[i,j] <- pairscores[j,i] <- score
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  
  # Correction
  apccorr <- mean(pairscores[upper.tri(pairscores)], na.rm=TRUE)
  pairscores <- pairscores - apccorr
  pairscores <- abs(pairscores)
  # Normalize
  denom <- max(pairscores, na.rm=TRUE)
  pairscores <- pairscores / ifelse(denom==0, 1, denom)
  diag(pairscores) <- 1
  #pairscores <- pairscores #because distance
  return(pairscores)
}

ProfileDCA.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, NumCores=1,
                           precalcProfs=NULL, precalcSubset=NULL, useAbs=TRUE, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, Verbose=Verbose)
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
  
  pairscores <- DCA_logRISE(pap, Verbose=Verbose, NumCores=NumCores, ...)
  rownames(pairscores) <- colnames(pairscores) <- n
  if (useAbs) pairscores <- abs(pairscores)
  if (max(pairscores) != 0)
    pairscores <- pairscores / max(abs(pairscores))
  
  
  return(pairscores)
}

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

Behdenna.ProtWeaver <- function(pw, Subset=NULL, Verbose=TRUE, 
                                MySpeciesTree=NULL, useSubtree=FALSE, 
                                useACCTRAN=TRUE, rawZScores=FALSE, 
                                precalcProfs=NULL, precalcSubset=NULL, ...){
  stopifnot('No species tree provided.'=(!is.null(MySpeciesTree)))
  stopifnot("Method 'Behdenna' requires a bifurcating tree"=CheckBifurcating(MySpeciesTree))
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  n <- names(pw)[uvals]
  if ( is.null(precalcProfs) ){
    if (Verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, Verbose=Verbose, speciesList=labels(MySpeciesTree))
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
  
  vals <- calc_SId_mat(fd, IdOnly=!useSubtree)
  if (useSubtree)
    M <- vals$S
  else
    M <- vals$Id
  Cmat <- vals$Cmat
  bl <- vals$blengths
  
  glmat <- abs(glmat)
  pairscores <- matrix(NA, nrow=l, ncol=l)
  
  ctr <- 0
  if (Verbose) cat('\n  Calculating pairscores:\n')
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l)){
    uval1 <- uvals[i]
    gl1 <- glmat[,i]
    n1 <- sum(gl1)
    for ( j in i:l){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        gl2 <- glmat[,j]
        n2 <- sum(gl2)
        score <- 0
        if ( n1*n2 != 0 ){
          score <- t(gl1) %*% M %*% gl2
          exp_mean <- n2 * (t(gl1) %*% M %*% bl)
          exp_var <- n2*t(gl1) %*% M %*% Cmat %*% t(M) %*% gl1
          score <- (score - exp_mean) / sqrt(abs(exp_var))
        }
        pairscores[i,j] <- pairscores[j,i] <- score
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')
  diag(pairscores) <- 0 #because z-scores
  
  if (!rawZScores){
    pairscores <- abs(pairscores)
    pairscores <- pairscores / ifelse(max(pairscores,na.rm=TRUE) != 0, 
                                      max(pairscores, na.rm=TRUE), 1)
    diag(pairscores) <- 1
  }
  rownames(pairscores) <- colnames(pairscores) <- n
  return(pairscores)
}

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

Ensemble.ProtWeaver <- function(pw,
                                Subset=NULL, Verbose=TRUE, MySpeciesTree=NULL,
                                PretrainedModel=NULL,
                                NoPrediction=FALSE, ...){
  
  flags <- rep(FALSE, 3)
  
  subs <- ProcessSubset(pw, Subset)
  n <- names(pw)
  uvals <- subs$uvals
  unames <- vapply(uvals, function(x) n[x], FUN.VALUE=character(1))
  splist <- NULL
  if (!is.null(MySpeciesTree)){
    splist <- labels(MySpeciesTree)
  }
    
  if (Verbose) cat('Calculating P/A profiles:\n')
  PAs <- PAProfiles(pw, uvals, Verbose=Verbose, speciesList=splist)
  CPs <- NULL
  takesCP <- c('MirrorTree') # Just using MirrorTree for prediction
  
  submodels <- c('ProfileDCA', 'Jaccard', 'Hamming', 'MutualInformation')
  if (attr(pw, 'useMT')){
    if (Verbose) cat('Calculating Cophenetic profiles:\n')
    CPs <- CophProfiles(pw, uvals, Verbose=Verbose, speciesList=splist)
    submodels <- c(submodels, takesCP)
  }
  
  if (!is.null(MySpeciesTree) && CheckBifurcating(MySpeciesTree)){
    submodels <- c(submodels, 'Behdenna')
  }

  if (attr(pw, 'useColoc')){
    submodels <- c(submodels, 'Coloc')
  }
  
  if (attr(pw, 'useResidue')){
    # Ensemble with this still needs to be trained
    #submodels <- c(submodels, 'useResidue')
    submodels <- submodels
  }
  
  if(!is.null(PretrainedModel)) {
    UseBuiltIns <- FALSE
    predictionmodel <- PretrainedModel
  } else {
    UseBuiltIns <- TRUE
  }
  
  results <- list()
  for ( model in submodels ){
    if (Verbose) cat('Running ', model, ':\n', sep='')
    if (model %in% takesCP) profs <- CPs
    else profs <- PAs
    results[[model]] <- predict(pw, model, Verbose=Verbose, 
                              ReturnRawData=TRUE, precalcProfs=profs,
                              precalcSubset=subs, 
                              MySpeciesTree=MySpeciesTree, ...)
  }
  
  cat('Calculating additional P/A Statistics...\n')
  results <- AdjMatToDf(results)
  pas <- PAStats(results, PAs) 
  results[,'AvgOcc'] <- pas$avg
  results[,'OccDiff'] <- pas$diff

  if (NoPrediction) return(list(res=results, noPostFormatting=TRUE))
  
  if (Verbose) cat('Predicting with Ensemble method...\n')
  if (UseBuiltIns){
    predictions <- predictWithBuiltins(results)
  } else if (is(predictionmodel, 'glm')) {
    predictions <- predict(predictionmodel, results[,-c(1,2)], type='response')
  } else {
    predictions <- predict(predictionmodel, results[,-c(1,2)])
  }
  outmat <- matrix(NA, nrow=length(uvals), ncol=length(uvals))

  for (i in seq_along(predictions)){
    i1 <- which(results[i,1] == unames)
    i2 <- which(results[i,2] == unames)
    pred <- predictions[i]
    outmat[i1,i2] <- outmat[i2,i1] <- pred
  }
  rownames(outmat) <- colnames(outmat) <- uvals
  diag(outmat) <- 1
  return(outmat)
}

########