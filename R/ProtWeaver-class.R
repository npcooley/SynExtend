##### -- ProtWeaver Class to find clusters of functionally linked genes -------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

## Defining S3 class for ProtWeaver Object
builtinmodels <- '../data/EnsembleModels.RData'

#### DEFINITION ####
# Class expects as input one of the following:
#   1. A list of dendrograms
#      a) Dendrograms with labels [ORG]_[INDEX]_[POS] will use co-loc code
#      b) Dendrograms without will not use co-loc
#   2. A list of character vectors defining COGs (will not use MirrorTree)
####################

#### METHODS ####
# predict()
# plot()
# print()
# str() (?)
#################


#### S3 Generic Definitions ####
PAProfiles <- function(x, ...) UseMethod('PAProfiles')
CophProfiles <- function(x, ...) UseMethod('CophProfiles')
Jaccard <- function(x, ...) UseMethod('Jaccard')
Hamming <- function(x, ...) UseMethod('Hamming')
MutualInformation <- function(x, ...) UseMethod('MutualInformation')
MirrorTree <- function(x, ...) UseMethod('MirrorTree')
ContextTree <- function(x, ...) UseMethod('ContextTree')
ProfileDCA <- function(x, ...) UseMethod('ProfileDCA')
Coloc <- function(x, ...) UseMethod('Coloc')
Behdenna <- function(x, ...) UseMethod('Behdenna')
Ensemble <- function(x, ...) UseMethod('Ensemble')
########


#### Class Definitions ####
new_ProtWeaver <- function(validatedInput){
  structure(validatedInput$ipt,
            allOrgs=validatedInput$allgenomes,
            useMT=validatedInput$flags$usemirrortree,
            useColoc=validatedInput$flags$usecoloc,
            class='ProtWeaver')
}

ProtWeaver <- function(ipt){
  vRes <- validate_ProtWeaver(ipt)
  new_ProtWeaver(vRes)
}

validate_ProtWeaver <- function(ipt){
  bitflags <- list(usecoloc=F, usemirrortree=F)
  stopifnot('ProtWeaver expects a list of dendrograms or character vectors as input.'=
              is(ipt, 'list'))
  
  stopifnot('Input has no groups!'=length(ipt)>0)
  checkdend <- sapply(ipt, is, class2='dendrogram')
  checkchar <- sapply(ipt, is, class2='character')
  stopifnot('Input list must be all vectors or all dendrograms'=
              !any(checkchar) || all(checkchar))
  stopifnot('Input list must be all vectors or all dendrograms'=
              !any(checkdend) || all(checkdend))
  
  stopifnot('Input list elements must be dendrograms or vectors of character'=
              all(checkdend) || all(checkchar))
  
  # Now we know that the input is either of type 'character' or 'dendrogram'
  if (all(checkchar)){
    warning(paste('Disabling MirrorTree-based algorithms.',
                  'Input list must include inputs of type "dendrogram"',
                  'for MT algorithms. Consult the documentation for more info.'))
    bitflags[['usemirrortree']] <- F
    allentries <- unique(unlist(ipt))
  } else {
    bitflags[['usemirrortree']] <- T
    allentries <- unique(unlist(sapply(ipt, labels)))
  }
  
  checkforcoloc <- grepl('.+_.+_[0-9]+', allentries)
  if ( all(checkforcoloc) ){
    bitflags[['usecoloc']] <- T
    allentries <- unique(gsub('(.+)_.+_[0-9]+', '\\1', allentries))
  }
  else{
    bitflags[['usecoloc']] <- F
    warning(paste('Co-localization disabled. Labels must be in the format',
              '[GENOME]_[INDEX]_[ORDER] to use co-localization',
            '(where ORDER is a numeric). Consult the documentation for more info.'))
  }
  
  if (is.null(names(ipt))){
    warning('Adding character labels to input data.')
    names(ipt) <- as.character(1:length(ipt))
  }
  if (any(names(ipt) == '')){
    warning('Adding labels to unnamed groups.')
    safe <- paste0('Unnamed_Grp_', as.character(1:length(ipt)))
    n <- names(ipt)
    names(ipt)[n==''] <- safe[n=='']
  }
  
  if (all(grepl('^[0-9]+$', allentries)))
    allentries <- allentries[order(as.integer(allentries))]
  else 
    allentries <- sort(allentries)
  
  return(list(ipt=ipt, allgenomes=allentries, flags=bitflags))
}

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
  newv <- suppressWarnings(validate_ProtWeaver(y[i]))
  new_ProtWeaver(newv)
}

predict.ProtWeaver <- function(pw, method='Jaccard', subset=NULL, verbose=TRUE, 
                               ensembleprediction=FALSE, ...){
  
  func <- getS3method(method, 'ProtWeaver')
  preds <- func(pw, subset=subset, verbose=verbose, ...)
  
  if (is(preds, 'list') && !is.null(preds$noPostFormatting))
    return(preds$res)
  if (ensembleprediction)
    return(preds)
  
  rs <- structure(preds,
                  method=method,
                  class='ProtWeb')
  
  return(rs)
}

########


#### Helper Functions ####

ProcessSubset <- function(pw, subset=NULL){
  pl <- length(pw)
  evalmap <- NULL
  uvals <- 1:pl
  if (!is.null(subset)){
    n <- names(pw)
    stopifnot("'subset' must be either character, numeric, or matrix"=
                (is(subset, 'character') || is(subset, 'numeric') || is(subset, 'matrix')))
    if (is(subset, 'matrix')){
      if( ncol(subset) != 2)
        stop('If subset is a matrix, it must have 2 columns')
      if( is(subset[1], 'character') ){
        subset <- matrix(vapply(c(subset), function(x) which(x==n), 0), ncol=2)
      }
      for ( i in 1:nrow(subset)){
        pos <- subset[i,]
        i1 <- as.character(min(pos))
        i2 <- max(pos)
        evalmap[[i1]] <- c(evalmap[[i1]], i2)
      }
      uvals <- unique(c(subset))
    } else {
      if (is(subset, 'character'))
        subset <- which(sapply(uvals, function(x) x == n))
      uvals <- unique(subset)
      evalmap <- lapply(uvals, function(x) uvals)
      names(evalmap) <- as.character(uvals)
    }
  }
  
  return(list(evalmap=evalmap, uvals=uvals))
}

AdjMatToDf <- function(preds){
  stopifnot(length(preds) > 0)
  expred <- preds[[1]]
  pair_locs <- upper.tri(expred)
  pairnames <- which(pair_locs, arr.ind=T)
  pairentry1 <- rownames(expred)[pairnames[,'row']]
  pairentry2 <- colnames(expred)[pairnames[,'col']]
  AdjDf <- data.frame(Gene1=pairentry1, Gene2=pairentry2)
  n <- names(preds)
  for ( i in seq_along(n))
    AdjDf[,n[i]] <- (preds[[i]])[pair_locs]
  
  nc <- ncol(AdjDf)
  rtk <- vapply(1:nrow(AdjDf), function(x) sum(is.na(AdjDf[x,3:nc])) < (nc/2 - 1),
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

DCA_logrise_run <- function(spins, links, regterm, printProgress=F, numCores=1){
  if (numCores != 1){
    require(parallel)
    availableCores <- detectCores()
    numCores <- ifelse(numCores < 0, availableCores, min(availableCores, numCores))
  }
  nnodes <- ncol(spins)
  if (printProgress){
    cat('  Finding topology...\n')
    charsperline <- 82
    charsvec <- diff(floor(seq(0, charsperline, by=charsperline/nnodes))) == 1
  }
  
  if (printProgress) cat('  |')
  links <- simplify2array(mclapply(1:nnodes, 
                                   function(i){
                                     probs <- links[,i]
                                     val <- optim(probs, DCA_minimize_fxn, 
                                                  gr=DCA_gradient_minimize_fxn,
                                                  method='BFGS',
                                                  control=list(reltol=1e-6),
                                                  R=regterm, spins=spins, i=i)$par
                                     if (printProgress && charsvec[i]) system('printf =')
                                     return(val)
                                   }, mc.cores=numCores, mc.preschedule = FALSE
  )
  )
  if (printProgress) cat('| ')
  for ( i in 1:(nrow(links)-1) ){
    for ( j in (i+1):ncol(links)){
      links[i,j] <- links[j,i] <- mean(links[i,j], links[j,i])
    }
  }
  
  if (printProgress) cat('Done.\n  Eliminating edges close to zero (may take a moment)...\n')
  # Shrink close to zero to zero
  vals <- links[upper.tri(links)]
  h <- hist(vals, breaks=40, plot=F)
  bin0 <- which(h$breaks==0)
  countsNorm <- h$counts == 0
  lb <- length(countsNorm)
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
  
  max_degree <- max(rowSums(links != 0))
  max_val <- max(links) * 2 #conservative estimate on link strength
  bounding <- exp(max_degree * max_val)
  if (printProgress) cat('  Refining edge weights...\n')
  
  if (printProgress) cat('  |')
  links <- simplify2array(
    mclapply(1:nnodes, 
             function(i) {
               probs <- links[,i]
               mask <- 1:nnodes == i
               nonzeros <- (probs != 0) | mask #have to include the self val
               adjustment <- ifelse(i==1, 0, sum(!nonzeros[1:(i-1)]))
               
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
               if (printProgress && charsvec[i]) system('printf =')
               return(probs)
             }, 
             mc.cores=numCores, mc.preschedule=FALSE))  
  if (printProgress) cat('| ')
  
  for ( i in 1:(nrow(links)-1) ){
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
                        numCores=1, zero_cutoff=0, verbose=TRUE, ...){
  
  mult <- 1/niter
  intPA <- PAProfiles + 0
  intPA[intPA==0] <- -1
  nc <- ncol(intPA)
  pp <- F
  
  if (verbose & niter==1) pp <- T  
  else if (verbose){
    cat('Running DCA with', iter, 'iterations:\n')
    pb <- txtProgressBar(max=niter, style=3)
  } 
  
  truelinks <- countsmat <- matrix(0, nrow=nc, ncol=nc)
  for ( i in 1:niter ){
    #initlinks <- matrix(0, nrow=nc, ncol=nc)
    initlinks <- matrix(rnorm(nc**2), nrow=nc)
    iterlink <- DCA_logrise_run(intPA, initlinks, reg_const, 
                                printProgress=pp, numCores=numCores)
    countsmat <- countsmat + (iterlink != 0)
    truelinks <- truelinks + mult * iterlink
    if(verbose & !pp) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\n')
  
  if (zero_cutoff > 0) {
    cutoff <- ifelse(zero_cutoff > 1, zero_cutoff, zero_cutoff * niter)
    countsmat[countsmat < cutoff] <- 0
    countsmat[countsmat > 0] <- 1
    truelinks <- truelinks * countsmat
  } 
  
  return(truelinks)
}

getBuiltInEnsembleModel <- function(pw, flags){
  # Key: (val is binary + 1)
  # 000 => 1: Jaccard, Hamming, MI, ProfileDCA (base)
  # 001 => 2: base and MT, CT
  # 010 => 3: base and Behdenna
  # 011 => 4: base and Behdenna, MT, CT
  # 100 => 5: base and Coloc
  # 101 => 6: base and Coloc, MT, CT
  # 110 => 7: base and Behdenna, Coloc
  # 111 => 8: base and Behdenna, Coloc, MT, CT
  model <- 4 * flags[3] + 2 * flags[2] + 1 * flags[1] + 1
  load(builtinmodels)
  return(BuiltInEnsembles[[model]])
}
########


#### Class-Specific Method Definitions ####

PAProfiles.ProtWeaver <- function(pw, toEval=NULL, verbose=TRUE, specieslist=NULL){
  cols <- names(pw)
  ao <- attr(pw, 'allOrgs')
  if (!is.null(specieslist)){
    stopifnot('Species list is missing species!'=all(ao %in% specieslist))
    allOrgs <- specieslist
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(pw, 'useColoc')
  useMT <- attr(pw, 'useMT')
  if (useMT)
    pw <- lapply(pw, labels)
  if (useColoc)
    pw <- lapply(pw, gsub, pattern='(.+)_.+_[0-9]+', replacement='\\1')
  
  skip <- F
  if ( !is.null(toEval) ){
    skip <- T
    locs <- unique(c(toEval))
  }
  profiles <- matrix(F, nrow=length(allOrgs), ncol=length(pw))
  rownames(profiles) <- allOrgs
  colnames(profiles) <- cols
  if (verbose) pb <- txtProgressBar(max=length(pw), style=3)
  for ( i in 1:length(pw) ){
    if( !skip || i %in% locs)
      profiles[pw[[i]],i] <- T
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) cat('\n')
  if (!is.null(toEval))
    profiles <- profiles[,locs]
  return(profiles)
}

CophProfiles.ProtWeaver <- function(pw, toEval=NULL, verbose=TRUE){
  ## TODO: Some way to handle paralogs
  cols <- names(pw)
  allOrgs <- attr(pw, 'allOrgs')
  useColoc <- attr(pw, 'useColoc')
  useMT <- attr(pw, 'useMT')
  
  stopifnot('ProtWeaver object must be initialized with dendrograms to run MirrorTree methods'=
              useMT)

  skip <- F
  if ( !is.null(toEval) ){
    skip <- T
    locs <- unique(c(toEval))
  }
  l <- length(allOrgs)
  num_entries <- (l * (l-1)) / 2
  outmat <- matrix(0, nrow=num_entries, ncol=length(pw))
  dummycoph <- matrix(NA, nrow=l, ncol=l)
  ut <- upper.tri(dummycoph)
  rownames(dummycoph) <- colnames(dummycoph) <- allOrgs
  if (verbose) pb <- txtProgressBar(max=length(pw), style=3)
  for ( i in seq_along(pw) ){
    if ( !skip || i %in% locs ){
      dummycoph[] <- NA
      cop <- NA
      # This is occasionally throwing errors that don't affect output for some reason
      cop <- as.matrix(Cophenetic(pw[[i]]))
      copOrgNames <- rownames(cop)
      if (useColoc){
        copOrgNames <- sapply(copOrgNames, gsub, pattern='(.+)_.+_[0-9]+', replacement='\\1')
        rownames(cop) <- colnames(cop) <- copOrgNames
      }
      dummycoph[copOrgNames, copOrgNames] <- cop
      outmat[,i] <- dummycoph[ut]
    }
    if (verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\n')
  colnames(outmat) <- cols
  if (!is.null(toEval)){
    outmat <- outmat[,locs]
    ltr <- sapply(1:nrow(outmat), function(x) all(is.na(outmat[x,])))
    outmat <- outmat[!ltr,]
  }
  return(outmat)
}

MirrorTree.ProtWeaver <- function(pw, correction=c(),
                                  subset=NULL, verbose=TRUE,
                                  mySpeciesTree=NULL, precalcProfs=NULL, ...){
  pl <- length(pw)
  subs <- ProcessSubset(pw, subset)
  evalmap <- subs$evalmap
  uvals <- subs$uvals
  
  if (is.null(precalcProfs)){
    if (verbose) cat('Pre-processing distance matrices...\n')
    CPs <- CophProfiles(pw, uvals, verbose=verbose)
  } else {
    CPs <- precalcProfs
  }
  l <- ncol(CPs)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }
  
  correction <- tolower(correction)
  if ('speciestree' %in% correction){
    if (verbose) cat('Correcting with species tree...\n')
    stopifnot('Missing mySpeciesTree'=!is.null(mySpeciesTree))
    stopifnot('mySpeciesTree must be a dendrogram'=is(mySpeciesTree, 'dendrogram'))
    corrvec <- as.matrix(Cophenetic(mySpeciesTree))
    corrvec <- corrvec[upper.tri(corrvec)]
    stopifnot('mySpeciesTree has incorrect number of leaf nodes'=
                length(corrvec)==nrow(CPs))
    CPs <- CPs - corrvec
  }
  if ('normalize' %in% correction){
    if (verbose) cat('Normalizing profiles...\n')
    means <- colMeans(CPs, na.rm=T)
    vars <- apply(CPs, MARGIN=2, var, na.rm=T)
    for ( i in 1:ncol(CPs) ){
      CPs[,i] <- (CPs[,i] - means[i]) / (ifelse(vars[i]!=0, sqrt(vars), 1))
    }
  }
  if ('satoaverage' %in% correction){
    means <- rowMeans(CPs, na.rm = T)
    if (verbose) cat('Calculating Sato projection vectors...\n')
    
    # Big profiles lead to space issues that crash R
    if (nrow(CPs)**2 < (2**28)){
      if (verbose) pb <- txtProgressBar(max=ncol(CPs), style=3)
      proj_op <- diag(nrow=nrow(CPs)) - (means %*% t(means))
      for ( i in 1:ncol(CPs) ){
        CPs[,i] <- c(CPs[,i] %*% proj_op)
        if ( verbose ) setTxtProgressBar(pb, i)
      }
    } else {
      if (verbose) pb <- txtProgressBar(max=(ncol(CPs)*nrow(CPs)), style=3)
      for ( i in 1:ncol(CPs) ){
        v <- projv <- CPs[,i]
        multv <- means * v
        for ( j in 1:nrow(CPs) ){
          if ( !is.na(v[j]) )
            projv[j] <- sum(multv * means[j])
          if ( verbose ) setTxtProgressBar(pb, (i-1)*nrow(CPs) + j)
        }
        CPs[,i] <- v - projv
      }
    }
    if (verbose) cat('\n')
  }
  
  pairscores <- matrix(NA, nrow=pl, ncol=pl)
  ctr <- 1
  if (verbose) pb <- txtProgressBar(max=(pl*(pl-1) / 2), style=3)
  for ( i in 1:(pl-1) ){
    acc1 <- which(i == uvals)
    if (length(acc1) == 0) acc1 <- 1
    v1 <- CPs[,acc1]
    for ( j in (i+1):pl ){
      if (evalmap[i,j]){
        acc2 <- which( j == uvals )
        v2 <- CPs[,acc2]
        val <- cor(v1, v2, use='na.or.complete', method='pearson')
        pairscores[i,j] <- pairscores[j,i] <- ifelse(is.na(val), 0, val)
      }
      ctr <- ctr + 1
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  for ( i in uvals )
    pairscores[i,i] <- 1
  if ('partialcorrelation' %in% correction){
    flag <- T
    if (!is.null(subset)){
      opsm <- pairscores
      pairscores <- pairscores[uvals, uvals]
      if (any(is.na(pairscores))){
        pairscores <- opsm
        warning('Partial correlation requires a square matrix. Skipping.')
        flag <- F
      }
    }
    if (flag){
      d <- det(pairscores)
      if (d == 0){
        warning('Matrix is exactly singular, cannot use partial correlation correction.')
      } else {
        inv <- solve(pairscores)
        cols <- matrix(diag(inv), nrow=nrow(inv), ncol=ncol(inv))
        rows <- matrix(diag(inv), nrow=nrow(inv), ncol=ncol(inv), byrow=T)
        divisor <- sqrt(cols * rows)
        pairscores <- (-1 * inv) / divisor
      }
      if ( !is.null(subset) ){
        opsm[uvals,uvals] <- pairscores
        pairscores <- opsm
      }
    }
  }
  for ( i in uvals )
    pairscores[i,i] <- 1
  rownames(pairscores) <- colnames(pairscores) <- names(pw)
  return(pairscores[uvals, uvals])
}

ContextTree.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, precalcProfs=NULL, 
                                   mySpeciesTree=NULL, ...){
  
  if ( !is.null(mySpeciesTree) && is(mySpeciesTree, 'dendrogram')){
    correction <- c('speciestree', 'normalize', 'partialcorrelation')
  } else { 
    correction <- c('normalize', 
                    'satoaverage', 
                    'partialcorrelation')
  }
  
  return(MirrorTree(pw, correction=correction,
                    verbose=verbose, 
                    precalcCProfs=precalcProfs,
                    mySpeciesTree=mySpeciesTree))
}

Jaccard.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE,
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap

  if ( is.null(precalcProfs) ){
    if (verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, verbose=verbose)
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
  ctr <- 1
  if (verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in 1:l ){
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
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  diag(pairscores) <- 0
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  pairscores <- 1 - pairscores # because distance
  return(pairscores)
}

Hamming.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, 
                               precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, verbose=verbose)
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
  ctr <- 1
  if (verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in 1:l ){
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
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  diag(pairscores) <- 0
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  mp <- max(pairscores, na.rm=T)
  pairscores <- (mp - pairscores) / mp #because distance
  return(pairscores)
}

MutualInformation.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, 
                                         precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if ( is.null(precalcProfs) ){
    if (verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, verbose=verbose)
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
  ctr <- 1
  if (verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in 1:l ){
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
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  n <- n[uvals]
  rownames(pairscores) <- colnames(pairscores) <- n
  
  # Correction
  denom <- mean(pairscores[upper.tri(pairscores)], na.rm=T)
  pairscores <- pairscores / ifelse(denom==0, 1, denom)
  
  # Normalize
  denom <- max(pairscores, na.rm=T)
  pairscores <- pairscores / ifelse(denom==0, 1, denom)
  diag(pairscores) <- 1
  #pairscores <- pairscores #because distance
  return(pairscores)
}

ProfileDCA.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, numCores=1,
                           precalcProfs=NULL, precalcSubset=NULL, useAbs=T, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  
  if ( is.null(precalcProfs) ){
    if (verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, verbose=verbose)
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
  
  pairscores <- DCA_logRISE(pap, verbose=verbose, numCores=numCores, ...)
  rownames(pairscores) <- colnames(pairscores) <- n
  if (useAbs) pairscores <- abs(pairscores)
  if (max(pairscores) != 0)
    pairscores <- pairscores / max(abs(pairscores))
  
  
  return(pairscores)
}

Coloc.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, 
                             precalcProfs=NULL, precalcSubset=NULL, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
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
  if (verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in 1:l){
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
        mult <- ifelse(length(shared)==0, 1, 1/length(shared))
        for ( k in seq_along(shared) ){
          spk <- shared[k]
          v1 <- lab1[l1sp==spk]
          v2 <- lab2[l2sp==spk]
          idx1 <- as.integer(gsub('.*_([0-9]*)$', '\\1', v1))
          idx2 <- as.integer(gsub('.*_([0-9]*)$', '\\1', v2))
          vals <- expand.grid(idx1, idx2)
          diffs <- exp(1 - abs(vals[,1] - vals[,2]))
          score <- score + sum(diffs) / nrow(vals)
        }
        score <- score * mult
        pairscores[i,j] <- pairscores[j,i] <- score
      }
      ctr <- ctr + 1
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  m <- ifelse(max(pairscores, na.rm=TRUE) != 0, max(pairscores,na.rm=TRUE), 1)
  pairscores <- pairscores / m
  diag(pairscores) <- 1
  rownames(pairscores) <- colnames(pairscores) <- n
  return(pairscores)
}

Behdenna.ProtWeaver <- function(pw, subset=NULL, verbose=TRUE, 
                                mySpeciesTree=NULL, useSubtree=FALSE, 
                                useACCTRAN=TRUE, rawZScores=FALSE, 
                                precalcProfs=NULL, precalcSubset=NULL, ...){
  stopifnot('No species tree provided.'=(!is.null(mySpeciesTree)))
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  n <- names(pw)[uvals]
  if ( is.null(precalcProfs) ){
    if (verbose) cat('Calculating PA Profiles...\n')
    pap <- PAProfiles(pw, uvals, verbose=verbose, specieslist=labels(mySpeciesTree))
  } else {
    pap <- precalcProfs
  }
  l <- length(uvals)
  stopifnot(nrow(pap) > 1)
  fd <- FastDend(mySpeciesTree)
  v1 <- abs(generateGainLossVec(fd, pap[,1], moveEventsUpward=useACCTRAN))
  glmat <- matrix(0, nrow=length(v1), ncol=ncol(pap))
  glmat[,1] <- v1
  if (verbose) cat('  Calculating gain/loss vectors:\n')
  if (verbose) pb <- txtProgressBar(max=ncol(pap), style=3)
  for ( i in 2:ncol(pap) ){
    glmat[,i] <- abs(generateGainLossVec(fd, pap[,i], moveEventsUpward=useACCTRAN))
    if (verbose) setTxtProgressBar(pb, i)
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
  if (verbose) cat('\n  Calculating pairscores:\n')
  if (verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in 1:l){
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
      if (verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (verbose) cat('\n')
  diag(pairscores) <- 0 #because z-scores
  
  if (!rawZScores){
    pairscores <- abs(pairscores)
    pairscores <- pairscores / ifelse(max(pairscores,na.rm=T) != 0, 
                                      max(pairscores, na.rm=T), 1)
    diag(pairscores) <- 1
  }
  rownames(pairscores) <- colnames(pairscores) <- n
  return(pairscores)
}

Ensemble.ProtWeaver <- function(pw,
                                subset=NULL, verbose=TRUE, mySpeciesTree=NULL,
                                pretrainedmodel=NULL,
                                noPrediction=FALSE, ...){
  
  flags <- rep(F, 3)
  
  subs <- ProcessSubset(pw, subset)
  uvals <- subs$uvals
  splist <- NULL
  if (!is.null(mySpeciesTree)){
    splist <- labels(mySpeciesTree)
  }
    
  if (verbose) cat('Calculating P/A profiles:\n')
  PAs <- PAProfiles(pw, uvals, verbose=verbose, specieslist=splist)
  CPs <- NULL
  takesCP <- c('ContextTree', 'MirrorTree')
  
  submodels <- c('ProfileDCA', 'Jaccard', 'Hamming', 'MutualInformation')
  if (attr(pw, 'useMT')){
    flags[1] <- T
    if (verbose) cat('Calculating Cophenetic profiles:\n')
    CPs <- CophProfiles(pw, uvals, verbose=verbose)
    submodels <- c(submodels, takesCP)
  }
  
  if (!is.null(mySpeciesTree)){
    flags[2] <- T
    submodels <- c(submodels, 'Behdenna')
  }

  if (attr(pw, 'useColoc')){
    flags[3] <- T
    submodels <- c(submodels, 'Coloc')
  }
  
  if(!is.null(pretrainedmodel)) {
    predictionmodel <- pretrainedmodel
  } else {
    predictionmodel <- getBuiltInEnsembleModel(pw, flags)
  }
  
  results <- list()
  for ( model in submodels ){
    if (verbose) cat('Running ', model, ':\n', sep='')
    if (model %in% takesCP) profs <- CPs
    else profs <- PAs
    results[[model]] <- predict(pw, model, verbose=verbose, 
                              ensembleprediction=T, precalcProfs=profs,
                              precalcSubset=subs, 
                              mySpeciesTree=mySpeciesTree, ...)
  }
  
  cat('Calculating additional P/A Statistics...\n')
  results <- AdjMatToDf(results)
  pas <- PAStats(results, PAs) 
  results[,'AvgOcc'] <- pas$avg
  results[,'OccDiff'] <- pas$diff
  # if (is.null(trainingset)){
  #   
  # }
  if (noPrediction) return(list(res=results, noPostFormatting=TRUE))
  
  if (verbose) cat('Predicting with Ensemble method...\n')
  predictions <- predict(predictionmodel, results[,-c(1,2)])
  outmat <- matrix(NA, nrow=length(uvals), ncol=length(uvals))
  for (i in seq_along(predictions)){
    i1 <- which(results[i,1] == uvals)
    i2 <- which(results[i,2] == uvals)
    pred <- predictions[i]
    outmat[i1,i2] <- outmat[i2,i1] <- pred
  }
  rownames(outmat) <- colnames(outmat) <- uvals
  
  minom <- -1*min(outmat, na.rm=TRUE)
  maxom <- max(outmat, na.rm=TRUE) + minom
  maxom <- ifelse(maxom==0, 1, maxom)
  
  outmat <- (outmat + minom) / maxom

  return(outmat)
}

########