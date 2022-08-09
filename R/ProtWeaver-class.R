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

#### IMPORTS ####
# Relies on functions in the following files:
#   - ProtWeaver-DMPreds.R          ( Distance Matrix Predictors )
#   - ProtWeaver-PAPreds.R          ( Pres/Abs Predictors)
#   - ProtWeaver-ColocPreds.R       ( Colocalization Predictors )
#   - ProtWeaver-ResiduePreds.R     ( Residue level Predictors )
#   - ProtWeaver-utils.R            ( Helper functions )
#################


#### S3 Generic Definitions ####
PAProfiles <- function(pw, ...) UseMethod('PAProfiles')
CophProfiles <- function(pw, ...) UseMethod('CophProfiles')
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