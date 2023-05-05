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
Ensemble <- function(pw, ...) UseMethod('Ensemble')
SpeciesTree <- function(pw, Verbose, Processors) UseMethod('SpeciesTree')
########


#### Class Constructor ####
new_ProtWeaver <- function(validatedInput){
  structure(validatedInput$ipt,
            allOrgs=validatedInput$allgenomes,
            useMT=validatedInput$flags$usemirrortree,
            useColoc=validatedInput$flags$usecoloc,
            useResidue=validatedInput$flags$useresidue,
            useStrand=validatedInput$flags$strandid,
            speciesTree=validatedInput$speciestree,
            class='ProtWeaver')
}

ProtWeaver <- function(ListOfData, MySpeciesTree=NULL, NoWarn=FALSE){
  stopifnot("MySpeciesTree should be NULL or an object of type 'dendrogram'"=
              is.null(MySpeciesTree) || is(MySpeciesTree,'dendrogram'))
  vRes <- validate_ProtWeaver(ListOfData, noWarn=NoWarn)
  if(!is.null(MySpeciesTree) && any(!(vRes$allgenomes %in% labels(MySpeciesTree)))){
    stop("MySpeciesTree is missing labels!")
  }
  if(!is.null(MySpeciesTree))
    vRes$allgenomes <- labels(MySpeciesTree)
  vRes$speciestree <- MySpeciesTree
  new_ProtWeaver(vRes)
}

validate_ProtWeaver <- function(ipt, noWarn=FALSE){
  bitflags <- list(usecoloc=FALSE, usemirrortree=FALSE, useresidue=FALSE)
  stopifnot('ProtWeaver expects a list of dendrograms or character vectors as input.'=
              is(ipt, 'list'))
  
  ipt <- ipt[!vapply(ipt, is.null, TRUE)]
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
    allentries <- character(0)
    for(tree in ipt){
      allentries <- unique(c(allentries, as.character(labels(tree))))
    }
  }
  
  if (bitflags[['usemirrortree']]){
    useResidueMI <- TRUE
    for ( tree in ipt ){
      useResidueMI <- dendrapply(tree, 
                                 \(x){
                                      sv <- !is.null(attr(x,'state'))
                                      if(is.leaf(x)) 
                                        return(sv)
                                      return(all(sv, unlist(x)))
                                    }, how='post.order')
      if (!useResidueMI){
        if (!noWarn) message('Disabling Residue methods. Input dendrograms must',
                             ' include ancenstral state reconstruction for residue',
                             ' methods. Consult the documentation for more info.\n')
        break
      }
    }
    bitflags[['useresidue']] <- useResidueMI
  }
  
  checkforcoloc <- grepl('[^_]+_[^_]+_.*[0-9]+', allentries)
  if ( all(checkforcoloc) ){
    bitflags[['usecoloc']] <- TRUE
    checkforstrand <- grepl('[^_]+_[^_]+_[01]_[0-9]+', allentries)
    bitflags[['strandid']] <- all(checkforstrand)
    allentries <- unique(gsub('([^_]+)_[^_]+_.*[0-9]+', '\\1', allentries))
  }
  else{
    bitflags[['usecoloc']] <- FALSE
    bitflags[['strandid']] <- FALSE
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
show.ProtWeaver <- function(object){
  if (length(object) == 1){
    cat(paste('a ProtWeaver object with', length(object),
              'group and', length(attr(object,'allOrgs')), 'genomes.\n'))
  } else {
    cat(paste('a ProtWeaver object with', length(object),
              'groups and', length(attr(object,'allOrgs')), 'genomes.\n'))
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

predict.ProtWeaver <- function(object, Method='Ensemble', Subset=NULL, Processors=1L,
                               MySpeciesTree=SpeciesTree(object), 
                               PretrainedModel=NULL,
                               RawZScores=FALSE, NoPrediction=FALSE, 
                               ReturnRawData=FALSE, Verbose=TRUE, ...){
  pw <- object
  func <- getS3method(Method, 'ProtWeaver')
  if(Verbose && !ReturnRawData) starttime <- Sys.time()
  
  preds <- func(pw, Subset=Subset, Verbose=Verbose, 
                MySpeciesTree=MySpeciesTree, Processors=Processors,
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
  if (Method=='TreeDistance'){
    pnames <- names(preds)
    for (i in seq_along(pnames)){
      names(preds[[i]]) <- n
      preds[[i]] <- structure(preds[[i]],
                              method=pnames[i],
                              class=c('ProtWeb', 'simMat'))
    }
    rs <- preds
    if (length(rs) == 1) rs <- rs[[1]]
  } else {
    names(preds) <- n
    rs <- structure(preds,
                    method=Method,
                    class=c('ProtWeb', 'simMat'))
  }
  
  rs
}

SpeciesTree.ProtWeaver <- function(pw, Verbose=TRUE, Processors=1L){
  tree <- attr(pw,'speciesTree')
  if(is.null(tree) && attr(pw, 'useMT'))
    tree <- findSpeciesTree(pw, Verbose=Verbose, Processors=Processors)
  
  tree
}

########

#### Internal S3 Methods ####

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
  
  cat('Combining results...\n')
  results <- AdjMatToDf(results, Verbose=Verbose)
  cat('Calculating additional P/A Statistics...\n')
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
  outmat <- simMat(NA_real_, length(unames), NAMES=unames)
  if (Verbose) pb <- txtProgressBar(max=length(predictions), style=3)
  for (i in seq_along(predictions)){
    i1 <- which(results[i,1] == unames)
    i2 <- which(results[i,2] == unames)
    pred <- predictions[i]
    outmat[i1, i2] <- pred
    if (Verbose) setTxtProgressBar(pb, i)
  }
  cat('\n')
  Diag(outmat) <- 1
  return(outmat)
}

########