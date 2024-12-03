##### -- EvoWeaver Class to find clusters of functionally linked genes -------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu
#
# TODO: ensure dendrogram branch lengths are numeric and not integer

#### DEFINITION ####
# Class expects as input one of the following:
#   1. A list of dendrograms
#      a) Dendrograms with labels [ORG]_[INDEX]_[POS] will use co-loc code
#      b) Dendrograms without will not use co-loc
#   2. A list of character vectors defining COGs (will not use RPMirrorTree)
#   3. A list of character vectors defining COGs with labels like 1(a)
####################

#### IMPORTS ####
# Relies on functions in the following files:
#   - EvoWeaver-PSPreds.R          ( Phylogenetic Structure Predictors )
#   - EvoWeaver-PPPreds.R          ( Phylogenetic Profiling Predictors)
#   - EvoWeaver-GOPreds.R          ( Gene Organization Predictors )
#   - EvoWeaver-SLPreds.R          ( Sequence-Level Predictors )
#   - EvoWeaver-utils.R            ( Helper functions )
#################


#### S3 Generic Definitions ####
Ensemble <- function(ew, ...) UseMethod('Ensemble')
SpeciesTree <- function(ew, Verbose, Processors) UseMethod('SpeciesTree')
########


#### Class Constructor ####
new_EvoWeaver <- function(validatedInput){
  structure(validatedInput$ipt,
            allOrgs=validatedInput$allgenomes,
            useMT=validatedInput$flags$usemirrortree,
            useColoc=validatedInput$flags$usecoloc,
            useResidue=validatedInput$flags$useresidue,
            useStrand=validatedInput$flags$strandid,
            speciesTree=validatedInput$speciestree,
            class='EvoWeaver')
}

EvoWeaver <- function(ListOfData, MySpeciesTree=NULL, NoWarn=FALSE){
  stopifnot("MySpeciesTree should be NULL or an object of type 'dendrogram'"=
              is.null(MySpeciesTree) || is(MySpeciesTree,'dendrogram'))
  vRes <- validate_EvoWeaver(ListOfData, noWarn=NoWarn)
  if(!is.null(MySpeciesTree) && any(!(vRes$allgenomes %in% labels(MySpeciesTree)))){
    stop("MySpeciesTree is missing labels!")
  }
  if(!is.null(MySpeciesTree))
    vRes$allgenomes <- labels(MySpeciesTree)
  vRes$speciestree <- MySpeciesTree
  new_EvoWeaver(vRes)
}

validate_EvoWeaver <- function(ipt, noWarn=FALSE){
  bitflags <- list(usecoloc=FALSE, usemirrortree=FALSE, useresidue=FALSE)
  stopifnot('EvoWeaver expects a list of dendrograms or character vectors as input.'=
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
    for(i in seq_along(ipt)){
      tree <- ipt[[i]]
      allentries <- unique(c(allentries, as.character(labels(tree))))

      # heights of trees must be present and numeric, if not then correct it
      if(is.null(attr(tree, 'height'))){
        stop("Input tree at index ", i, " is missing a height!")
      }
      if(is.integer(attr(tree, 'height')) || any(!is.character(rapply(tree, \(x) attr(x, 'label'))))){
        ipt[[i]] <- dendrapply(tree, \(x){
          if(is.leaf(x)) attr(x, 'label') <- as.character(attr(x, 'label'))
          attr(x, 'height') <- as.numeric(attr(x, 'height'))
          x
        })
      }
    }
  }

  if (bitflags[['usemirrortree']]){
    useResidueMI <- TRUE
    for ( tree in ipt ){
      # useResidueMI <- dendrapply(tree,
      #                            \(x){
      #                                 sv <- !is.null(attr(x,'state'))
      #                                 if(is.leaf(x))
      #                                   return(sv)
      #                                 return(all(sv, unlist(x)))
      #                               }, how='post.order')
      useResidueMI <- all(rapply(tree, \(x) !is.null(attr(x, 'state'))))
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
                         '[GENOME]_[CONTIG]_[INDEX]_[ORDER] to use co-localization ',
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

validate_EvoWeaver_methods <- function(methodnames){
  if(!is.character(methodnames)){
    stop('Argument "Method" requires input of type "character"')
  }
  pp_preds <- c("ExtantJaccard","Hamming","CorrGL","GLMI",
    "ProfileDCA","Behdenna","GLDistance","PAJaccard",
    "PAOverlap","PAPV")

  ps_preds <- c("RPMirrorTree", "RPContextTree", "TreeDistance")

  go_preds <- c("GeneDistance", "MoransI", "OrientationMI")

  sl_preds <- c("SequenceInfo", "GeneVector", "Ancestral")


  all_predictors <- c(pp_preds, ps_preds, go_preds, sl_preds, "Ensemble")
  meta_categories <- c("PhylogeneticProfiling", "PhylogeneticStructure",
    "GeneOrganization", "SequenceLevel")
    ## keyword inputs subset to the method categories used in manuscript\
  if("all" %in% methodnames){
    methodnames <- c(methodnames, meta_categories)
    methodnames <- methodnames[methodnames != 'all']
  }
  if("PhylogeneticProfiling" %in% methodnames)
    methodnames <- c(methodnames, pp_preds[c(4,7:9)])
  if ("PhylogeneticStructure" %in% methodnames)
    methodnames <- c(methodnames, ps_preds)
  if ("GeneOrganization" %in% methodnames)
    methodnames <- c(methodnames, go_preds)
  if ("SequenceLevel" %in% methodnames)
    methodnames <- c(methodnames, sl_preds[seq_len(2L)])
  methodnames <- unique(methodnames)
  methodnames <- methodnames[!methodnames %in% meta_categories]

  p_notin <- match(methodnames, all_predictors, nomatch=0) == 0
  if(any(p_notin)){
    unrecognized_methods <- methodnames[p_notin]
    unrecognized_methods <- paste(unrecognized_methods, collapse='", "')
    unrecognized_methods <- paste0('"', unrecognized_methods, '"')
    stop("Invalid methods provided to 'Method' argument. The following are unrecognized: ", unrecognized_methods)
  }

  ## Things for precalculation
  precalc_pa <- FALSE
  if(sum(match(methodnames, pp_preds, nomatch=0) > 0)){
    precalc_pa <- TRUE
  }

  list(Method=methodnames, precalcPA=precalc_pa)
}

Standardize_Subset <- function(Subset, ew){
  n <- names(ew)
  if(is.matrix(Subset)){
    Subset <- as.data.frame(Subset)
  }
  if(is.data.frame(Subset)){
    if(ncol(Subset) != 2)
      stop("If 'Subset' is of type 'data.frame', it must have exactly two columns.")
  } else if (is.character(Subset)){
    v <- unique(Subset)
    vp <- match(v, n, nomatch=0)
    if(any(vp==0)){
      pos_missing <- v[vp==0]
      line <- paste(pos_missing, collapse='", "')
      line <- paste0('"', line, '"')
      stop("Invalid gene group names in 'Subset'. The following were not found: ", line)
    }
    Subset <- vp
  } else if (is.numeric(Subset)){
    v <- unique(Subset)
    vp <- as.integer(v)
    if(any(vp != v))
      warning("Some values in Subset were changed during conversion from numeric to integer.")
    if(any(vp > length(n) | vp < 1)){
      stop("'Subset' has values that are out of range for the EvoWeaver object.")
    }
    Subset <- vp
  } else {
    stop("'Subset' must be either an integer vector, a character vector, or a data.frame/matrix with two columns.")
  }

  if(!is.data.frame(Subset)){
    ## by now, Subset must be either a data.frame or a 1D vector
    ## in the latter case, it must be an integer
    l <- length(n)
    vp <- sort(Subset)
    l2 <- length(vp)
    ## first pair has l-1 pairs, second has l-2, etc...
    ## this equals sum_{i in [1:l2]}(l - i)
    ## = sum(l) - sum(i)
    ## = l*l2 - TRI(l2), where TRI(n) the n'th triangular number
    rep_times <- rep(l, l2) - seq_len(l2)
    total_pairs <- rep_times
    v1 <- rep(vp, times=rep_times)
    v2 <- integer(sum(rep_times))
    all_v <- seq_len(l)
    rep_times <- cumsum(c(0,rep_times))
    for(i in seq_len(l2))
      v2[seq(rep_times[i]+1, rep_times[i+1])] <- all_v[-vp[seq_len(i)]]
    Subset <- data.frame(Gene1=n[v1], Gene2=n[v2])
  } else {
    for(i in seq_len(2)){
      if(is.numeric(Subset[,i])){
        vp <- as.integer(Subset[,i])
        if(any(vp != Subset[,i])){
          warning("Some values in Subset were changed during conversion from numeric to integer.")
        }
        Subset[,i] <- vp
      }
      if(is.integer(Subset[,i])){
        v <- unique(unlist(Subset[,i]))
        if(any(v > length(n) | v < 1)){
          stop("Column ", i, " has integer values that are out of range for the EvoWeaver object.")
        }
        Subset[,i] <- n[Subset[,i]]
      } else if(is.character(Subset[,i])){
        v <- unique(unlist(Subset[,i]))
        vp <- match(v, n, nomatch=0)
        if(any(v==0)){
          pos_missing <- v[vp==0]
          line <- paste(pos_missing, collapse='", "')
          line <- paste0('"', line, '"')
          stop("Invalid gene group names in 'Subset'. The following were not found: ", line)
        }
      }
    }
  colnames(Subset) <- c("Gene1", "Gene2")
  }
  Subset
}

########

#### User-Exposed S3 Methods ####
show.EvoWeaver <- function(object){
  if (length(object) == 1){
    cat(paste('a EvoWeaver object with', length(object),
              'group and', length(attr(object,'allOrgs')), 'genomes.\n'))
  } else {
    cat(paste('a EvoWeaver object with', length(object),
              'groups and', length(attr(object,'allOrgs')), 'genomes.\n'))
  }
}

print.EvoWeaver <- function(x, ...){
  if (length(x) == 1){
    cat(paste('a EvoWeaver object with', length(x),
              'group and', length(attr(x,'allOrgs')), 'genomes.\n'))
  } else {
    cat(paste('a EvoWeaver object with', length(x),
              'groups and', length(attr(x,'allOrgs')), 'genomes.\n'))
  }
}

`[.EvoWeaver` <- function(x, i){
  y <- unclass(x)
  newv <- validate_EvoWeaver(y[i], noWarn=TRUE)
  newv$speciestree <- attr(x, 'speciesTree')
  if(!is.null(newv$speciestree))
    newv$allgenomes <- labels(newv$speciestree)
  new_EvoWeaver(newv)
}

predict.EvoWeaver <- function(object, Method='Ensemble', Subset=NULL, Processors=1L,
                               MySpeciesTree=SpeciesTree(object, Verbose=Verbose),
                               PretrainedModel="KEGG",
                               NoPrediction=FALSE,
                               ReturnDataFrame=TRUE,
                               Verbose=interactive(),
                               CombinePVal=TRUE,...){
  ew <- object
  if(!is.null(Subset))
    Subset <- Standardize_Subset(Subset, ew)
  USEENSEMBLE <- FALSE
  if(is.character(Method) && "Ensemble" %in% Method){
    if(length(Method) > 1) stop("Method='Ensemble' cannot be accompanied by other methods.")
    USEENSEMBLE <- TRUE
    Method <- c("PhylogeneticProfiling")
    if(attr(ew, "useMT"))
      Method <- c(Method, "PhylogeneticStructure")
    if(attr(ew, "useColoc"))
      Method <- c(Method, "GeneOrganization")
    if(attr(ew, "useResidue"))
      Method <- c(Method, "SequenceLevel")
  }
  validatedMethod <- validate_EvoWeaver_methods(Method)
  precalc_pa <- validatedMethod$precalcPA
  Method <- validatedMethod$Method

  subs <- ProcessSubset(ew, Subset)
  n <- names(ew)
  uvals <- subs$uvals
  unames <- vapply(uvals, function(x) n[x], FUN.VALUE=character(1))
  splist <- NULL
  if (!is.null(MySpeciesTree)){
    splist <- labels(MySpeciesTree)
  }

  if(precalc_pa){
    if (Verbose) cat('Calculating P/A profiles:\n')
    PAs <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=splist)
  } else {
    PAs <- NULL
  }

  multiplepredictors <- length(Method)!=1
  methodnames <- character(0L)
  lst <- list()
  ctr <- 1L
  for(methodtype in Method){
    if(Verbose) cat("\nMethod '", methodtype, "':\n", sep='')
    func <- getS3method(methodtype, 'EvoWeaver')
    if(Verbose) starttime <- Sys.time()

    preds <- func(ew, Subset=Subset, Verbose=Verbose,
                  MySpeciesTree=MySpeciesTree, Processors=Processors,
                  PretrainedModel=PretrainedModel,
                  NoPrediction=NoPrediction,
                  precalcSubset=subs,
                  precalcProfs=PAs,
                  CombinePVal=CombinePVal,...)

    if (Verbose){
      cat('Done.\nTime difference of',
          round(difftime(Sys.time(), starttime, units = 'secs'), 2),
          'seconds.\n\n')
    }
    if (is(preds, 'list') && !is.null(preds$noPostFormatting))
      return(invisible(preds$res))

    if (methodtype=='TreeDistance'){
      pnames <- names(preds)
      for (i in seq_along(pnames)){
        #names(preds[[i]]) <- n
        preds[[i]] <- structure(preds[[i]],
                                method=pnames[i],
                                class=c('EvoWeb', 'simMat'))
      }
      rs <- preds
      for(treemethod in pnames){
        lst[[ctr]] <- rs[[treemethod]]
        ctr <- ctr + 1L
      }
      if (length(rs) == 1) rs <- rs[[1]]
      methodnames <- c(methodnames, pnames)
    } else {
      #names(preds) <- n
      rs <- structure(preds,
                      method=methodtype,
                      class=c('EvoWeb', 'simMat'))
      methodnames <- c(methodnames, methodtype)
      lst[[ctr]] <- rs
      ctr <- ctr + 1L
    }
  }
  names(lst) <- methodnames

  if(ReturnDataFrame || USEENSEMBLE){
    if(Verbose) cat("Building Dataframe:\n")
    lst <- AdjMatToDf(lst, Verbose=Verbose, Subset=Subset, CombinePVal)
    if(Verbose) cat("Done.\n\n")
  }

  if(USEENSEMBLE && !NoPrediction){
    lst <- cbind(lst, Ensemble=EvoWeaverEnsemblePrediction(lst, PretrainedModel, CombinePVal))
  }
  lst
}

SpeciesTree.EvoWeaver <- function(ew, Verbose=TRUE, Processors=1L){
  tree <- attr(ew,'speciesTree')
  if(is.null(tree) && attr(ew, 'useMT'))
    tree <- findSpeciesTree(ew, Verbose=Verbose, Processors=Processors)

  tree
}

########

#### Internal S3 Methods ####

EvoWeaverEnsemblePrediction <- function(preds, PretrainedModel="KEGG", CombinePVal){
  if(is.null(PretrainedModel)){
    stop("No model provided to 'PretrainedModel'. Try specifying PretrainedModel=\"KEGG\" or \"CORUM\"")
  }
  if(!CombinePVal){
    pn <- vapply(strsplit(colnames(preds), '.', fixed=TRUE), .subset, character(1L), 1L)
    pn_more <- table(pn) > 1
    pn_keep <- names(pn_more)[!pn_more]
    opreds <- preds[,pn%in%pn_keep]
    colnames(opreds) <- pn[pn%in%pn_keep]
    pn_more <- names(pn_more)[pn_more]
    for(n in pn_more){
      p <- which(n == pn)
      opreds[[n]] <- do.call(`*`, preds[,p])
    }
    preds <- opreds
  }
  if (is.character(PretrainedModel)){
    predictions <- predictWithBuiltins(preds, PretrainedModel)
  } else if(is(PretrainedModel, 'glm')) {
    predictions <- predict(PretrainedModel, preds, type='response')
  } else {
    predictions <- predict(PretrainedModel, preds)
  }
  predictions
}

OldEnsemble.EvoWeaver <- function(ew,
                                Subset=NULL, Verbose=TRUE, MySpeciesTree=NULL,
                                PretrainedModel=NULL,
                                NoPrediction=FALSE, ...){

  flags <- rep(FALSE, 3)

  subs <- ProcessSubset(ew, Subset)
  n <- names(ew)
  uvals <- subs$uvals
  unames <- vapply(uvals, function(x) n[x], FUN.VALUE=character(1))
  splist <- NULL
  if (!is.null(MySpeciesTree)){
    splist <- labels(MySpeciesTree)
  }

  if (Verbose) cat('Calculating P/A profiles:\n')
  PAs <- PAProfiles(ew, uvals, Verbose=Verbose, speciesList=splist)
  CPs <- NULL
  takesCP <- c('RPMirrorTree') # Just using MirrorTree for prediction

  submodels <- c('ProfileDCA', 'ExtantJaccard', 'Hamming', 'GLMI')
  if (attr(ew, 'useMT')){
    if (Verbose) cat('Calculating Cophenetic profiles:\n')
    CPs <- CophProfiles(ew, uvals, Verbose=Verbose, speciesList=splist)
    submodels <- c(submodels, takesCP)
  }

  if (!is.null(MySpeciesTree) && CheckBifurcating(MySpeciesTree)){
    submodels <- c(submodels, 'Behdenna')
  }

  if (attr(ew, 'useColoc')){
    submodels <- c(submodels, 'GeneDistance')
  }

  if (attr(ew, 'useResidue')){
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
    results[[model]] <- predict(ew, model, Verbose=Verbose,
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
  if(Verbose) cat('\n')
  Diag(outmat) <- 1
  return(outmat)
}

########
