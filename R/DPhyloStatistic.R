DPhyloStatistic <- function(dend, PAProfile, NumIter=1000L){
  # Error check simple variables
  stopifnot("'dend' must be of type 'dendrogram'"=is(dend, 'dendrogram'))
  if(is(NumIter, 'numeric')) NumIter <- as.integer(NumIter)
  stopifnot("NumIter must be a positive integer!"=is.integer(NumIter) && NumIter > 0)
  
  allLabels <- labels(dend)
  if (!is.character(allLabels)){
    dend <- dendrapply(dend, \(x){
      if (!is.null(attr(x, 'leaf'))){
        attr(x, 'label') <- as.character(attr(x, 'label'))
      }
      return(x)
    })
  }
  allLabels <- labels(dend)
  # Error check PAProfile
  if(is(PAProfile, 'integer') || is(PAProfile, 'numeric')){
    stopifnot("PAProfile length not equal to number of labels!"=
                length(PAProfile) == length(allLabels))
    if (any(PAProfile != 1 && PAProfile != 0)){
      stop("PAProfile has values not in (0,1)!")
    }
    PAProfile <- allLabels[PAProfile==1]
  } else if (is(PAProfile, 'character')){
    if (any(vapply(PAProfile, \(x) !(x %in% allLabels), logical(1)))){
      warning("PAProfile contains labels not in dendrogram!")
      PAProfile <- intersect(allLabels, PAProfile)
    }
  } else if (is(PAProfile, 'logical')){
    stopifnot("PAProfile length not equal to number of labels!"=
              length(PAProfile) == length(allLabels))
    PAProfile <- allLabels[PAProfile]
  } else {
    stop("PAProfile must be a vector of character, integer, or logical.")
  }
  
  if (length(PAProfile) == 0){
    warning("No elements present, returning 0")
    return(0)
  }
  
  y <- .Call("initCDend", dend)
  on.exit(rm(y))
  
  Dobs <- .Call('calcDValue', y, PAProfile)
  Dr <- .Call('calcDRandValue', y, allLabels, length(PAProfile), NumIter)
  Db <- .Call('calcDBrownValue', y, allLabels, NumIter, length(PAProfile) / length(allLabels), 0.5, length(PAProfile) / length(allLabels))
  if (Db - Dr == 0){
    warning("Denominator is zero!")
  }

  Dstatistic <- (Dobs - Db) / (Dr - Db)
  return(Dstatistic)
}