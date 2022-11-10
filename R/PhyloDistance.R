GeneralizedRF <- function(val, RawScore=FALSE){
  if (RawScore){
    retval <- val
    names(retval) <- c("Similarity", "dend1.Entropy", "dend2.Entropy")
    if (val[1] == 0) retval[1:2] <- c(NA, NA)
    return(retval)
  }
  
  if (val[1] == 0) return(1)
  maxval <- 0.5*(val[2] + val[3])
  retval <- (maxval - val[1]) / maxval
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

RFDist <- function(val, RawScore=FALSE){
  if (RawScore){
    retval <- val
    names(retval) <- c("UniqueSplits", "dend1.Splits", "dend2.Splits")
    return(retval)
  }
  
  maxval <- val[2] + val[3]
  retval <- val[1] / maxval
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

KFDist <- function(val){
  maxval <- val[2]
  retval <- val[1] / maxval
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

JRFDist <- function(val, RawScore=FALSE){
  if (RawScore){
    retval <- val
    names(retval) <- c("Distance", "dend1.NumSplits", "dend2.NumSplits")
    if (val[1] == 0) retval[1:2] <- c(NA, NA)
    return(retval)
  }
  if (val[1] == 0) return(1)
  maxval <- (val[2] + val[3])
  retval <- 1 - (val[1] / maxval)
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

PhyloDistance <- function(dend1, dend2, Method="GRF", RawScore=FALSE, JRFExp=2){
  stopifnot("inputs must both be dendrograms!"=is(dend1, 'dendrogram') && is(dend2, 'dendrogram'))
  if (is.integer(JRFExp)) JRFExp <- as.numeric(JRFExp) 
  stopifnot("ExpVal must be numeric or integer"=is.numeric(JRFExp))
  
  if (!is(labels(dend1), 'character')){
    dend1 <- rapply(dend1, \(x){
      if (!is.null(attr(x, 'leaf')))
        attr(x, 'label') <- as.character(attr(x, 'label'))
      return(x)
    }, how='replace')
  }
  if (!is(labels(dend2), 'character')){
    dend2 <- rapply(dend2, \(x){
      if (!is.null(attr(x, 'leaf')))
        attr(x, 'label') <- as.character(attr(x, 'label'))
      return(x)
    }, how='replace')
  }
  incommonLabs <- intersect(labels(dend1), labels(dend2))
  if (length(incommonLabs) == 0){ 
    val = c(0, NA, NA)
  } else {
    tree1ptr <- .Call("initCDend", dend1)
    on.exit(rm(tree1ptr))
    tree2ptr <- .Call("initCDend", dend2)
    on.exit(rm(tree2ptr))
    
    if (Method == 'GRF'){
      val <- .Call("GRFInfo", tree1ptr, tree2ptr, incommonLabs, FALSE, 0)
      return(GeneralizedRF(val, RawScore))
    } else if (Method == 'JRF'){
      val <- .Call("GRFInfo", tree1ptr, tree2ptr, incommonLabs, TRUE, JRFExp)
      return(JRFDist(val, RawScore))
    } else if (Method == 'RF'){
      val <- .Call("RFDist", tree1ptr, tree2ptr, incommonLabs)
      return(RFDist(val, RawScore))
    } else if (Method == 'KF') {
      val <- .Call("KFDist", tree1ptr, tree2ptr, incommonLabs)
      return(KFDist(val))
    }
    else
      stop("Method not recognized!")
  }
  
  return(val)
}
