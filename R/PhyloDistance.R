CIDist <- function(val, incommonLabs, RawScore=FALSE){
  # This method is now being called "Clustering Information Distance"
  # It's called with `Method='CI'`
  if (RawScore){
    retval <- val
    CI_DISTANCE_INTERNAL <- NULL
    data('CIDist_NullDist', package="SynExtend", envir=environment())
    if (val[1] == 0){
      retval[c(2,3)] <- c(NA, NA)
    } 
    if(length(incommonLabs) < 4){
      retval[4] <- NA
    } else if (!is.null(CI_DISTANCE_INTERNAL)) {
      maxval <- 0.5*(val[2] + val[3])
      s <- (maxval - val[1]) / maxval
      leninter <- min(length(incommonLabs)-3L,197)
      rowtocheck <- c(0, CI_DISTANCE_INTERNAL[2:10,leninter], 1)
      pvals <- c(0,1,5,10,25,50,25,10,5,1,0)/100
      # move score to right side
      pvalind <- ifelse(s > 0.5, 1-s, s)
      ploc <- which(pvalind < rowtocheck)[1]
      if(ploc == 1){
        pv <- 0
      } else {
        pv <- (s - rowtocheck[ploc-1]) / (rowtocheck[ploc] - rowtocheck[ploc-1])
        pv <- pv * abs(pvals[ploc] - pvals[ploc-1]) + min(pvals[ploc-1], pvals[ploc])
        pv <- pv * 2
      }
      retval[4] <- pv
    } else {
      retval[4] <- NA
    }
    names(retval) <- c("Similarity", "dend1.Entropy", "dend2.Entropy", "p.value")
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
    if (val[1] == 0) retval[c(1,2)] <- c(NA, NA)
    return(retval)
  }
  if (val[1] == 0) return(1)
  maxval <- (val[2] + val[3])
  retval <- 1 - (val[1] / maxval)
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

PhyloDistance <- function(dend1, dend2, Method=c("CI", "RF", "KF", "JRF"), RawScore=FALSE, JRFExp=2){
  Method <- match.arg(Method)
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
    val <- c(0, NA, NA)
  } else {
    tree1ptr <- .Call("initCDend", dend1, PACKAGE="SynExtend")
    on.exit(rm(tree1ptr))
    tree2ptr <- .Call("initCDend", dend2, PACKAGE="SynExtend")
    on.exit(rm(tree2ptr))
    
    if (Method == 'CI'){
      val <- .Call("GRFInfo", tree1ptr, tree2ptr, 
                   incommonLabs, FALSE, 0, PACKAGE="SynExtend")
      return(CIDist(val, incommonLabs, RawScore))
    } else if (Method == 'JRF'){
      val <- .Call("GRFInfo", tree1ptr, tree2ptr, 
                   incommonLabs, TRUE, JRFExp, PACKAGE="SynExtend")
      return(JRFDist(val, RawScore))
    } else if (Method == 'RF'){
      val <- .Call("RFDist", tree1ptr, tree2ptr, 
                   incommonLabs, PACKAGE="SynExtend")
      return(RFDist(val, RawScore))
    } else if (Method == 'KF') {
      val <- .Call("KFDist", tree1ptr, tree2ptr, 
                   incommonLabs, PACKAGE="SynExtend")
      return(KFDist(val))
    }
    else
      stop("Method not recognized!")
  }
  
  return(val)
}
