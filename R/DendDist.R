DendDist <- function(dend1, dend2, Method, RawScore=FALSE){
  stopifnot("inputs must both be dendrograms!"=is(dend1, 'dendrogram') && is(dend2, 'dendrogram'))
  
  if (!is(labels(dend1), 'character')){
    dend1 <- dendrapply(dend1, \(x){
      if (!is.null(attr(x, 'leaf')))
        attr(x, 'label') <- as.character(attr(x, 'label'))
      return(x)
    })
  }
  if (!is(labels(dend2), 'character')){
    dend2 <- dendrapply(dend2, \(x){
      if (!is.null(attr(x, 'leaf')))
        attr(x, 'label') <- as.character(attr(x, 'label'))
      return(x)
    })
  }
  incommonLabs <- intersect(labels(dend1), labels(dend2))
  if (length(incommonLabs) == 0){ 
    val = c(0, NA, NA)
  } else {
    tree1ptr <- .Call("initCDend", dend1)
    on.exit(rm(tree1ptr))
    tree2ptr <- .Call("initCDend", dend2)
    on.exit(rm(tree2ptr))
    
    if (Method == 'GRF')
      val <- .Call("GRFInfo", tree1ptr, tree2ptr, incommonLabs)
    else if (Method == 'RF')
      val <- .Call("RFDist", tree1ptr, tree2ptr, incommonLabs)
    else
      stop("Method not recognized!")
  }
  
  return(val)
}

GeneralizedRF <- function(dend1, dend2, RawScore=FALSE){
  val <- DendDist(dend1, dend2, Method="GRF")
  if (RawScore){
    retval <- val
    names(retval) <- c("Similarity", "dend1.Entropy", "dend2.Entropy")
    if (val[1] == 0) retval[1:2] <- c(NA, NA)
    return(retval)
  }
  
  if (val[1] == 0) return(1)
  maxval <- (val[2] + val[3])
  retval <- (maxval - val[1]) / maxval
  if (maxval == 0)
    retval <- as.integer(val[1] != 0)
  return(retval)
}

RFDist <- function(dend1, dend2, RawScore=FALSE){
  val <- DendDist(dend1, dend2, Method="RF")
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
