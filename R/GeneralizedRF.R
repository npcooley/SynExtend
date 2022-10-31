GeneralizedRF <- function(dend1, dend2){
  stopifnot("inputs must both be dendrograms!"=is(dend1, 'dendrogram') && is(dend2, 'dendrogram'))
  
  incommonLabs <- intersect(labels(dend1), labels(dend2))
  if (length(incommonLabs) == 0) return(0)
  
  tree1ptr <- .Call("initCDend", dend1)
  on.exit(rm(tree1ptr))
  tree2ptr <- .Call("initCDend", dend2)
  on.exit(rm(tree2ptr))
  
  val <- .Call("GRFInfo", incommonLabs)
  if (val[0] == 0) return(1)
  maxval <- (val[2] + val[3]) / 2
  retval <- maxval - val[1]
  return(retval)
}