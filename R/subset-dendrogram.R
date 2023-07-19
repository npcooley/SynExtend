subset.dendrogram <- function(x, subset, invert=FALSE, ...){
  if(!inherits(x, 'dendrogram'))
    stop("'x' is not a dendrogram!")
  if(invert){
    treelabs <- labels(x)
    subset <- treelabs[!(treelabs%in%subset)]
  }
  if(length(subset) == 0L){
    warning("Subsetting resulted in no leaves")
    x <- list()
    class(x) <- class(x)
    attr(x, 'members') <- 0L
    attr(x, 'height') <- 0.0
    return(x)
  }
  
  dendrapply(x, \(y){
    ## Leaves
    if(is.leaf(y)){
      if(attr(y, 'label') %in% subset) 
        return(y)
      else 
        return(NULL)
    }
    
    ## Internal Nodes
    nonNull <- which(!vapply(y, is.null, logical(1L)))
    if(length(nonNull) == 0L){
      return(NULL)
    } else if(length(nonNull) == 1L){
      return(y[[nonNull]]) 
    } else {
      nmemb <- sum(vapply(y, attr, which='members', integer(1L)))
      attr(y, 'members') <- nmemb
      attr(y, 'midpoint') <- (nmemb-1) / 2
      return(y)
    }
  }, how='post.order')
}
