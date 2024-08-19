subset.dendrogram <- function(x, subset, invert=FALSE, ...){
  if(!inherits(x, 'dendrogram'))
    stop("'x' is not a dendrogram!")
  if(invert){
    treelabs <- labels(x)
    subset <- treelabs[!(treelabs%in%subset)]
  }
  if(length(subset) == 0L){
    warning("Subsetting resulted in no leaves")
    tmp <- list()
    class(tmp) <- class(x)
    attr(tmp, 'members') <- 0L
    attr(tmp, 'height') <- 0.0
    return(tmp)
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
      nmemb <- vapply(y, attr, which='members', integer(1L))
      attr(y, 'members') <- sum(nmemb)

      l <- length(y)
      if(is.leaf(y[[1]]) && is.leaf(y[[l]])){
        mp <- (sum(nmemb) - 1) / 2
      } else if(is.leaf(y[[1]])){
        mp <- (sum(nmemb[-l]) + attr(y[[l]], 'midpoint')) / 2
      } else if(is.leaf(y[[l]])){
        mp <- (attr(y[[1]], 'midpoint') + sum(nmemb[-1])) / 2 + attr(y[[1]], 'midpoint')
      } else {
        mp <- (sum(nmemb[-l]) + attr(y[[1]], 'midpoint') + attr(y[[l]], 'midpoint')) / 2
      }
      attr(y, 'midpoint') <- mp

      return(y)
    }
  }, how='post.order')
}
