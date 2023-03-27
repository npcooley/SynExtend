dendrapply <- function(X, FUN, ..., how=c("pre.order", "post.order")){
  ##
  ## dendrapply: applies function recursively to dendrogram object
  ## -------------
  ## Author: Aidan Lakshman (AHL27@pitt.edu), Date: 2/28/2023
  ## Original function by Martin Maechler, 2004
  ##
  FUN <- match.fun(FUN)
  apply_method <- match.arg(how)
  travtype <- switch(apply_method,
                     pre.order=0L,
                     post.order=1L)
  objclass <- class(X)
  ## Free allocated memory in case of early termination
  on.exit(.C("free_dendrapply_list"))
  if( !inherits(X, "dendrogram") ) stop("'X' is not a dendrogram")
  wrapper <- function(node) {
    ## VECTOR_ELT always unclasses the object
    ## This solution works as long as all nodes 
    ## have the same class, although indexing
    ## node[[i]] will use [[.list
    class(node) <- objclass
    res<-FUN(node, ...)
    if(travtype == 0L && !is.leaf(node)){
      ## catch for dendrapply(d, labels)
      if(!is.list(res)){
        res <- as.list(res)
      }
      ## catch for dendrapply(d, \(x) list())
      if(length(res) < (n <- length(node))){
        res <- vector('list', n)
      }
      
      res[seq_len(n)] <- node
    }
    res
  }
  # If we only have one node, it'll hang
  # We can get around this by just applying the function to the leaf
  # and returning--no need for C code here.
  if(is.leaf(X)){
    return(wrapper(X))
  }
  return(.Call("do_dendrapply", X, wrapper, parent.frame(), travtype))
}