##### -- FastDend Class, internally used for dendrogram manipulation  ---------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

setClass('FastDend',
          slots = list(branches = 'matrix', leaflabels = 'character', 
                       nleaves = 'numeric', nbranches = 'numeric',
                       root = 'numeric', dendrogram = 'ANY'))

setMethod(show, 'FastDend', function(object) {
  cat("'FastDend' object with", object@nleaves, 'leaves and', 
      object@nbranches, 'branches.\n')
})

setMethod('plot', 'FastDend', function(x) plot(x@dendrogram))


setGeneric('branches', function(x) standardGeneric('branches'))
setMethod('branches', 'FastDend', function(x) {
  return(x@branches)
})

setGeneric('leaves', function(x) standardGeneric('leaves'))
setMethod('leaves', 'FastDend', function(x) {
  return(x@leaflabels)
})

setGeneric('root', function(x) standardGeneric('root'))
setMethod('root', 'FastDend', function(x) {
  return(x@root)
})

setMethod('labels', 'FastDend', function(object) leaves(object))

setGeneric('subtree', function(x, y) standardGeneric('subtree'))
setMethod('subtree', 'FastDend', function(x, y) {
  tree <- x@dendrogram
  path <- rootPath(x, y)
  path <- path[-length(path)]
  subtree <- tree
  for (element in path) subtree <- subtree[[element]]
  return(FastDend(subtree))
})

setGeneric('parentTree', function(fDend, leafNode, height) standardGeneric('parentTree'))
setMethod('parentTree', 'FastDend', function(fDend, leafNode, height) {
  stopifnot("leafNode must be type 'character' or 'numeric'." = (is(leafNode, 'character') || is(leafNode, 'numeric')))
  stopifnot("'height' must be a non-negative integer." = (is(height, 'numeric') && height >= 0 && as.integer(height) == height))
  if (is(leafNode, 'character')){
    leafNode <- which(fDend@leaflabels == leafNode)
  }
  stopifnot('Leaf does not exist in tree.' = (length(leafNode) == 1 && leafNode %in% seq_len(fDend@nleaves)))
  
  b <- fDend@branches
  root <- fDend@root
  curnode <- leafNode
  while(height > 0) {
    if(curnode == root){
      warning('Height greater than maximal possible subtree.')
      return(fDend@dendrogram)
    }
    curnode <- b[curnode, 'parent']
    height <- height - 1
  }
  
  return(subtree(fDend, curnode))
})

setGeneric('rootPath', function(x, y) standardGeneric('rootPath'))
setMethod('rootPath', 'FastDend', function(x, y) {
  tree <- x@dendrogram
  stopifnot(is(tree, 'dendrogram'))
  stopifnot("'y' must be type 'character' or 'numeric'." = (is(y, 'character') || is(y, 'numeric')))
  if (is(y, 'character')) {
    y <- which(x@leaflabels == y)
  }
  b <- x@branches
  stopifnot('Branch does not exist in tree.' = (length(y) == 1 && y %in% seq_len(nrow(b))))
  path <- pathnames <- rep(NA, x@nbranches+1)
  root <- x@root
  curpos <- y
  i <- x@nbranches
  path[i] <- 0
  pathnames[i] <- curpos
  i <- i - 1
  
  while(curpos != root){
    parent <- b[curpos, 'parent']
    path[i] <- ifelse(b[parent, 'left'] == curpos, 1, 2)
    pathnames[i] <- parent
    curpos <- parent
    i <- i - 1
  }
  path <- path[!is.na(path)]
  pathnames <- pathnames[!is.na(pathnames)]
  names(path) <- pathnames
  return(path)
})

FastDend <- function(dend){
  dfs_helper <- function(node, mat, labels, curval){
    left <- node[[1]]
    right <- node[[2]]
    leftlength <- attr(node, 'height') - attr(left, 'height')
    rightlength <- attr(node, 'height') - attr(right, 'height')
    
    if (is.null(attr(left, 'leaf'))){
      nextval <- which.min(mat[,'label'])
      mat[curval, 'left'] <- nextval
      mat[nextval, 'parent'] <- curval
      mat[nextval, 'label'] <- nextval
      mat[nextval, 'length'] <- leftlength
      mat <- dfs_helper(left, mat, labels, nextval)
    } else {
      leaflabel <- which(labels == attr(left, 'label'))
      mat[curval, 'left'] <- leaflabel
      mat[leaflabel, 'parent'] <- curval
      mat[leaflabel, 'length'] <- leftlength
    }
    
    if (is.null(attr(right, 'leaf'))){
      nextval <- which.min(mat[,'label'])
      mat[curval, 'right'] <- nextval
      mat[nextval, 'parent'] <- curval
      mat[nextval, 'label'] <- nextval
      mat[nextval, 'length'] <- rightlength
      mat <- dfs_helper(right, mat, labels, nextval)
    } else {
      leaflabel <- which(labels == attr(right, 'label'))
      mat[curval, 'right'] <- leaflabel
      mat[leaflabel, 'parent'] <- curval
      mat[leaflabel, 'length'] <- rightlength
    }
    
    return(mat)
  }
  
  
  labels <- labels(dend)
  nleaves <- length(labels)
  leaflabels <- seq_len(nleaves)
  names(labels) <- leaflabels
  nbranches <- 2*nleaves + 3
  mat <- matrix(c(-1, 0, 0, 0, 0, 0, 1), nrow=nbranches, ncol = 7, byrow = TRUE)
  colnames(mat) <- c('label', 'isLeaf', 'left', 'right', 'parent', 'length', 'support')
  mat[leaflabels, 'isLeaf'] <- 1
  mat[leaflabels, 'label'] <- leaflabels
  rootval <- nleaves + 1
  mat[rootval, 'label'] <- 0
  mat[rootval, 'parent'] <- -1
  mat[rootval, 'length'] <- attr(dend, 'height')
  mat <- dfs_helper(dend, mat, labels, rootval)
  mat <- mat[mat[,'label'] >= 0,]
  
  d <- new('FastDend', branches = mat, leaflabels = labels,
           nleaves = nleaves, nbranches = nrow(mat)-1,
           root = rootval, dendrogram = dend)
  return(d)
}

# Extra Methods
fitch_up <- function(fDend, ov){
  ov <- c(ov, rep(-1, nrow(branches(fDend)) - length(ov)))
  helperfunc <- function(branchmat, node){
    left <- branchmat[node, 'left']
    if (ov[left] < 0){
      helperfunc(branchmat, left)
    }
    right <- branchmat[node, 'right']
    if (ov[right] < 0){
      helperfunc(branchmat, right)
    }
    
    lv <- ov[left]
    rv <- ov[right]
    fillval <- 2
    if (lv != rv && (lv == 2 || rv == 2))
      fillval <- min(lv, rv)
    else
      fillval <- ifelse(xor(lv, rv), 2, lv)
    ov[node] <<- fillval
  }
  helperfunc(branches(fDend), root(fDend))
  return(ov)
}

fitch_down <- function(fDend, ov, pushUp=TRUE, rootFill = 0){
  helperfunc <- function(branchesmat, node){
    if (branchesmat[node, 'isLeaf'] != 1) {
      if (ov[node] == 2){
        if (branchesmat[node, 'parent'] == -1){
          ov[node] <<- rootFill
        } else {
          pval <- ov[branchesmat[node, 'parent']]
          ov[node] <<- ifelse(pushUp, abs(1-pval), pval)
        }
      } 
      helperfunc(branchesmat, branchesmat[node, 'left'])
      helperfunc(branchesmat, branchesmat[node, 'right'])
    }
  }
  
  helperfunc(branches(fDend), root(fDend))
  return(ov)
}

fitch_parsimony<- function(fDend, occurrenceVec, moveEventsUpward=TRUE){
  ovf <- fitch_up(fDend, occurrenceVec)
  ovf <- fitch_down(fDend, ovf, moveEventsUpward)
  return(ovf)
}

gain_loss_vec <- function(fDend, occurrenceVec){
  glvec <- rep(0, length(occurrenceVec))
  helperfunc <- function(bmat, branch){
    parent <- bmat[branch, 'parent']
    if (parent != -1){ # not the root
      glvec[branch] <<- occurrenceVec[branch] - occurrenceVec[parent]
    }
    
    if (bmat[branch, 'isLeaf'] != 1){
      helperfunc(bmat, bmat[branch, 'left'])
      helperfunc(bmat, bmat[branch, 'right'])
    }
  }
  
  helperfunc(branches(fDend), root(fDend))
  return(glvec)
}

calc_SId_mat <- function(fd, IdOnly=FALSE){
  fdbranch <- fd@branches
  branchmat <- matrix(0, ncol=nrow(fdbranch), nrow=nrow(fdbranch))
  blengths <- numeric(length=nrow(fdbranch))
  for (i in seq_len(nrow(fdbranch))){
    r <- fdbranch[i,]
    v <- ifelse(r[1] == 0, 302, r[1])
    blengths[v] <- r[6]
    if (r[2] == 1)
      next
    if(r[3] != 0)
      branchmat[v, r[3]] <- 1
    if(r[4] != 0)
      branchmat[v, r[4]] <- 1
  }
  blengths <- blengths / sum(blengths)
  Cmat <- matrix(0, nrow=length(blengths), ncol=length(blengths))
  for (i in seq_along(blengths)){
    for (j in seq_along(blengths)){
      if (i == j)
        Cmat[i,j] <- blengths[i] * (1-blengths[j])
      else
        Cmat[i,j] <- -1*blengths[i]*blengths[j]
    }
  }
  
  S1 <- branchmat
  Id <- diag(nrow=nrow(S1), ncol=ncol(S1))
  S <- S_c <- S1
  if (!IdOnly){
    repeat{
      S_c <- S_c %*% S1
      if (all(S_c == 0))
        break
      S <- S + S_c
    }
  }
  
  return(list(S=S, Id=Id, blengths=blengths, Cmat=Cmat))
}

generateGainLossVec <- function(fDend, PAvec, moveEventsUpward=TRUE){
  stopifnot(is(fDend, 'FastDend'))
  stopifnot(length(leaves(fDend)) == length(PAvec))
  occvec <- fitch_parsimony(fDend, PAvec, moveEventsUpward)
  glv <- gain_loss_vec(fDend, occvec)
  return(glv)
}