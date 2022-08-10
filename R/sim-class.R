sim <- function(x=NA_real_, nelem, NAMES=NULL, DIAG=FALSE){
  if (missing(nelem)){
    return(as.sim(x, NAMES=NAMES, DIAG=DIAG))
  }

  if (DIAG)
    num_vals <- nelem*(nelem+1) / 2
  else
    num_vals <- nelem*(nelem-1) / 2
  if (num_vals %% length(x) != 0)
    warning("number of elements provided (", length(x), 
            ") is not a multiple of requested length (", num_vals, ")")
  v <- rep(x, length.out=num_vals)
  return(as.sim(v, NAMES=NAMES, DIAG=DIAG))
}

as.sim <- function(x, ...) UseMethod('as.sim')

as.sim.default <- function(x, ...){
  if (is(x, 'matrix')){
    return(as.sim.matrix(x, ...))
  } else if (is(x, 'vector')){
    return(as.sim.vector(x, ...))
  } else {
    stop("Input must be of type 'matrix' or type 'vector'.")
  }
}

as.sim.vector <- function(x, NAMES=NULL, DIAG=TRUE){
  # We can just solve quadratic formula if input a vector
  # n^2 - n - 2*length(x) = 0
  # a=1, b=-1, c=-2*length(x)
  b <- ifelse(DIAG, -1, 1)
  l <- 2*length(x)
  deter <- 1 + 4*l
  # Has to be the positive case to be valid
  val <- (sqrt(deter) + b) / 2
  stopifnot('x is not a valid length!'=val%%1==0)
  stopifnot('incorrect NAMES length!'=(is.null(NAMES) || length(NAMES)==val))
  
  if ( !DIAG ){
    total_num <- (val) * (val+1) / 2
    outvec <- numeric(total_num)
    svals <- (c(0, cumsum(val:1)) + 1)[-(val+1)]
    outvec[-svals] <- x
    outvec[svals] <- 1
    x <- outvec
  }
  
  if (is.null(NAMES)){
    NAMES <- as.character(seq_len(val))
  }
  
  structure(x, 
            nrow=val,
            NAMES=NAMES,
            class='sim')
}

as.sim.matrix <- function(x, ...){
  nr <- nrow(x)
  nc <- ncol(x)
  stopifnot("Matrix must be square"=nr==nc)
  if (all(x[upper.tri(x)] != t(x)[upper.tri(x)])){
    warning("Matrix is not symmetric, using upper triangle.")
  }
  
  if (is.null(rownames(x)) || is.null(colnames(x))){
    NAMES <- as.character(1:nr)
  } else if(!setequal(rownames(x),colnames(x))){
    warning("Matrix has different row and column names, using default index.")
    NAMES <- as.character(1:nr)
  } else {
    NAMES <- rownames(x)
  }
  v <- t(x)[lower.tri(x, diag=TRUE)]
  return(as.sim.vector(v, NAMES=NAMES, DIAG=TRUE))
}

show.sim <- function(s, n=10, ...){
  nr <- attr(s, 'nrow')
  CUTOFF <- FALSE
  if (n < 3) n <- 3
  if (nr > n){
    CUTOFF <- TRUE
    
  }
  ns <- attr(s, 'NAMES')
  s <- unclass(s)
  format_vals <- function(x) paste0(format(x, justify='right', width=6L), collapse='')
  nstr <- ns
  if (CUTOFF){
    nstr <- c(ns[1:(n-1)], '', ns[nr]) 
  }
  outstr <- c(format_vals(c(' ', nstr)), '\n')
  ctr <- 1L
  for ( i in seq_len(nr) ){
    if (!CUTOFF || i < (n-1) || i == nr){
      linestr <- c(ns[i], rep(' ', i-1L))
      for ( j in seq_len(nr-i+1L) ){
        linestr <- c(linestr, sprintf('%.2f', s[ctr]))
        ctr <- ctr + 1L
      }
      
      if (CUTOFF){
        linestr <- c(linestr[1:n], '…', linestr[nr+1]) 
      }
      if (CUTOFF && i == nr){
        linestr[n+1] <- ''
        linestr[n+2] <- sprintf('%.2f', s[length(s)])
      }
    } else if (CUTOFF && i == (n)){
      linestr <- c('⋮', rep('⋮', n-1), '', '⋮')
    } else {
      next
    }
    outstr <- c(outstr, format_vals(linestr), '\n')
  }
  statusstr <- paste0('\n     Similarity matrix with ', nr, ' members.\n')
  outstr <- c(outstr, statusstr)
  # invisible to quietly return NULL
  invisible(cat(outstr))
}

print.sim <- function(x, ...) show.sim(x, ...)

as.matrix.sim <- function(s){
  nr <- attr(s, 'nrow')
  ns <- attr(s, 'NAMES')
  outmat <- diag(1, nrow=nr)
  colnames(outmat) <- rownames(outmat) <- ns
  
  # have to do it this way due to column-wise ordering
  outmat[lower.tri(outmat, diag=TRUE)] <- s
  outmat <- t(outmat)
  outmat[lower.tri(outmat, diag=TRUE)] <- s
  return(t(outmat))
}

`[.sim` <- function(x, i, j){
  nr <- attr(x, 'nrow')
  ns <- attr(x, 'NAMES')
  svals <- c(0, cumsum(nr:1)) + 1
  class(x) <- 'vector'
  hasi <- !missing(i)
  hasj <- !missing(j)
  COLOUT <- FALSE
  
  if (hasi){
    if (is(i, 'character')){
      i <- sapply(i, \(ii) which(ii==ns)[1])
    }
    if (any(is.na(i))){
      stop('Incorrect indices provided')
    }
    if (any(i > nr | i < 1)){
      stop('Indices out of bounds')
    }
  }
  
  if (hasj){
    if (is(j, 'character')){
      j <- sapply(j, \(ji) which(j==ns)[1])
    }
    if (any(is.na(j))){
      stop('Incorrect indices provided')
    }
    if (any(j > nr | j < 1)){
      stop('Indices out of bounds')
    }
  }
  
  if (hasj && !hasi){
    i <- j
    hasj <- FALSE
    hasi <- TRUE
    COLOUT <- TRUE
  }
  if (!hasj){
    outmat <- matrix(NA_real_, nrow=length(i), ncol=nr)
    for ( idxi in seq_along(i) ){
      idx <- i[idxi]
      outvec <- numeric(nr)
      outvec[idx:nr] <- x[svals[idx]:(svals[idx+1]-1)]
      if (idx != 1)
        outvec[1:(idx-1)] <- x[svals[1:(idx-1)] + (idx-1):1]
      outmat[idxi,] <- outvec
    }
    rownames(outmat) <- ns[i]
    colnames(outmat) <- ns
    if (COLOUT) outmat <- t(outmat)
    return(outmat)
  }
  all_accessed <- expand.grid(i, j)
  idxvec <- integer(nrow(all_accessed))
  namevec <- character(nrow(all_accessed))
  for (idx in seq_len(nrow(all_accessed))){
    i1 <- all_accessed[idx, 1]
    i2 <- all_accessed[idx, 2]
    if ( i1 > i2 ){
      namevec[idx] <- paste(ns[i2], ns[i1], sep=',')
      idxvec[idx] <- svals[i2] + (i1 - i2)
    } else if ( i1 == i2 ){
      namevec[idx] <- paste(ns[i1], ns[i1], sep=',')
      idxvec[idx] <- svals[i1] 
    } else {
      namevec[idx] <- paste(ns[i1], ns[i2], sep=',')
      idxvec[idx] <- svals[i1] + (i2 - i1)
    }
  }
  outvec <- x[idxvec]
  names(outvec) <- namevec
  return(outvec)
}

`[<-.sim` <- function(x, i, j, value){
  nr <- attr(x, 'nrow')
  ns <- attr(x, 'NAMES')
  svals <- c(0, cumsum(nr:1)) + 1
  class(x) <- c('vector')
  hasi <- !missing(i)
  hasj <- !missing(j)
  
  if (hasi){
    if (is(i, 'character')){
      i <- sapply(i, \(ii) which(ii==ns)[1])
    }
    if (any(is.na(i))){
      stop('Incorrect indices provided')
    }
    if (any(i > nr | i < 1)){
      stop('Indices out of bounds')
    }
  }
  
  if (hasj){
    if (is(j, 'character')){
      j <- sapply(j, \(ji) which(j==ns)[1])
    }
    if (any(is.na(j))){
      stop('Incorrect indices provided')
    }
    if (any(j > nr | j < 1)){
      stop('Indices out of bounds')
    }
  }
  
  # Symmetric so we can do this to simplify
  if (hasj && !hasi){
    i <- j
    hasi <- TRUE
    hasj <- FALSE
  }
  
  ## Checking for valid inputs
  # 3 cases: none given, only i given, i and j both given
  num_given <- length(value)
  if (!hasi && !hasj){
    if (length(x) %% num_given != 0){
      warning("number of items to replace is not a multiple of replacement length")
    } 
    value <- rep(value, length.out=length(x))
    return(as.sim(value, NAMES=ns))
  }
  li <- length(i)
  if (hasi && !hasj){
    num_req <- nr * li
  } else {
    lj <- length(j)
    num_req <- li * lj
  }
  
  if (num_req %% num_given != 0){
    warning("number of items to replace is not a multiple of replacement length")
  } 
  value <- rep(value, length.out=num_req)
  
  if (!hasj){
    for ( idxi in seq_along(i) ){
      idx <- i[idxi]
      offset <- nr * (idxi-1)
      x[svals[idx]:(svals[idx+1]-1)] <- value[(idx:nr) + offset]
      if (idx != 1)
        x[svals[1:(idx-1)] + (idx-1):1] <- value[(1:(idx-1)) + offset]
    }
    class(x) <- 'sim'
    return(x)
  }
  
  all_accessed <- expand.grid(i, j)
  idxvec <- integer(nrow(all_accessed))
  for (idx in seq_len(nrow(all_accessed))){
    i1 <- all_accessed[idx, 1]
    i2 <- all_accessed[idx, 2]
    if ( i1 > i2 ){
      idxvec[idx] <- svals[i2] + (i1 - i2)
    } else if ( i1 == i2 ){
      idxvec[idx] <- svals[i1] 
    } else {
      idxvec[idx] <- svals[i1] + (i2 - i1)
    }
  }
  
  x[idxvec] <- value
  class(x) <- 'sim'
  return(x)
  #return(as.sim(v, NAMES=ns, DIAG=TRUE))
  
}

diag <- function(x, ...) UseMethod('diag')

diag.default <- base::diag

diag.sim <- function(x, nrow, ncol, names=TRUE){
  nr <- attr(x, 'nrow')
  svals <- c(0, cumsum(nr:2)) + 1
  v <- unclass(x)
  return(v[svals])
}

`diag<-` <- function(x, value, ...) UseMethod('diag<-')

`diag<-.default` <- base::`diag<-`

`diag<-.sim` <- function(x, value, ...){
  nr <- attr(x, 'nrow')
  if (nr %% length(value) != 0){
    warning("number of items to replace is not a multiple of replacement length")
  }
  value <- rep(value, length.out=nr)
  svals <- c(0, cumsum(nr:2)) + 1
  class(x) <- 'vector'
  x[svals] <- value
  class(x) <- 'sim'
  return(x)
}

names.sim <- function(x) {
  return(attr(x, 'NAMES'))
}

`names<-.sim` <- function(x, value){
  ns <- attr(x, 'NAMES')
  stopifnot('Incorrect number of names provided.'=length(value) == length(ns))
  attr(x, 'NAMES') <- value
  return(x)
}

as.data.frame.sim <- function(x) {
  l <- length(x)
  nr <- attr(x, 'nrow')
  i1 <- i2 <- rep(NA_integer_, l)
  v <- rep(NA_real_, l)
  ctr <- 1
  for ( i in seq_len(nr) ){
    for ( j in i:nr ) {
      i1[ctr] <- i
      i2[ctr] <- j
      v[ctr] <- x[i,j]
      ctr <- ctr + 1
    }
  }
  
  return(data.frame(row=i1, col=i2, value=v))
}