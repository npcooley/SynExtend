ShuffleC <- function(vec, n=length(vec), replace=FALSE){
  stopifnot(n > 0 && (replace || n <= length(vec)))
  if (replace)
    return(ShuffleWithRecomb(vec, n))
  
  trimlen <- seq_len(n)
  if(is(vec, 'integer')){
    return(.C("shuffleRInt", vec, n)[[1]][trimlen])
  } else {
    return(vec[.C("shuffleRInt", seq_len(length(vec)), n)[[1]]][trimlen])
  }
}

ShuffleWithRecomb <- function(vec, n=length(vec)){
  n <- as.integer(n)
  l <- length(vec)
  idxes <- .C("shuffleRRepl", seq_len(n), n)[[1L]] %% l + 1L
  return(vec[idxes])
}
