ShuffleC <- function(vec, n=length(vec)){
  stopifnot(n > 0 && n <= length(vec))
  trimlen <- seq_len(n)
  if(is(vec, 'integer')){
    return(.C("shuffleRInt", vec, n)[[1]][trimlen])
  } else {
    return(vec[.C("shuffleRInt", seq_len(length(vec)), n)[[1]]][trimlen])
  }
}