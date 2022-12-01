ShuffleC <- function(vec){
  l <- length(vec)
  stopifnot(l > 0)
  if(l==1) return(vec)
  
  if(is(vec, 'integer')){
    return(.C("shuffleRInt", vec, l)[[1]])
  } else {
    return(vec[.C("shuffleRInt", seq_len(l), l)[[1]]])
  }
}