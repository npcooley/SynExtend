HungarianAssignment <- function(m){
  on.exit(.C("hungarianCleanup", PACKAGE="SynExtend"))
  if (min(m) < 0){
    m <- m - min(m)
  }
  res <- .Call("HungarianAssignment", c(t(m)), dim(m))
  names(res) <- seq_along(res)
  res <- res + 1
  return(res)
}