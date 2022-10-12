flatdendrapply <- function(dend, NODEFUN=NULL, LEAFFUN=NODEFUN, 
                           INCLUDEROOT=TRUE, ...){
  stopifnot("flatdendrapply only works on objects of class 'dendrogram'"=
              is(dend, 'dendrogram'))
  if (!is(NODEFUN, 'function') && !is(LEAFFUN, 'function'))
    stop("At least one of NODEFUN and LEAFFUN must be a function!")
  
  val <- lapply(dend, 
                \(x){
                  if (is.null(attr(x, 'leaf'))){
                    if (!is(NODEFUN, 'function'))
                      v <- list()
                    else
                      v <- list(NODEFUN(x, ...))
                    for ( child in x ) v <- c(v, Recall(child))
                    return(v)
                  } 
                  else if (!is(LEAFFUN, 'function'))
                    return(list())
                  else 
                    return(list(LEAFFUN(x, ...)))
                }
  )
  retval <- unlist(val, recursive=FALSE)
  if (!INCLUDEROOT)
    retval[[1]] <- NULL
  
  lens <- vapply(retval, length, FUN.VALUE=0L)
  atom <- vapply(retval, is.atomic, FUN.VALUE=TRUE)
  
  if(all(lens[1] == 1L) && all(atom))
    retval <- unlist(retval)
  
  return(retval)
}
