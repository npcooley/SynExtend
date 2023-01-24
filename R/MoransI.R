##### -- Moran's I to estimate Phylogenetic Autocorrelation -------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

######
# T     := the tree we are considering
# D     := some measure of distance between individuals in the tree
# Mij   := M[i][j]
# Dij   := distance between entry i and j
# n     := number of individuals in the tree
# sum_i := sum over i from 1:n
#
# Moran's I is defined as the following:
# (n / sum(D)) * (sum_i{ sum_j { Dij(xi - mean(x))(xj - mean(x)) } }) / sum_i{ (xi-mean(x))**2 }
#
# Note that there is some redundancy in this measure since
# we add both i,j and j,i in all cases.
# Note also that Dij := 0 when i=j.
#####

MoransI <- function(values, weights, alternative='two.sided'){
  hyptest <- pmatch(alternative, c('two.sided', 'less', 'greater'))
  if(is.na(alternative) || length(alternative) > 1){
    stop("'alternative' must be an unambiguous abbreviation of one of the
         following: 'two.sided', 'less', 'greater'")
  }
  if(!is(weights, 'dist')){
    if (is(weights, 'matrix')){
      warning("'weights' is matrix rather than 'dist', attempting to convert...")
      weights <- as.dist(weights)
    } else {
      stop("'weights' must be of class 'dist'!")
    }
  }
  if(length(values) != attr(weights, "Size")){
    stop("'values' and 'weights' do not imply the same number of individuals!")
  }
  if(!is(values, 'atomic')){
    stop("'vals' must be a numeric vector")
  }
  
  if (length(values) <= 1){
    warning("Correlation with a single value is meaningless.")
    return(list(observed=values[1], expected=NA_real_, sd=Inf, p.value=1))
  }

  res <- .Call('MoransI', as.double(values), as.double(c(weights)), length(values))
  if(is.null(res[3])){
    # This really only happens when all the values are zero,
    # or like really really close to zero (on the order of <1e-295)
    return(list(observed=0, expected=-1/(length(values)-1), sd=0, p.value=1))
  }
  retval <- list(observed=res[1], expected=res[2], sd=sqrt(res[3]))
  if (length(values) <= 3){
    warning('Fewer than 3 values, variance is infinite!')
    retval$sd <- Inf
    retval$p.value <- 1 
    return(retval)
  }  
  denom <- ifelse(retval$sd==0, 1, retval$sd)
  p <- NULL
  if (hyptest == 1){
    p <- 2*pnorm(abs(retval$observed - retval$expected) / denom, lower.tail=FALSE)
  } else if (hyptest==2){
    p <- pnorm((retval$observed - retval$expected) / denom, lower.tail=TRUE)
  } else {
    p <- pnorm((retval$observed - retval$expected) / denom, lower.tail=FALSE)
  }
  retval$p.value <- p 
  return(retval)
}
