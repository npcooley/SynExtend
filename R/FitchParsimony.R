FitchParsimony <- function(dend, num_traits, traits_list,
                           initial_state=rep(0L,num_traits), 
                           fill_ambiguous=TRUE){
  if(is.numeric(num_traits))
    num_traits <- as.integer(num_traits)
  if(!is.integer(num_traits) || num_traits < 1L){
    stop("'num_traits' must be a positive integer")
  }
  if(!inherits(dend, 'dendrogram'))
    stop("'dend' must be an object of class 'dendrogram'")
  if(length(traits_list) != num_traits)
    stop("Length of traits_list does not match num_traits!")
  if(!CheckBifurcating(dend)){
    stop("Implementation is only supported for bifurcating trees")
  }
  if(length(initial_state) != num_traits && !is.null(initial_state)){
    stop("'initial_state' must be NULL or a vector of length `num_traits`")
  }
  # upward fitch
  dend <- dendrapply(dend, \(x){
    s <- rep(0L, num_traits)
    if(is.leaf(x)){
      a <- attr(x, 'label')
      s <- vapply(seq_len(num_traits), 
                  \(i) as.integer(a %in% traits_list[[i]]), 
                  integer(1L))
    } else {
      s1 <- attr(x[[1]], 'FitchState')
      s2 <- attr(x[[2]], 'FitchState')
      for(i in seq_len(num_traits)){
        cv <- sort(c(s1[i], s2[i]))
        s[i] <- ifelse((cv[2]==2) || (cv[1]==cv[2]), cv[1], 2L)
      }
    }
    attr(x, 'FitchState') <- s
    x
  }, how='post.order')
  
  # set root to 0,0
  if(!is.null(initial_state))
    attr(dend, 'FitchState') <- initial_state
  
  # downward fitch
  dend <- dendrapply(dend, \(x){
    if(!is.leaf(x)){
      sp <- attr(x, 'FitchState')
      for(i in seq_along(x)){
        sc <- attr(x[[i]], 'FitchState')
        for(j in seq_len(num_traits))
          if(sc[j] == 2L) sc[j] <- sp[j]
        attr(x[[i]], 'FitchState') <- sc
      }
    }
    x
  }, how='post.order') 
  
  if(fill_ambiguous){
    dend <- dendrapply(dend, \(x){
      s <- attr(x, 'FitchState')
      for(i in seq_len(num_traits)){
        if(s[i] == 2L)
          s[i] <- sample(c(0L,1L), 1L)
      }
      attr(x, 'FitchState') <- s
      x
    }, how='post.order')
  }
  
  dend
}
