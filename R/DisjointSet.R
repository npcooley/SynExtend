# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

DisjointSet <- function(Pairs,
                        Verbose = FALSE) {
  
  # grab p1 and p2 columns and perform all the necessary overhead for the 
  # FindSets() function
  # expect to take in other parameters in the future and run an implementation
  # of MCL
  
  if (Verbose) {
    pBar <- txtProgressBar(style = 1L)
    TimeStart <- Sys.time()
  }
  
  if (is(object = Pairs,
         class2 = "PairSummaries")) {
    
    if (nrow(Pairs) < 1L) {
      stop ("No rows present?")
    }
    
    # collect IDs for factorization
    p1 <- Pairs$p1
    p2 <- Pairs$p2
    P1Ids <- unique(p1)
    P2Ids <- unique(p2)
    AllIds <- unique(c(P1Ids,
                       P2Ids))
    F1 <- as.integer(factor(p1,
                            levels = AllIds))
    F2 <- as.integer(factor(p2,
                            levels = AllIds))
    
    # run FindSets
    IntRes <- FindSets(p1 = F1,
                       p2 = F2,
                       Verbose = Verbose)
    
    UClusts <- unique(IntRes[, 2L])
    L <- length(UClusts)
    FinRes <- vector(mode = "list",
                     length = L)
    
    cat("\nAssigning single linkage clusters.\n")
    for (m1 in seq_along(UClusts)) {
      # IntRes[, 1L] is the integer representations of character ids
      # present in p1
      # IntRes[, 2L] is the integer representation of the "parent" node for
      # a given cluster
      # for each parent node select all associated nodes
      # then extract the character ids of those nodes
      FinRes[[m1]] <- IntRes[IntRes[, 2L, drop = TRUE] == UClusts[m1], 1L, drop = TRUE]
      # grab IDs from p1
      # grap IDs from p2
      # unique them
      FinRes[[m1]] <- unique(c(p1[F1 %in% FinRes[[m1]]],
                               p2[F1 %in% FinRes[[m1]]]))
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / L)
      }
    }
    if (Verbose) {
      cat("\n")
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    
    names(FinRes) <- seq(length(FinRes))
    
  } else {
    stop ("No method exists for object.")
  }
  
  # return a list of character identifiers for disjoint sets
  
  return(FinRes)
  
}