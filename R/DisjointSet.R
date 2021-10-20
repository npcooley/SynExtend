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
    F1 <- factor(p1,
                 levels = AllIds)
    F2 <- factor(p2,
                 levels = AllIds)
    
    FI1 <- as.integer(F1)
    FI2 <- as.integer(F2)
    FC1 <- as.character(F1)
    FC2 <- as.character(F2)
    
    FInts <- c(FI1, FI2)
    FChars <- c(FC1, FC2)
    Key1 <- !duplicated(FInts)
    IntMap <- cbind("FactorRep" = FInts[Key1],
                    "UniqueIDs" = FChars[Key1])
    
    # run FindSets
    IntRes <- FindSets(p1 = FI1,
                       p2 = FI2,
                       Verbose = TRUE)
    
    rm(list = c("FI1",
                "FI2",
                "FC1",
                "FC2",
                "FInts",
                "FChars",
                "Key1"))
    
    if (Verbose) {
      cat("\nAssigning single linkage clusters.\n")
    }
    
    # all nodes
    Nodes <- IntRes[, 1L, drop = TRUE]
    # all representatives
    Reps <- IntRes[, 2L, drop = TRUE]
    # set the order to the unique reps
    AllIds <- IntMap[match(x = IntMap[, 1L],
                           table = Nodes), 2L]
    
    UClusts <- split(x = AllIds,
                     f = Reps)
    
    if (Verbose) {
      cat("Assignments complete.\n")
    }
    if (Verbose) {
      cat("\n")
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    
    names(UClusts) <- seq(length(UClusts))
    
  } else {
    stop ("No method exists for object.")
  }
  
  # return a list of character identifiers for disjoint sets
  
  return(UClusts)
  
}