# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

FindSets <- function(p1,
                     p2,
                     Verbose = FALSE) {
  
  # p1 and p2 here are the integer representatives of identifiers
  # it is likely these could be converted to numerics in cases of Very Large
  # Data, but that seems unnecessary for now
  # a constraint based on construction is that all(p2 > p1) should be true,
  # though this seems unnecessary from the point of view of the algorithm
  if (!is.integer(p1) |
      !is.integer(p2)) {
    stop ("Nodes must be represented by integers.")
  }
  if (length(p1) != length(p2)) {
    stop ("Unpaired nodes present.")
  }
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    L1 <- length(p1)
    cat("\nPass 1:\n")
  }
  Origins <- sort(unique(c(unique(p1),
                           unique(p2))))
  Nodes <- seq(length(Origins))
  L <- length(Nodes)
  Rank <- vector(mode = "integer",
                 length = L)
  
  # a visual map of algo progress:
  # Res <- matrix(data = 0L,
  #               nrow = L,
  #               ncol = (1L + length(p1) + L))
  # Res[, 1L] <- Nodes
  for (m1 in seq_along(p1)) {
    
    # find root of p1
    x <- p1[m1]
    while (x != Nodes[x]) {
      x <- Nodes[x]
    }
    
    # find root of p2
    y <- p2[m1]
    while (y != Nodes[y]) {
      y <- Nodes[y]
    }
    
    if (x == y) {
      # do nothing
    } else if (Rank[x] < Rank[y]) {
      Nodes[x] <- y
    } else if (Rank[x] > Rank[y]) {
      Nodes[y] <- x
    } else {
      Nodes[y] <- x
      Rank[x] <- Rank[x] + 1L
    }
    # visually map the progress: 
    # Res[, m1 + 1L] <- Nodes
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / L1)
    }
  }
  if (Verbose) {
    TimeEnd <- Sys.time()
    cat("\n")
    print(TimeEnd - TimeStart)
    TimeStart <- Sys.time()
    cat("\nPass 2:\n")
    L2 <- length(Nodes)
  }
  # end first pass, nodes are pointed upward only so far as paths have been explored
  # while scrolling through pairs

  # nodes whose pointer is towards another node that itself points elsewhere must be re-rooted
  # i.e. scroll through nodes and ask if parents are roots, if they are not, chase
  # the known paths till you get to a root
  for (m1 in seq_along(Nodes)) {
    if (Nodes[m1] != m1) {
      # node is not it's own root
      while (Nodes[Nodes[m1]] != Nodes[m1]) {
        Nodes[m1] <- Nodes[Nodes[m1]]
      }
    }
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / L2)
    }
    # visually map the progress
    # Res[, m1 + length(p1) + 1L] <- Nodes
  }
  if (Verbose) {
    cat("\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  
  # return(Res)
  return(cbind(Origins,
               Nodes))
}
