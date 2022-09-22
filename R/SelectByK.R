###### -- drop predicted pairs based on user PID confidence -------------------
# version 1:
# users can select cluster number total by moving the quantile of within
# cluster sum of squares to select with
# user confidence selects for the PID centroid of the clusters that users are comfortable with
# so you will retain pairs whose PID is below that centroid, but appear closer in distance
# to the lowest acceptable centroid than to the highest unacceptable centroid

SelectByK <- function(Pairs,
                      UserConfidence = 0.5,
                      ClusterSelect = 0.4,
                      MaxClusters = 15L,
                      ReturnAllCommunities = FALSE,
                      Verbose = FALSE,
                      ShowPlot = FALSE) {
  # require both Score and PID
  if (!is(object = Pairs,
          class2 = "PairSummaries")) {
    stop ("Pairs must be an object of class 'PairSummaries'.")
  }
  COLS <- colnames(Pairs)
  if (!all(c("SCORE", "PID") %in% COLS)) {
    stop("PairSummaries object must contain calculated PIDs and Scores.")
  }
  if (MaxClusters >= nrow(Pairs)) {
    stop ("User has requested more max clusters than rows exist in 'PairSummaries' object.")
  }
  if (MaxClusters <= 2L) {
    warning ("User has requested a number of clusters that may not provide adequate group separation.")
  }
  
  # normalize score,
  # convert absolute matches to relative
  dat1 <- data.frame("RelativeMatch" = Pairs$ExactMatch * 2L / (Pairs$p1FeatureLength + Pairs$p2FeatureLength),
                     "Consensus" = Pairs$Consensus,
                     "PID" = Pairs$PID,
                     "SCORE" = Pairs$SCORE / max(Pairs$SCORE),
                     "TetDist" = Pairs$TetDist)
  
  # cluster out to user specified total clusters
  NClust <- 2:MaxClusters
  kmc <- vector(mode = "list",
                length = length(NClust))
  for (m1 in seq_along(NClust)) {
    kmc[[m1]] <- kmeans(x = dat1,
                        centers = NClust[m1],
                        iter.max = 25L,
                        nstart = 25L)
    
    if (Verbose) {
      print(m1)
    }
  }
  
  wss <- sapply(X = kmc,
                FUN = function(x) {
                  x$tot.withinss
                })
  
  # select the cluster number based on quantile of within cluster sum of squares
  # that users input, default is 50%, seems to be fine for tests of initial version
  EvalClust <- max(which(wss >= quantile(wss, ClusterSelect))) + 1L
  if (EvalClust >= MaxClusters) {
    warning("Evaluated clusters may be insufficient for this task. Retry with a larger 'MaxClusters'.")
    EvalClust <- MaxClusters
  }
  # EvalClust <- max(which(diff(wss) <= quantile(diff(wss), ClusterSelect))) + 1L
  
  # print(wss)
  # print(quantile(wss, ClusterSelect))
  # print(EvalClust)
  
  res1 <- kmc[[EvalClust]]$centers[, "PID"]
  res1 <- names(res1)[res1 >= UserConfidence]
  
  res2 <- kmc[[EvalClust]]$cluster %in% as.integer(res1)
  
  res3 <- Pairs[res2, ]
  res4 <- vector(mode = "list",
                 length = nrow(kmc[[EvalClust]]$centers))
  for (m2 in seq_along(res4)) {
    # kmeans increments clusters by 1 from 1 by default, though numerically
    # they my not be ordered in increasing PID
    res4[[m2]] <- Pairs[kmc[[EvalClust]]$cluster %in% m2, ]
    attr(res4[[m2]], "GeneCalls") <- NULL
  }
  if (ShowPlot) {
    plot(x = 0,
         y = 0,
         type = "n",
         xlab = "Frequency",
         ylab = "PID",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = "Cluster CDFs")
    for (m2 in seq_along(res4)) {
      points(x = sort(res4[[m2]]$PID),
             y = seq_len(nrow(res4[[m2]])) / nrow(res4[[m2]]),
             pch = 1,
             col = m2)
    }
  }
  
  if (ReturnAllCommunities) {
    # return(list(kmc,
    #             res3,
    #             EvalClust,
    #             res4))
    return(list("RetainedPairs" = res3,
                "PairsByGroup" = res4))
  } else {
    # return(list(kmc,
    #             res3,
    #             EvalClust))
    return(res3)
  }
}

