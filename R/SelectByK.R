###### -- drop predicted pairs based on user PID confidence -------------------
# version 2:
# users select a confidence that is the PID centroid a cluster must be at or above
# to be kept
# so you will retain pairs whose PID is below that centroid, but appear closer in distance
# to the lowest acceptable centroid than to the highest unacceptable centroid
# users can also move the ClusterScalar argument, this adjusts the decision making for
# how many clusters are going to be used for selection
# cluster number selection if performed by fitting the total within-cluster sum of squares
# to a right hyperbola / one-site binding isotherm, the half-max / Kd becomes the target of
# the scalar argument, when it is left with the default of 1, ceiling(Kd) selects the clustering regime to be used
# as scalar is adjusted, ceiling(Kd * ClusterScalar) continues to select for the clustering regime

SelectByK <- function(Pairs,
                      UserConfidence = 0.5,
                      ClusterScalar = 1,
                      MaxClusters = 15L,
                      ReturnAllCommunities = FALSE,
                      Verbose = FALSE,
                      ShowPlot = FALSE,
                      RetainHighest = TRUE) {
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
  # dat1 <- data.frame("RelativeMatch" = Pairs$ExactMatch * 2L / (Pairs$p1FeatureLength + Pairs$p2FeatureLength),
  #                    "Consensus" = Pairs$Consensus,
  #                    "PID" = Pairs$PID,
  #                    "SCORE" = Pairs$SCORE / max(Pairs$SCORE),
  #                    "TetDist" = Pairs$TetDist)
  
  # treating score differently
  # this doesn't seem to really change anything though
  s1 <- Pairs$SCORE
  s2 <- min(s1[s1 > 0])
  if (s2 > 1) {
    s2 <- 0.1
  }
  if (any(s1 < 0)) {
    s1[s1 < 0] <- s2
  }
  s1 <- log(s1)
  s2 <- max(s1)
  s1 <- s1 / s2
  dat1 <- data.frame("RelativeMatch" = Pairs$ExactMatch * 2L / (Pairs$p1FeatureLength + Pairs$p2FeatureLength),
                     "Consensus" = Pairs$Consensus,
                     "PID" = Pairs$PID,
                     # "SCORE" = Pairs$SCORE,
                     # "TetDist" = Pairs$TetDist,
                     # "SeqDiff" = abs(Pairs$p1FeatureLength - Pairs$p2FeatureLength),
                     "SCORE" = s1)
  
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
  
  OneSite <- function(X,
                      Bmax,
                      Kd) {
    Y <- (Bmax * X) / (Kd + X)
  }
  
  # first column defaults to the x axis
  n <- NClust
  dat2 <- cbind("n" = n,
                "wss" = wss)
  # transform data to make the 2 cluster data the origin
  dat3 <- cbind("n" = n - 2L,
                "wss" = abs(wss - wss[1]))
  # plot(dat3)
  
  FitA <- nls(dat3[, 2L]~OneSite(X = dat3[, 1L],
                                 Bmax,
                                 Kd),
              start = list(Bmax = max(dat3[, 2L]),
                           Kd = unname(quantile(dat3[, 1L], .25))))
  FitASum <- summary(FitA)
  FitASum$coefficients
  
  # fit is offset by -2L to plot and fit correctly, re-offset by +1 to select the correct
  # list position
  EvalClust <- ceiling((FitASum$coefficients["Kd", "Estimate"] + 1L) * ClusterScalar)
  if (EvalClust >= MaxClusters) {
    warning("Evaluated clusters may be insufficient for this task. Retry with a larger 'MaxClusters'.")
    EvalClust <- MaxClusters
  }
  if (EvalClust < 1L) {
    warning("Scalar selection requested a number of clusters less than 2, defaulting to 2 clusterings.")
    EvalClust <- 1L
  }
  
  SEQ <- seq(from = dat3[1, 1L],
             to = dat3[nrow(dat3), 1L],
             by = 0.05)
  CURVE <- OneSite(X = SEQ,
                   Kd = FitASum$coefficients[2, 1],
                   Bmax = FitASum$coefficients[1, 1])
  
  # original method, select by quantile
  # select the cluster number based on quantile of within cluster sum of squares
  # that users input, default is 50%, seems to be fine for tests of initial version
  # EvalClust <- max(which(wss >= quantile(wss, ClusterSelect))) + 1L
  # if (EvalClust >= MaxClusters) {
  #   warning("Evaluated clusters may be insufficient for this task. Retry with a larger 'MaxClusters'.")
  #   EvalClust <- MaxClusters
  # }
  # EvalClust <- max(which(diff(wss) <= quantile(diff(wss), ClusterSelect))) + 1L
  
  # print(wss)
  # print(quantile(wss, ClusterSelect))
  # print(EvalClust)
  
  res1 <- kmc[[EvalClust]]$centers[, "PID"]
  res1 <- names(res1)[res1 >= UserConfidence]
  if (RetainHighest & length(res1) == 0) {
    res1 <- names(which.max(kmc[[EvalClust]]$centers[, "PID"]))
  }
  
  res2 <- kmc[[EvalClust]]$cluster %in% as.integer(res1)
  
  res3 <- Pairs[res2, ]
  res4 <- vector(mode = "list",
                 length = nrow(kmc[[EvalClust]]$centers))
  for (m2 in seq_along(res4)) {
    # kmeans increments clusters by 1 from 1 by default, though numerically
    # they may not be ordered in increasing PID
    res4[[m2]] <- Pairs[kmc[[EvalClust]]$cluster %in% m2, ]
    attr(res4[[m2]], "GeneCalls") <- NULL
  }
  if (ShowPlot) {
    layout(mat = matrix(data = c(1,2,3,4),
                        nrow = 2,
                        byrow = TRUE))
    par(mar = c(2.8,2.8,2.8,1),
        mgp = c(1.55,0.5,0))
    plot(x = dat3[, 1L],
         y = dat3[, 2L],
         main = "fit Bmax Kd",
         xaxt = "n",
         xlab = "Cluster number",
         ylab = "Transformed WSS val")
    axis(side = 1,
         at = seq(from = 0,
                  to = max(dat3[, 1L]),
                  by = 1L),
         labels = seq(from = 0,
                      to = max(dat3[, 1L]),
                      by = 1L) + 2L)
    lines(x = SEQ,
          y = CURVE)
    abline(v = EvalClust - 1,
           lty = 2,
           col = "blue")
    plot(density(Pairs$PID),
         main = "density")
    plot(x = 0,
         y = 0,
         type = "n",
         ylab = "Frequency",
         xlab = "PID",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = "Cluster CDFs")
    for (m2 in seq_along(res4)) {
      u1 <- res4[[m2]]$PID
      u2 <- mean(u1)
      points(x = sort(u1),
             y = seq_along(u1) / length(u1),
             col = m2,
             pch = if (u2 >= UserConfidence) {
               1
             } else {
               4
             })
    }
    resgroups <- kmc[[EvalClust]]$cluster
    pchkey1 <- tapply(X = Pairs$PID,
                      INDEX = resgroups,
                      FUN = mean)
    pchkey2 <- names(pchkey1)[pchkey1 >= UserConfidence]
    if (RetainHighest & length(pchkey2) == 0) {
      pchkey2 <- names(which.max(kmc[[EvalClust]]$centers[, "PID"]))
    }
    pchkey2 <- as.integer(pchkey2)
    plot(x = Pairs$PID,
         y = Pairs$SCORE,
         xlim = c(0, 1),
         ylim = range(Pairs$SCORE),
         col = resgroups,
         pch = ifelse(test = resgroups %in% pchkey2,
                      yes = 1,
                      no = 4),
         main = "Pairs",
         xlab = "PID",
         ylab = "SCORE")
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

