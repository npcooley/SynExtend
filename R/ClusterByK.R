###### -- k means based clustering of predicted pairs -------------------------
# author: nicholas cooley
# contact: npc19@pitt.edu

ClusterByK <- function(SynExtendObject,
                       UserConfidence = list("PID" = 0.3),
                       ClusterScalar = 4,
                       MaxClusters = 15L,
                       ColSelect = c("p1featurelength",
                                     "p2featurelength",
                                     "TotalMatch",
                                     "Consensus",
                                     "PID",
                                     "Score"),
                       ColNorm = "Score",
                       ShowPlot = FALSE,
                       Verbose = FALSE) {
  # start with timing
  if (Verbose) {
    pBar <- txtProgressBar(style = 1)
    FunctionTimeStart <- Sys.time()
    PBAR <- MaxClusters - 1L
  }
  
  # normalize score
  # NormScore <- function(vec) {
  #   s1 <- min(vec[vec > 0])
  #   if (s1 > 1) {
  #     s1 <- 0.1
  #   }
  #   if (any(vec <= 0)) {
  #     vec[vec <= 0] <- s1
  #   }
  #   ret <- log(vec)
  #   ret <- ret / max(ret)
  #   return(ret)
  # }
  NormScore <- function(vec) {
    vec <- vec / sqrt(sum(vec^2))
    return(vec)
  }
  
  # right hyperbola fit
  OneSite <- function(X,
                      Bmax,
                      Kd) {
    Y <- (Bmax * X) / (Kd + X)
  }
  
  # overhead checking
  if (!is(object = SynExtendObject,
          class2 = "PairSummaries")) {
    stop ("SynExtendObject must be an object of class 'PairSummaries'.")
  }
  if (nrow(SynExtendObject) < MaxClusters) {
    MaxClusters <- nrow(SynExtendObject)
  }
  if (MaxClusters < 2) {
    cat("Too few rows present.\n")
    SynExtendObject <- cbind(SynExtendObject,
                             "ClusterID" = 1L)
    
    return(SynExtendObject)
  }
  if (length(UserConfidence) != 1L |
      is.null(names(UserConfidence)) |
      length(UserConfidence[[1]]) != 1L) {
    stop ("UserConfidence must be a named list of length 1 with only 1 value.")
  }
  if (!(names(UserConfidence) %in% colnames(SynExtendObject))) {
    stop ("UserConfidence must be named after a column present in the 'SynExtendObject' object.")
  }
  UserConfSelect <- names(UserConfidence)
  UserConfVal <- UserConfidence[[1L]]
  if (length(ColSelect) < 1L) {
    stop ("ColSelect must be a character vector of appropriate column names.")
  }
  if (any(!(ColSelect %in% colnames(SynExtendObject)))) {
    stop ("ColSelect must be a character vector of appropriate column names.")
  }
  
  # create dat1 from the available columns
  # if total match is present with the feature lengths
  # create a relative match feature
  if (all(c("p1featurelength",
            "p2featurelength",
            "TotalMatch") %in% ColSelect)) {
    ColSelect <- ColSelect[!(ColSelect %in% c("p1featurelength",
                                              "p2featurelength",
                                              "TotalMatch"))]
    GreaterMatch <- SynExtendObject$TotalMatch / ifelse(test = SynExtendObject$p1featurelength > SynExtendObject$p2featurelength,
                                                        yes = SynExtendObject$p1featurelength,
                                                        no = SynExtendObject$p2featurelength)
    RelMatch <- (SynExtendObject$TotalMatch * 2L) / (SynExtendObject$p1featurelength + SynExtendObject$p2featurelength)
    dat1 <- cbind(# "RelMatch" = RelMatch,
                  "GreaterMatch" = GreaterMatch,
                  SynExtendObject[, names(SynExtendObject) %in% ColSelect])
  } else {
    dat1 <- SynExtendObject[, names(SynExtendObject) %in% ColSelect]
  }
  
  # normalize columns if requested
  if (!is.null(ColNorm)) {
    for (m1 in seq_along(ColNorm)) {
      if (ColNorm[m1] %in% colnames(dat1)) {
        dat1[[ColNorm[m1]]] <- NormScore(dat1[[ColNorm[m1]]])
      }
    }
  }
  
  # the old function just truncates the object, here we're going to just append
  # a new column and an attribute of PID means
  # cluster out to user specified total clusters
  # run the k means clusterings
  NClust <- 2:MaxClusters
  kmc <- vector(mode = "list",
                length = length(NClust))
  for (m1 in seq_along(NClust)) {
    kmc[[m1]] <- suppressWarnings(kmeans(x = dat1,
                                         centers = NClust[m1],
                                         iter.max = 25L,
                                         nstart = 25L))
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / PBAR)
    }
  }
  if (Verbose) {
    close(pBar)
  }
  
  # get the within cluster sum of squares
  wss <- sapply(X = kmc,
                FUN = function(x) {
                  x$tot.withinss
                })
  # first column defaults to the x axis
  n <- NClust
  dat2 <- cbind("n" = n,
                "wss" = wss)
  # transform data to make the 2 cluster data the origin
  dat3 <- cbind("n" = n - 2L,
                "wss" = abs(wss - wss[1]))
  
  FitA <- nls(dat3[, 2L]~OneSite(X = dat3[, 1L],
                                 Bmax,
                                 Kd),
              start = list(Bmax = max(dat3[, 2L]),
                           Kd = unname(quantile(dat3[, 1L], .25))))
  FitASum <- summary(FitA)
  
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
  
  UserW <- which(colnames(kmc[[EvalClust]]$centers) == UserConfSelect)
  res1 <- kmc[[EvalClust]]$centers[, UserW]
  res4 <- vector(mode = "list",
                 length = nrow(kmc[[EvalClust]]$centers))
  for (m2 in seq_along(res4)) {
    # kmeans increments clusters by 1 from 1 by default, though numerically
    # they may not be ordered in increasing PID
    res4[[m2]] <- SynExtendObject[kmc[[EvalClust]]$cluster %in% m2, ]
    attr(res4[[m2]], "GeneCalls") <- NULL
  }
  
  res <- cbind(SynExtendObject, 
               "ClusterID" = kmc[[EvalClust]]$cluster)
  
  attr(res,
       "UserCriteria") <- UserConfidence
  attr(res,
       "centers") <- kmc[[EvalClust]]$centers
  attr(res,
       "Retain") <- res1 >= UserConfVal
  attr(res,
       "GeneCalls") <- attr(SynExtendObject,
                            "GeneCalls")
  attr(x = res,
       which = "AlignmentFunction") <- attr(x = SynExtendObject,
                                            which = "AlignmentFunction")
  attr(x = res,
       which = "UserConfidence") <- UserConfidence
  attr(x = res,
       which = "KVal") <- attr(x = SynExtendObject,
                               which = "KVal")
  class(res) <- c("data.frame",
                  "PairSummaries")
  
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
    hist(res$PID,
         breaks = seq(from = 0,
                      to = 1,
                      by = 0.01))
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
      lines(x = sort(u1),
            y = seq_along(u1) / length(u1),
            col = m2,
            lty = if (u2 >= UserConfidence) {
              1
            } else {
              4
            })
    }
    resgroups <- kmc[[EvalClust]]$cluster
    pchkey1 <- tapply(X = res$PID,
                      INDEX = resgroups,
                      FUN = mean)
    pchkey2 <- names(pchkey1)[pchkey1 >= UserConfidence]
    pchkey2 <- as.integer(pchkey2)
    plot(x = res$PID,
         y = res$Score,
         xlim = c(0, 1),
         ylim = range(res$Score),
         col = resgroups,
         pch = ifelse(test = resgroups %in% pchkey2,
                      yes = 1,
                      no = 4),
         main = "Pairs",
         xlab = "PID",
         ylab = "SCORE")
  }
  
  if (Verbose) {
    FunctionTimeEnd <- Sys.time()
    FunctionTimeTotal <- FunctionTimeEnd - FunctionTimeStart
    print(FunctionTimeTotal)
  }
  
  return(res)
  
}

