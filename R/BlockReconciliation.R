###### -- Reject predictions that conflict with strong blocks -----------------
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

BlockReconciliation <- function(Pairs,
                                ConservativeRejection = TRUE,
                                Precedent = "Size",
                                PIDThreshold = NULL,
                                SCOREThreshold = NULL,
                                Verbose = FALSE) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  if (!is(object = Pairs,
          class2 = "PairSummaries")) {
    stop ("Pairs must be an object of class 'PairSummaries'.")
  }
  if (!is.null(PIDThreshold)) {
    if (length(PIDThreshold) > 1L) {
      stop ("PIDThreshold must be a numeric of length 1.")
    }
  }
  if (!is.null(SCOREThreshold)) { 
    if (length(SCOREThreshold) > 1L) {
      stop ("SCOREThreshold must be a numeric of length. 1")
    }
  }
  if (!is.null(PIDThreshold)) {
    if (is.na(PIDThreshold)) {
      stop ("PIDThreshold cannot be NA.")
    }
  }
  if (!is.null(SCOREThreshold)) {
    if (is.na(SCOREThreshold)) {
      stop ("SCOREThreshold cannot be NA.")
    }
  }
  if (length(Precedent) != 1L) {
    stop ("Precedent must be a character vector of either 'Size' or 'Metric'.")
  }
  if (!is.character(Precedent)) {
    stop ("Precedent must be a character vector of either 'Size' or 'Metric'.")
  }
  if (!(Precedent %in% c("Size",
                         "Metric"))) {
    stop ("Precedent must be a character vector of either 'Size' or 'Metric'.")
  }
  # for each genome to genome comparison
  # build a matrix in the form:
  # rowval | colval | diagval | antidiagval | PID/PredictedPID | SIMPredictedSIM
  
  POIDs <- paste(Pairs$p1,
                 Pairs$p2,
                 sep = "_")
  FeaturesMat <- do.call(rbind,
                         strsplit(x = POIDs,
                                  split = "_",
                                  fixed = TRUE))
  if (is.null(PIDThreshold) &
      is.null(SCOREThreshold)) {
    FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                              "i1" = as.integer(FeaturesMat[, 2L]),
                              "f1" = as.integer(FeaturesMat[, 3L]),
                              "g2" = as.integer(FeaturesMat[, 4L]),
                              "i2" = as.integer(FeaturesMat[, 5L]),
                              "f2" = as.integer(FeaturesMat[, 6L]),
                              "keyval" = seq(nrow(FeaturesMat)))
  } else if (is.null(PIDThreshold) &
             !is.null(SCOREThreshold) &
             "SCORE" %in% colnames(Pairs)) {
    FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                              "i1" = as.integer(FeaturesMat[, 2L]),
                              "f1" = as.integer(FeaturesMat[, 3L]),
                              "g2" = as.integer(FeaturesMat[, 4L]),
                              "i2" = as.integer(FeaturesMat[, 5L]),
                              "f2" = as.integer(FeaturesMat[, 6L]),
                              "keyval" = seq(nrow(FeaturesMat)),
                              "SCORE" = Pairs$SCORE)
  } else if (!is.null(PIDThreshold) &
             is.null(SCOREThreshold) &
             "PID" %in% colnames(Pairs)) {
    FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                              "i1" = as.integer(FeaturesMat[, 2L]),
                              "f1" = as.integer(FeaturesMat[, 3L]),
                              "g2" = as.integer(FeaturesMat[, 4L]),
                              "i2" = as.integer(FeaturesMat[, 5L]),
                              "f2" = as.integer(FeaturesMat[, 6L]),
                              "keyval" = seq(nrow(FeaturesMat)),
                              "PID" = Pairs$PID)
  } else if (!is.null(PIDThreshold) &
             !is.null(SCOREThreshold) &
             "PID" %in% colnames(Pairs) &
             "SCORE" %in% colnames(Pairs)) {
    FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                              "i1" = as.integer(FeaturesMat[, 2L]),
                              "f1" = as.integer(FeaturesMat[, 3L]),
                              "g2" = as.integer(FeaturesMat[, 4L]),
                              "i2" = as.integer(FeaturesMat[, 5L]),
                              "f2" = as.integer(FeaturesMat[, 6L]),
                              "keyval" = seq(nrow(FeaturesMat)),
                              "PID" = Pairs$PID,
                              "SCORE" = Pairs$SCORE)
  } else {
    FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                              "i1" = as.integer(FeaturesMat[, 2L]),
                              "f1" = as.integer(FeaturesMat[, 3L]),
                              "g2" = as.integer(FeaturesMat[, 4L]),
                              "i2" = as.integer(FeaturesMat[, 5L]),
                              "f2" = as.integer(FeaturesMat[, 6L]),
                              "keyval" = seq(nrow(FeaturesMat)))
  }
  
  
  GMat <- unique(FeaturesMat[, c(1, 4)])
  rownames(GMat) <- NULL
  
  Res <- vector(mode = "list",
                length = nrow(GMat))
  
  for (m1 in seq(nrow(GMat))) {
    # subset to current genomes
    w1 <- GMat[m1, 1L]
    w2 <- GMat[m1, 2L]
    
    w3 <- FeaturesMat[, 1L] == w1
    w4 <- FeaturesMat[, 4L] == w2
    CMat <- FeaturesMat[w3 & w4, ]
    rownames(CMat) <- NULL
    # separate index pairs
    # blocks cannot span indices
    IMat <- split(x = CMat,
                  f = list(CMat$i1,
                           CMat$i2),
                  drop = TRUE)
    Res[[m1]] <- vector(mode = "list",
                        length = length(IMat))
    for (m2 in seq_along(IMat)) {
      
      # build in the diagonal ranks
      # munge the data such that the matrix is sorted by largest
      # block to smallest block,
      # i.e. largest continuous diag rank
      dr1 <- IMat[[m2]][, 3L] - IMat[[m2]][, 6L]
      dr2 <- IMat[[m2]][, 3L] + IMat[[m2]][, 6L]
      IMat[[m2]] <- cbind(IMat[[m2]],
                          "ID" = seq(nrow(IMat[[m2]])),
                          "rank1" = dr1,
                          "rank2" = dr2)
      # go through both rankings and break out non-contiguous blocks
      dr3 <- unname(split(x = IMat[[m2]],
                          f = dr1,
                          drop = TRUE))
      dr4 <- unname(split(x = IMat[[m2]],
                          f = dr2,
                          drop = TRUE))
      for (m3 in seq_along(dr3)) {
        # if the current rank has more than one pair
        if (nrow(dr3[[m3]]) > 1L) {
          # build a dummy vector for the split
          sp1 <- vector(mode = "integer",
                        length = nrow(dr3[[m3]]))
          # extract the feature positions that are always increasing
          sp2 <- dr3[[m3]][, "f1"]
          # construct iterators
          it1 <- 1L
          it2 <- sp2[1L]
          # loop through dummy vector, and assign the split iterator
          # updating it when a gap larger than the tolerance appears
          for (m4 in seq_along(sp1)) {
            it3 <- sp2[m4]
            if (it3 - it2 > 1L) {
              it1 <- it1 + 1L
            }
            sp1[m4] <- it1
            it2 <- it3
          } # end loop through split map
          
          if (it1 > 1L) {
            # the split map was updated, split the matrix
            dr3[[m3]] <- unname(split(x = dr3[[m3]],
                                      f = sp1))
          } else {
            dr3[[m3]] <- dr3[m3]
          }
          
        } else {
          dr3[[m3]] <- dr3[m3]
        }
      }
      
      for (m3 in seq_along(dr4)) {
        # if the current rank has more than one pair
        if (nrow(dr4[[m3]]) > 1L) {
          # build a dummy vector for the split
          sp1 <- vector(mode = "integer",
                        length = nrow(dr4[[m3]]))
          # extract the feature positions that are always increasing
          sp2 <- dr4[[m3]][, "f1"]
          # construct iterators
          it1 <- 1L
          it2 <- sp2[1L]
          # loop through dummy vector, and assign the split iterator
          # updating it when a gap larger than the tolerance appears
          for (m4 in seq_along(sp1)) {
            it3 <- sp2[m4]
            if (it3 - it2 > 1L) {
              it1 <- it1 + 1L
            }
            sp1[m4] <- it1
            it2 <- it3
          } # end loop through split map
          
          if (it1 > 1L) {
            # the split map was updated, split the matrix
            dr4[[m3]] <- unname(split(x = dr4[[m3]],
                                      f = sp1))
          } else {
            dr4[[m3]] <- dr4[m3]
          }
          
        } else {
          dr4[[m3]] <- dr4[m3]
        }
      }
      
      # lists now split into blocks
      dr3 <- unlist(dr3,
                    recursive = FALSE)
      dr4 <- unlist(dr4,
                    recursive = FALSE)
      dr5 <- c(dr3,
               dr4)
      
      o1 <- sapply(dr5,
                   function(x) nrow(x),
                   simplify = TRUE)
      o2 <- order(o1,
                  decreasing = TRUE)
      dr5 <- dr5[o2]
      # o1 <- sapply(dr3,
      #              function(x) nrow(x),
      #              simplify = TRUE)
      # o2 <- sapply(dr4,
      #              function(x) nrow(x),
      #              simplify = TRUE)
      # o3 <- order(c(o1, o2),
      #             decreasing = TRUE)
      # dr5 <- c(dr3,
      #          dr4)
      # 
      # # list of matrices that represent blocks
      # # sorted by size
      # dr5 <- dr5[o3]
      # setID <- rep(seq(length(o3)),
      #              sapply(dr5,
      #                     function(x) nrow(x)))
      # setID <- rep(seq(length(o1)),
      #              sapply(dr5,
      #                     function(x) nrow(x)))
      setID <- rep(seq(length(o1)),
                   o1[o2])
      dr5 <- do.call(rbind,
                     dr5)
      rownames(dr5) <- NULL
      dr5 <- cbind(dr5,
                   "SETID" = setID)
      # hypothetically dr5 is now set up in descending block orders
      # drop duplicated appearances of the same line,
      # i.e. lines can only be evaluated against blocks that are larger than they are
      dr5 <- dr5[!duplicated(dr5[, "ID"]), ]
      
      # loop through the setIDs that remain, removing pairs that interfere with
      # a block with higher precedence
      
      dr6 <- unname(split(x = dr5,
                          f = dr5[, "SETID"],
                          drop = FALSE))
      # reset setID with probable new block sizes and with additional sorting
      # from the mean PID / SCORE of the block IF PID / SCORE is available
      if ("PID" %in% colnames(dr5)) {
        o1 <- sapply(dr6,
                     function(x) nrow(x),
                     simplify = TRUE)
        o2 <- sapply(X = dr6,
                     FUN = function(x) {
                       mean(x$PID)
                     },
                     simplify = TRUE)
        if (Precedent == "Size") {
          o3 <- order(o1,
                      o2,
                      decreasing = TRUE)
        } else if (Precedent == "Metric") {
          o3 <- order(o2,
                      o1,
                      decreasing = TRUE)
        }
        dr6 <- dr6[o3]
        setID <- rep(seq(length(o3)),
                     o1[o3])
      } else if ("SCORE" %in% colnames(dr5)) {
        o1 <- sapply(dr6,
                     function(x) nrow(x),
                     simplify = TRUE)
        o2 <- sapply(X = dr6,
                     FUN = function(x) {
                       mean(x$SCORE)
                     },
                     simplify = TRUE)
        if (Precedent == "Size") {
          o3 <- order(o1,
                      o2,
                      decreasing = TRUE)
        } else if (Precedent == "Metric") {
          o3 <- order(o2,
                      o1,
                      decreasing = TRUE)
        }
        dr6 <- dr6[o3]
        setID <- rep(seq(length(o3)),
                     o1[o3])
      } else {
        o1 <- sapply(dr6,
                     function(x) nrow(x),
                     simplify = TRUE)
        o2 <- order(o1,
                    decreasing = TRUE)
        dr6 <- dr6[o2]
        setID <- rep(seq(length(o2)),
                     o1[o2])
      }
      
      # if (Precedent == "Metric") {
      #   dr6 <- dr6[sapply(X = dr6,
      #                     FUN = function(x) {
      #                       nrow(x) > 1L
      #                     },
      #                     simplify = TRUE)]
      # }
      
      dr6 <- do.call(rbind,
                     dr6)
      dr6[, "SETID"] <- setID
      # return(dr6)
      # loop through the unique set ids that are present
      SETS <- dr6[, "SETID"]
      PRSETS <- unique(SETS)
      KEEP <- rep(TRUE,
                  nrow(dr6))
      F1Vec <- dr6[, "f1"]
      F2Vec <- dr6[, "f2"]
      if ("SCORE" %in% colnames(dr6)) {
        F3Vec <- dr6[, "SCORE"]
      }
      if ("PID" %in% colnames(dr6)) {
        F4Vec <- dr6[, "PID"]
      }
      # return(dr6)
      if (length(PRSETS) > 1L) {
        # no need to evaluate last set
        PRSETS <- PRSETS[seq_len(length(PRSETS) - 1L)]
        for (m3 in seq_along(PRSETS)) {
          # which line - call the current range
          ph1 <- which(SETS == PRSETS[m3])
          if (length(ph1) == 1L) {
            # exit loop upon first occurrence of a block of size 1
            if (Precedent == "Size") {
              break
            } else {
              next
            }
          }
          # evaluate from the lowest row + 1L onward
          EVALSTART <- max(ph1) + 1L
          EVALEND <- nrow(dr6)
          # find bounds of current range
          Range1 <- range(F1Vec[ph1])
          Range2 <- range(F2Vec[ph1])
          if (ConservativeRejection) {
            # reject within the block only
            if ("SCORE" %in% colnames(dr6) &
                "PID" %in% colnames(dr6)) {
              # both PID and SCORE thresholds
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which(F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                            F1Vec[EVALSTART:EVALEND] <= Range1[2L] &
                            F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                            F2Vec[EVALSTART:EVALEND] <= Range2[2L] &
                            (F3Vec[EVALSTART:EVALEND] < SCOREThreshold |
                               F4Vec[EVALSTART:EVALEND] < PIDThreshold))
            } else if ("SCORE" %in% colnames(dr6) &
                       !("PID" %in% colnames(dr6))) {
              # SCORE thresholds only
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which(F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                            F1Vec[EVALSTART:EVALEND] <= Range1[2L] &
                            F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                            F2Vec[EVALSTART:EVALEND] <= Range2[2L] &
                            F3Vec[EVALSTART:EVALEND] < SCOREThreshold)
            } else if (!("SCORE" %in% colnames(dr6)) &
                       "PID" %in% colnames(dr6)) {
              # PID thresholds only
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which(F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                            F1Vec[EVALSTART:EVALEND] <= Range1[2L] &
                            F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                            F2Vec[EVALSTART:EVALEND] <= Range2[2L] &
                            F4Vec[EVALSTART:EVALEND] < PIDThreshold)
            } else {
              # no thresholds
              w1 <- which(F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                            F1Vec[EVALSTART:EVALEND] <= Range1[2L] &
                            F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                            F2Vec[EVALSTART:EVALEND] <= Range2[2L])
            }
            
          } else {
            if ("SCORE" %in% colnames(dr6) &
                "PID" %in% colnames(dr6)) {
              # both PID and SCORE thresholds
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which((F3Vec[EVALSTART:EVALEND] < SCOREThreshold |
                             F4Vec[EVALSTART:EVALEND] < PIDThreshold) &
                            ((F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                                F1Vec[EVALSTART:EVALEND] <= Range1[2L]) |
                               (F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                                  F2Vec[EVALSTART:EVALEND] <= Range2[2L])))
              # if (any(dr6[(EVALSTART:EVALEND)[w1], "keyval"] %in% c(613, 615))) {
              #   return(list(dr6,
              #               w1,
              #               Range1,
              #               Range2,
              #               F1Vec,
              #               F2Vec,
              #               F3Vec,
              #               F4Vec,
              #               EVALSTART,
              #               EVALEND))
              # }
              
            } else if ("SCORE" %in% colnames(dr6) &
                       !("PID" %in% colnames(dr6))) {
              # SCORE thresholds only
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which(F3Vec[EVALSTART:EVALEND] < SCOREThreshold &
                            ((F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                                F1Vec[EVALSTART:EVALEND] <= Range1[2L]) |
                               (F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                                  F2Vec[EVALSTART:EVALEND] <= Range2[2L])))
            } else if (!("SCORE" %in% colnames(dr6)) &
                       "PID" %in% colnames(dr6)) {
              # PID thresholds only
              # w1 is the rejection set, so vals below the thresholds are TRUE
              w1 <- which(F4Vec[EVALSTART:EVALEND] < PIDThreshold &
                            ((F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                                F1Vec[EVALSTART:EVALEND] <= Range1[2L]) |
                               (F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                                  F2Vec[EVALSTART:EVALEND] <= Range2[2L])))
            } else {
              # no thresholds
              w1 <- which((F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                             F1Vec[EVALSTART:EVALEND] <= Range1[2L]) |
                            (F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                               F2Vec[EVALSTART:EVALEND] <= Range2[2L]))
            }
          }
          
          if (length(w1) > 0L) {
            KEEP[(EVALSTART:EVALEND)[w1]] <- FALSE
          }
        }
      } else {
        # do nothing no evaluation takes place, no removals take place
      }
      # return(list(dr6,
      #             KEEP))
      dr7 <- dr6[KEEP, ]
      Res[[m1]][[m2]] <- dr7[, "keyval"]
      # return(list(IMat[[m2]],
      #             dr3,
      #             dr4,
      #             dr5,
      #             dr6,
      #             dr7))
      
    } # end loop through index pairs within genome pairs
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / nrow(GMat))
    }
    
  } # end loop through genome-genome comparisons
  
  if (Verbose) {
    cat("\n")
    close(pBar)
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  Res <- unlist(Res,
                recursive = TRUE)
  Res <- Pairs[sort(Res), ]
  rownames(Res) <- NULL
  return(Res)
}


