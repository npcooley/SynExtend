###### -- Reject predictions that conflict with strong blocks -----------------
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

BlockReconciliation <- function(Pairs,
                                ConservativeRejection = TRUE,
                                AlwaysReject = TRUE,
                                Verbose = FALSE) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
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
  FeaturesMat <- data.frame("g1" = as.integer(FeaturesMat[, 1L]),
                            "i1" = as.integer(FeaturesMat[, 2L]),
                            "f1" = as.integer(FeaturesMat[, 3L]),
                            "g2" = as.integer(FeaturesMat[, 4L]),
                            "i2" = as.integer(FeaturesMat[, 5L]),
                            "f2" = as.integer(FeaturesMat[, 6L]),
                            "keyval" = seq(nrow(FeaturesMat)))
  
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
          sp2 <- dr3[[m3]][, 3L]
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
          sp2 <- dr4[[m3]][, 3L]
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
      
      o1 <- sapply(dr3,
                   function(x) nrow(x),
                   simplify = TRUE)
      o2 <- sapply(dr4,
                   function(x) nrow(x),
                   simplify = TRUE)
      o3 <- order(c(o1, o2),
                  decreasing = TRUE)
      dr5 <- c(dr3,
               dr4)
      
      # list of matrices that represent blocks
      # sorted by size
      dr5 <- dr5[o3]
      setID <- rep(seq(length(o3)),
                   sapply(dr5,
                          function(x) nrow(x)))
      dr5 <- do.call(rbind,
                     dr5)
      rownames(dr5) <- NULL
      dr5 <- cbind(dr5,
                   "SETID" = setID)
      # hypothetically dr5 is now set up in descending block orders
      dr5 <- dr5[!duplicated(dr5[, "ID"]), ]
      
      # loop through the setIDs that remain, removing pairs that interfer with
      # a block with higher precedence
      
      dr6 <- unname(split(x = dr5,
                          f = dr5[, "SETID"],
                          drop = FALSE))
      o1 <- sapply(dr6,
                   function(x) nrow(x),
                   simplify = TRUE)
      o2 <- order(o1,
                  decreasing = TRUE)
      dr6 <- dr6[o2]
      dr6 <- do.call(rbind,
                     dr6)
      
      # loop through the unique set ids that are present
      SETS <- dr6[, "SETID"]
      PRSETS <- unique(SETS)
      KEEP <- rep(TRUE,
                  nrow(dr6))
      F1Vec <- dr6[, 3L]
      F2Vec <- dr6[, 6L]
      if (length(PRSETS) > 1L) {
        # no need to evaluate last set
        PRSETS <- PRSETS[1:(length(PRSETS) - 1L)]
        for (m3 in seq_along(PRSETS)) {
          # which line call the current range
          ph1 <- which(SETS == PRSETS[m3])
          if (length(ph1) == 1L) {
            # exit loop upon first occurrence of a block of size 1
            break
          }
          # evaluate from the lowest row + 1L onward
          EVALSTART <- max(ph1) + 1L
          EVALEND <- nrow(dr6)
          # find bounds of current range
          Range1 <- range(F1Vec[ph1])
          Range2 <- range(F2Vec[ph1])
          if (ConservativeRejection) {
            w1 <- which(F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                          F1Vec[EVALSTART:EVALEND] <= Range1[2L] &
                          F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                          F2Vec[EVALSTART:EVALEND] <= Range2[2L])
          } else {
            w1 <- which((F1Vec[EVALSTART:EVALEND] >= Range1[1L] &
                           F1Vec[EVALSTART:EVALEND] <= Range1[2L]) |
                          (F2Vec[EVALSTART:EVALEND] >= Range2[1L] &
                             F2Vec[EVALSTART:EVALEND] <= Range2[2L]))
          }
          
          if (length(w1) > 0L) {
            KEEP[(EVALSTART:EVALEND)[w1]] <- FALSE
          }
        }
      } else {
        # do nothing no evaluation takes place, no removals take place
      }
      
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


