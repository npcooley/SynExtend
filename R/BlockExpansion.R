###### -- Expand blocks of predicted pairs ------------------------------------
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

BlockExpansion <- function(Pairs,
                           GapTolerance = 4L,
                           DropSingletons = FALSE,
                           Criteria = "PID",
                           Floor = 0.5,
                           NewPairsOnly = TRUE,
                           DBPATH,
                           Verbose = FALSE) {
  
  if (!("PID" %in% colnames(Pairs))) {
    stop ("PairSummaries Object must have PIDs calculated.")
  }
  if (missing(DBPATH)) {
    stop("DBPATH must be supplied.")
  }
  # given the gap tolerance, split diags where appropriate
  if (is.null(GapTolerance)) {
    GapTolerance <- 2L
  }
  if (GapTolerance <= 1L) {
    stop ("GapTolerance defines the diff() of features within a block. It cannot be <= 1.")
  }
  
  if (Verbose) {
    TimeStart <- Sys.time()
  }
  
  MAT1 <- get(data("HEC_MI1",
                   package = "DECIPHER",
                   envir = environment()))
  MAT2 <- get(data("HEC_MI2",
                   package = "DECIPHER",
                   envir = environment()))
  structureMatrix <- matrix(c(0.187, -0.8, -0.873,
                              -0.8, 0.561, -0.979,
                              -0.873, -0.979, 0.221),
                            3,
                            3,
                            dimnames=list(c("H", "E", "C"),
                                          c("H", "E", "C")))
  substitutionMatrix <- matrix(c(1.5, -2.134, -0.739, -1.298,
                                 -2.134, 1.832, -2.462, 0.2,
                                 -0.739, -2.462, 1.522, -2.062,
                                 -1.298, 0.2, -2.062, 1.275),
                               nrow = 4,
                               dimnames = list(DNA_BASES, DNA_BASES))
  
  # break PairSummaries object down into a workable format
  # build overhead data in a way that makes sense
  GeneCalls <- attr(x = Pairs,
                    which = "GeneCalls")
  GCIDs <- as.integer(names(GeneCalls))
  L <- length(GeneCalls)
  L2 <- (L * (L - 1L)) / 2L
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
                            "f2" = as.integer(FeaturesMat[, 6L]))
  
  GMat <- unique(FeaturesMat[, c(1, 4)])
  rownames(GMat) <- NULL
  
  Res <- vector(mode = "list",
                length = nrow(GMat))
  
  for (m1 in seq(nrow(GMat))) {
    
    # subset to current genomes
    w1 <- GMat[m1, 1L]
    w2 <- GMat[m1, 2L]
    
    if (Verbose) {
      cat(paste0("### Genome pair: ",
                 w1,
                 " - ",
                 w2,
                 " ###\n"))
    }
    
    w3 <- FeaturesMat[, 1L] == w1
    w4 <- FeaturesMat[, 4L] == w2
    CMat <- FeaturesMat[w3 & w4, ]
    rownames(CMat) <- NULL
    # separate index pairs
    IMat <- split(x = CMat,
                  f = list(CMat$i1,
                           CMat$i2),
                  drop = TRUE)
    i1 <- which(GCIDs == w1)
    i2 <- which(GCIDs == w2)
    i3 <- unique(GeneCalls[[i1]]$Index)
    i4 <- unique(GeneCalls[[i2]]$Index)
    i5 <- sapply(i3,
                 function(x) {
                   which(GeneCalls[[i1]]$Index == x)
                 },
                 simplify = FALSE)
    i6 <- sapply(i4,
                 function(x) {
                   which(GeneCalls[[i2]]$Index == x)
                 },
                 simplify = FALSE)
    i1lower <- sapply(i5,
                      function(x) {
                        min(x)
                      },
                      simplify = TRUE)
    i1upper <- sapply(i5,
                      function(x) {
                        max(x)
                      },
                      simplify = TRUE)
    i2lower <- sapply(i6,
                      function(x) {
                        min(x)
                      },
                      simplify = TRUE)
    i2upper <- sapply(i6,
                      function(x) {
                        max(x)
                      },
                      simplify = TRUE)
    
    Genome1 <- SearchDB(dbFile = DBPATH,
                        identifier = as.character(w1),
                        nameBy = "description",
                        type = "DNAStringSet",
                        verbose = FALSE)
    
    Genome2 <- SearchDB(dbFile = DBPATH,
                        identifier = as.character(w2),
                        nameBy = "description",
                        type = "DNAStringSet",
                        verbose = FALSE)
    
    # for each Index Pair Matrix
    Res[[m1]] <- vector(mode = "list",
                        length = length(IMat))
    
    for (m2 in seq_along(IMat)) {
      
      if (Verbose) {
        cat(paste0("Index pair ",
                   IMat[[m2]][1L, 2L],
                   " - ",
                   IMat[[m2]][1L, 5L],
                   "\n"))
      }
      # current offsets
      ci1lower <- i1lower[IMat[[m2]][1L, 2L]]
      ci1upper <- i1upper[IMat[[m2]][1L, 2L]]
      ci2lower <- i2lower[IMat[[m2]][1L, 5L]]
      ci2upper <- i2upper[IMat[[m2]][1L, 5L]]
      if (Verbose) {
        cat(paste0("Feature set ",
                   ci1lower,
                   " - ",
                   ci1upper,
                   " and ",
                   ci2lower,
                   " - ",
                   ci2upper,
                   ":\n"))
      }
      
      # extract seqs for features one on the current index for w1
      w5 <- GeneCalls[[w1]]$Index == IMat[[m2]][1L, 2L]
      w6 <- GeneCalls[[w1]]$Coding[w5] & GeneCalls[[w1]]$Type[w5] == "gene"
      w7 <- unique(GeneCalls[[w1]]$Translation_Table[w5])
      w7 <- w7[!is.na(w7)]
      w7 <- w7[1L]
      z1 <- unname(GeneCalls[[w1]]$Range[w5])
      z2 <- lengths(z1)
      z1 <- unlist(z1,
                   recursive = FALSE)
      NTFeatures01 <- extractAt(x = Genome1[[IMat[[m2]][1L, 2L]]],
                                at = z1)
      
      CollapseCount <- 0L
      w <- which(z2 > 1L)
      # if no collapsing needs to occur, do not enter loop
      if (length(w) > 0L) {
        # if collapsing must take place build a placeholder of positions to remove
        # once collapsing correct positions has occurred
        remove <- vector(mode = "integer",
                         length = sum(z2[w]) - length(w))
        for (m4 in w) {
          NTFeatures01[[m4 + CollapseCount]] <- unlist(NTFeatures01[m4:(m4 + z2[m4] - 1L) + CollapseCount])
          remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
          CollapseCount <- CollapseCount + z2[m4] - 1L
        }
        # return(list(w,
        #             z2[w],
        #             remove))
        NTFeatures01[remove] <- NULL
      }
      
      FlipMe <- GeneCalls[[w1]]$Strand[w5] == 1L
      if (any(FlipMe)) {
        NTFeatures01[FlipMe] <- reverseComplement(NTFeatures01[FlipMe])
      }
      
      names(NTFeatures01) <- paste(rep(w1,
                                       sum(w5)),
                                   rep(IMat[[m2]][1L, 2L],
                                       sum(w5)),
                                   seq(from = ci1lower,
                                       to = ci1upper,
                                       by = 1L),
                                   sep = "_")
      if (!is.na(w7)) {
        AAFeatures01 <- translate(x = NTFeatures01[w6],
                                  genetic.code = getGeneticCode(id_or_name2 = w7,
                                                                as.data.frame = FALSE,
                                                                full.search = FALSE),
                                  if.fuzzy.codon = "solve")
      } else {
        AAFeatures01 <- AAStringSet()
      }
        names(AAFeatures01) <- names(NTFeatures01)[w6]
        
        Features01Match <- match(x = names(NTFeatures01),
                                 table = names(AAFeatures01))
        Features01Key <- names(NTFeatures01) %in% names(AAFeatures01)
      
      # extract seqs for features one on the current index for w2
      w5 <- GeneCalls[[w2]]$Index == IMat[[m2]][1L, 5L]
      w6 <- GeneCalls[[w2]]$Coding[w5] & GeneCalls[[w2]]$Type[w5] == "gene"
      w7 <- unique(GeneCalls[[w2]]$Translation_Table[w5])
      w7 <- w7[!is.na(w7)]
      w7 <- w7[1L]
      z1 <- unname(GeneCalls[[w2]]$Range[w5])
      z2 <- lengths(z1)
      z1 <- unlist(z1,
                   recursive = FALSE)
      NTFeatures02 <- extractAt(x = Genome2[[IMat[[m2]][1L, 5L]]],
                                at = z1)
      
      CollapseCount <- 0L
      w <- which(z2 > 1L)
      # if no collapsing needs to occur, do not enter loop
      if (length(w) > 0L) {
        # if collapsing must take place build a placeholder of positions to remove
        # once collapsing correct positions has occurred
        remove <- vector(mode = "integer",
                         length = sum(z2[w]) - length(w))
        for (m4 in w) {
          NTFeatures02[[m4 + CollapseCount]] <- unlist(NTFeatures02[m4:(m4 + z2[m4] - 1L) + CollapseCount])
          remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
          CollapseCount <- CollapseCount + z2[m4] - 1L
        }
        # return(list(w,
        #             z2[w],
        #             remove))
        NTFeatures02[remove] <- NULL
      }
      
      FlipMe <- GeneCalls[[w2]]$Strand[w5] == 1L
      if (any(FlipMe)) {
        NTFeatures02[FlipMe] <- reverseComplement(NTFeatures02[FlipMe])
      }
      
      names(NTFeatures02) <- paste(rep(w2,
                                       sum(w5)),
                                   rep(IMat[[m2]][1L, 5L],
                                       sum(w5)),
                                   seq(from = ci2lower,
                                       to = ci2upper,
                                       by = 1L),
                                   sep = "_")
      if (!is.na(w7)) {
        AAFeatures02 <- translate(x = NTFeatures02[w6],
                                  genetic.code = getGeneticCode(id_or_name2 = w7,
                                                                as.data.frame = FALSE,
                                                                full.search = FALSE),
                                  if.fuzzy.codon = "solve")
      } else {
        AAFeatures02 <- AAStringSet()
      }
        names(AAFeatures02) <- names(NTFeatures02)[w6]
        
        Features02Match <- match(x = names(NTFeatures02),
                                 table = names(AAFeatures02))
        Features02Key <- names(NTFeatures02) %in% names(AAFeatures02)
      
      # return(list(NTFeatures01,
      #             NTFeatures02,
      #             AAFeatures01,
      #             AAFeatures02,
      #             Features01Key,
      #             Features01Match,
      #             Features02Key,
      #             Features02Match))
      
      dr1 <- IMat[[m2]][, 3L] - IMat[[m2]][, 6L]
      dr2 <- IMat[[m2]][, 3L] + IMat[[m2]][, 6L]
      IMat[[m2]] <- cbind(IMat[[m2]],
                          "ID" = seq(nrow(IMat[[m2]])),
                          "rank1" = dr1,
                          "rank2" = dr2)
      
      # given the diagonal and anti-diagonal `ranks` that have been assigned
      # build 3 ledgers:
      # one for the diagonal, one for the ant-diagonal, and one for singletons
      # you can be in both the diagonal and the anti-diagonal ledger at the same tiem
      # but you cannot be in the singleton ledger and either diag or anti-diag
      
      # at this point we have 4 points describing each pair
      # 1: g1 feature ID number
      # 2: g2 features ID number
      # 3: diag rank
      # 4: anti-diag rank
      
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
            if (it3 - it2 > GapTolerance) {
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
            if (it3 - it2 > GapTolerance) {
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
      
      # return(list(IMat[[m2]],
      #             dr3,
      #             dr4))
      
      dr3 <- unlist(dr3,
                    recursive = FALSE)
      # dr3 IDs singletons
      dr3a <- unlist(sapply(dr3,
                            function(x) {
                              if (nrow(x) == 1L) {
                                x$ID
                              } else {
                                NA
                              }
                            },
                            simplify = FALSE))
      # dr3 IDs blocks
      dr3b <- unlist(sapply(dr3,
                            function(x) {
                              if (nrow(x) > 1L) {
                                x$ID
                              }
                            },
                            simplify = FALSE))
      dr4 <- unlist(dr4,
                    recursive = FALSE)
      # dr4 IDs singletons
      dr4a <- unlist(sapply(dr4,
                            function(x) {
                              if (nrow(x) == 1L) {
                                x$ID
                              } else {
                                NA
                              }
                            },
                            simplify = FALSE))
      # dr4 IDs blocks
      dr4b <- unlist(sapply(dr4,
                            function(x) {
                              if (nrow(x) > 1L) {
                                x$ID
                              }
                            },
                            simplify = FALSE))
      
      # drop dr3 singleton positions that are present in dr4
      dr3c <- dr3a %in% dr4b
      # drop dr4 singleton positions that are present in dr3
      dr4c <- dr4a %in% dr3b
      # drop singletons from one rank set that are present in a block from the
      # other set
      
      dr3 <- dr3[!dr3c]
      dr4 <- dr4[!dr4c]
      
      dr5 <- c(dr3, dr4)
      dr5 <- unique(dr5)
      if (DropSingletons) {
        checkrows <- sapply(dr5,
                            function(x) nrow(x),
                            simplify = TRUE)
        dr5 <- dr5[checkrows > 1L]
        if (length(dr5) == 0L) {
          # break out of m2 position without assigning Res[[m1]][[m2]] anything
          next
        }
      }
      
      # loop through dr5
      # if the position is a singleton pair, assign all possible expansions
      # if the position is a blocked set of pairs assign the two expansions
      
      dr6 <- vector(mode = "list",
                    length = length(dr5))
      for (m3 in seq_along(dr6)) {
        f1s <- dr5[[m3]][, 3L]
        f2s <- dr5[[m3]][, 6L]
        if (length(f1s) == 1L) {
          # singleton pair
          dr6[[m3]] <- data.frame("f1" = c(f1s - 1L, f1s - 1L, f1s + 1L, f1s + 1L),
                                  "f2" = c(f2s - 1L, f2s + 1L, f2s + 1L, f2s - 1L),
                                  "direction" = c(1L, 2L, 3L, 4L))
        } else {
          # a contiguous block of pairs
          if (length(unique(dr5[[m3]][, 8L])) == 1L) {
            # the regular diagonal
            f1s <- dr5[[m3]][, 3L]
            f2s <- dr5[[m3]][, 6L]
            f1f <- seq(from = min(f1s) - 1L,
                       to = max(f1s) + 1L,
                       by = 1L)
            f2f <- seq(from = min(f2s) - 1L,
                       to = max(f2s) + 1L,
                       by = 1L)
            f1f <- f1f[!(f1f %in% f1s)]
            f2f <- f2f[!(f2f %in% f2s)]
            
            if (length(f1f) > 2L) {
              # gaps in the block assign as no-expanding checks
              dr6[[m3]] <- data.frame("f1" = f1f,
                                      "f2" = f2f,
                                      "direction" = c(1L,
                                                      rep(0L,
                                                          length(f1f) - 2L),
                                                      3L))
            } else {
              # no gaps in the block
              dr6[[m3]] <- data.frame("f1" = f1f,
                                      "f2" = f2f,
                                      "direction" = c(1L,
                                                      3L))
            }
          } else {
            # the anti diagonal
            f1s <- dr5[[m3]][, 3L]
            f2s <- dr5[[m3]][, 6L]
            f1f <- seq(from = min(f1s) - 1L,
                       to = max(f1s) + 1L,
                       by = 1L)
            f2f <- seq(from = max(f2s) + 1L,
                       to = min(f2s) - 1L,
                       by = -1L)
            f1f <- f1f[!(f1f %in% f1s)]
            f2f <- f2f[!(f2f %in% f2s)]
            
            if (length(f1f) > 2L) {
              # gaps in the block assign as no-expanding checks
              dr6[[m3]] <- data.frame("f1" = f1f,
                                      "f2" = f2f,
                                      "direction" = c(2L,
                                                      rep(0L,
                                                          length(f1f) - 2L),
                                                      4L))
            } else {
              # no gaps in the block
              dr6[[m3]] <- data.frame("f1" = f1f,
                                      "f2" = f2f,
                                      "direction" = c(2L,
                                                      4L))
            }
          }
        } # end row check
      } # end dr6 loop
      
      dr6 <- do.call(rbind,
                     dr6)
      
      dr6 <- dr6[dr6[, 1L] >= ci1lower &
                   dr6[, 1L] <= ci1upper &
                   dr6[, 2L] >= ci2lower &
                   dr6[, 2L] <= ci2upper, , drop = FALSE]
      
      if (nrow(dr6) < 1L) {
        next
      }
      
      dr6 <- dr6[order(dr6[, 1L]), ]
      
      # return(dr6)
      # for every line in dr6
      # n + ? alignments will be attempted
      # for every n alignments that pass some threshold a new `Pair` is recorded
      # to add to pair summaries
      L <- nrow(dr6)
      VSize <- L * 2L
      
      p1placeholder <- p2placeholder <- p1FeatureLength <- p2FeatureLength <- rep(NA_integer_,
                                                                                  times = VSize)
      PIDVector <- SCOREVector <- rep(NA_real_,
                                      times = VSize)
      AType <- rep(NA_character_,
                   times = VSize)
      
      # return(list(dr6,
      #             ci1lower,
      #             ci1upper,
      #             ci2lower,
      #             ci2upper,
      #             Features01Key,
      #             Features02Key,
      #             Features01Match,
      #             Features02Match,
      #             NTFeatures01,
      #             AAFeatures01,
      #             NTFeatures02,
      #             AAFeatures02))
      
      Count <- 1L
      Continue <- TRUE
      if (Verbose) {
        pBar <- txtProgressBar(style = 1L)
      }
      
      for (m3 in seq(nrow(dr6))) {
        
        # each line contains feature coordinates for an alignment
        # and the direction to expand, if the alignment passes a threshold
        # with each expansion, check whether the expansion is within bounds
        f1 <- dr6[m3, 1L]
        f2 <- dr6[m3, 2L]
        advanceID <- dr6[m3, 3L]
        
        while (Continue) {
          # first ask if f1 and f2 are within bounds
          # if they are, do not attempt alignment
          # exit the while loop and move to the next position in dr6
          if (f1 > ci1upper |
              f1 < ci1lower |
              f2 > ci2upper |
              f2 < ci2lower) {
            # Continue <- FALSE
            break
          }
          # third ask if these sequences have already been aligned
          # CURRENTLY NOT IMPLEMENTED
          # third ask if both f1 and f2 are both coding
          if (Features01Key[f1 - ci1lower + 1L] &
              Features02Key[f2 - ci2lower + 1L]) {
            # both are coding
            
            ali <- AlignProfiles(pattern = AAFeatures01[Features01Match[f1 - ci1lower + 1L]],
                                 subject = AAFeatures02[Features02Match[f2 - ci2lower + 1L]])
            PID <- 1 - DistanceMatrix(myXStringSet = ali,
                                      type = "matrix",
                                      includeTerminalGaps = TRUE,
                                      verbose = FALSE)[1L, 2L]
            # UW <- unique(width(ali))
            SCORE <- ScoreAlignment(myXStringSet = ali,
                                    structures = PredictHEC(ali,
                                                            type="probabilities",
                                                            HEC_MI1 = MAT1,
                                                            HEC_MI2 = MAT2),
                                    structureMatrix = structureMatrix)
            CType <- "AA"
            # if (m2 >= 2L) {
            #   print(f1)
            #   print(f2)
            #   print(PID)
            #   print(CType)
            # }
          } else {
            # at least one is not coding
            ali <- AlignProfiles(pattern = NTFeatures01[f1 - ci1lower + 1L],
                                 subject = NTFeatures02[f2 - ci2lower + 1L])
            PID <- 1 - DistanceMatrix(myXStringSet = ali,
                                      type = "matrix",
                                      includeTerminalGaps = TRUE,
                                      verbose = FALSE)[1L, 2L]
            # UW <- unique(width(ali))
            SCORE <- ScoreAlignment(myXStringSet = ali,
                                    substitutionMatrix = substitutionMatrix)
            CType <- "NT"
            # if (m2 >= 2L) {
            #   print(f1)
            #   print(f2)
            #   print(PID)
            #   print(CType)
            # }
          } # end if else statement for coding / non
          # assign the dist 
          # update the count
          if (Criteria == "PID") {
            CHECK <- PID
          } else if (Criteria == "Score") {
            CHECK <- SCORE
          }
          if (CHECK < Floor) {
            # the current PID does not meet inclusion criteria
            # exit the while loop and move to the next position in dr6
            Continue <- FALSE
          } else {
            # PID meets inclusion criteria
            # assign integers to result vectors
            p1placeholder[Count] <- f1
            p2placeholder[Count] <- f2
            p1FeatureLength[Count] <- width(NTFeatures01[f1 - ci1lower + 1L])
            p2FeatureLength[Count] <- width(NTFeatures02[f2 - ci2lower + 1L])
            PIDVector[Count] <- PID
            SCOREVector[Count] <- SCORE
            AType[Count] <- CType
            
            # update f1 and f2
            if (advanceID == 0L) {
              # interior alignment, no forward movement allowed
              Continue <- FALSE
            } else if (advanceID == 1L) {
              # advancing down the diagonal
              f1 <- f1 - 1L
              f2 <- f2 - 1L
            } else if (advanceID == 2L) {
              # advancing up the anti-diagonal
              f1 <- f1 - 1L
              f2 <- f2 + 1L
            } else if (advanceID == 3L) {
              # advancing up the diagonal
              f1 <- f1 + 1L
              f2 <- f2 + 1L
            } else if (advanceID == 4L) {
              # advancing down the anti-diagonal
              f1 <- f1 + 1L
              f2 <- f2 - 1L
            } # end advancement if elses 
            
            # if (f1 == 1572 &
            #     f2 == 5750) {
            #   return(list(f1,
            #               f2,
            #               Count,
            #               m3,
            #               p1placeholder[Count],
            #               p2placeholder[Count],
            #               dr6))
            # }
            
            Count <- Count + 1L
            if (Count >= VSize) {
              # if Count exceeds VSize, increase size
              VSize <- VSize * 2L
              p1placeholder <- c(p1placeholder,
                                 rep(NA_integer_,
                                     times = VSize))
              p2placeholder <- c(p2placeholder,
                                 rep(NA_integer_,
                                     times = VSize))
              p1FeatureLength <- c(p1FeatureLength,
                                   rep(NA_integer_,
                                       times = VSize))
              p2FeatureLength <- c(p2FeatureLength,
                                   rep(NA_integer_,
                                       times = VSize))
              PIDVector <- c(PIDVector,
                             rep(NA_real_,
                                 times = VSize))
              AType <- c(AType,
                         rep(NA_character_,
                             times = VSize))
            } # end count size if statement
            
          } # end PID check
        } # end while loop
        Continue <- TRUE
        
        if (Verbose) {
          setTxtProgressBar(pb = pBar,
                            value = m3 / L)
        }
      } # end of m3 loop through dr6
      
      if (Verbose) {
        cat("\n")
        close(pBar)
      }
      
      if (any(!is.na(PIDVector))) {
        L2 <- max(which(!is.na(PIDVector)))
        
        p1placeholder <- p1placeholder[1:L2]
        p2placeholder <- p2placeholder[1:L2]
        p1FeatureLength <- p1FeatureLength[1:L2]
        p2FeatureLength <- p2FeatureLength[1:L2]
        PIDVector <- PIDVector[1:L2]
        SCOREVector <- SCOREVector[1:L2]
        AType <- AType[1:L2]
        
        # if (m2 == 2L) {
        #   return(list(IMat,
        #               w1,
        #               w2,
        #               L2,
        #               p1placeholder,
        #               p2placeholder,
        #               p1FeatureLength,
        #               p2FeatureLength,
        #               PIDVector,
        #               SCOREVector,
        #               AType))
        # }
        
        if ("SCORE" %in% colnames(Pairs)) {
          Res[[m1]][[m2]] <- data.frame("p1" = paste(rep(w1,
                                                         times = L2),
                                                     rep(IMat[[m2]][1L, 2L],
                                                         times = L2),
                                                     p1placeholder,
                                                     sep = "_"),
                                        "p2" = paste(rep(w2,
                                                         times = L2),
                                                     rep(IMat[[m2]][1L, 5L],
                                                         times = L2),
                                                     p2placeholder,
                                                     sep = "_"),
                                        "ExactMatch" = rep(0L,
                                                           times = L2),
                                        "TotalKmers" = rep(0L,
                                                           times = L2),
                                        "MaxKmer" = rep(0L,
                                                        times = L2),
                                        "Consensus" = rep(0,
                                                          times = L2),
                                        "p1FeatureLength" = p1FeatureLength,
                                        "p2FeatureLength" = p2FeatureLength,
                                        "Adjacent" = rep(0L,
                                                         times = L2),
                                        "TetDist" = rep(0,
                                                        times = L2),
                                        "PID" = PIDVector,
                                        "SCORE" = SCOREVector,
                                        "PIDType" = AType,
                                        "PredictedPID" = rep(NA_real_,
                                                             times = L2),
                                        stringsAsFactors = FALSE)
        } else {
          Res[[m1]][[m2]] <- data.frame("p1" = paste(rep(w1,
                                                         times = L2),
                                                     rep(IMat[[m2]][1L, 2L],
                                                         times = L2),
                                                     p1placeholder,
                                                     sep = "_"),
                                        "p2" = paste(rep(w2,
                                                         times = L2),
                                                     rep(IMat[[m2]][1L, 5L],
                                                         times = L2),
                                                     p2placeholder,
                                                     sep = "_"),
                                        "ExactMatch" = rep(0L,
                                                           times = L2),
                                        "TotalKmers" = rep(0L,
                                                           times = L2),
                                        "MaxKmer" = rep(0L,
                                                        times = L2),
                                        "Consensus" = rep(0,
                                                          times = L2),
                                        "p1FeatureLength" = p1FeatureLength,
                                        "p2FeatureLength" = p2FeatureLength,
                                        "Adjacent" = rep(0L,
                                                         times = L2),
                                        "TetDist" = rep(0,
                                                        times = L2),
                                        "PID" = PIDVector,
                                        "PIDType" = AType,
                                        "PredictedPID" = rep(NA_real_,
                                                             times = L2),
                                        stringsAsFactors = FALSE)
        }
        
      } else {
        # do nothing
      }
      
    } # end m2 loop
    
    
  }
  
  Res <- unlist(Res,
                recursive = FALSE)
  Res <- do.call(rbind,
                 Res)
  
  # check for duplicates
  if (nrow(Res) > 1L) {
    IDS <- paste(Res$p1,
                 Res$p2,
                 sep = "_")
    check1 <- !duplicated(IDS)
    Res <- Res[check1, ]
    if (nrow(Res) > 0L) {
      IDS <- paste(Res$p1,
                   Res$p2,
                   sep = "_")
      check2 <- !(IDS %in% POIDs)
      Res <- Res[check2, ]
    }
    
  }
  if (nrow(Res) < 1L) {
    return(NULL)
  }
  
  if (NewPairsOnly) {
    if (Verbose) {
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    return(Res)
  } else {
    if (nrow(Res) > 0L) {
      FeaturesMat2 <- do.call(rbind,
                              strsplit(x = paste(Res$p1,
                                                 Res$p2,
                                                 sep = "_"),
                                       split = "_",
                                       fixed = TRUE))
      FeaturesMat2 <- data.frame("g1" = as.integer(FeaturesMat2[, 1L]),
                                 "i1" = as.integer(FeaturesMat2[, 2L]),
                                 "f1" = as.integer(FeaturesMat2[, 3L]),
                                 "g2" = as.integer(FeaturesMat2[, 4L]),
                                 "i2" = as.integer(FeaturesMat2[, 5L]),
                                 "f2" = as.integer(FeaturesMat2[, 6L]))
      FeaturesMat3 <- rbind(FeaturesMat,
                            FeaturesMat2)
      
      o1 <- order(FeaturesMat3[, 4L],
                  FeaturesMat3[, 1L],
                  FeaturesMat3[, 3L])
      
      Res2 <- rbind(Pairs,
                    Res)
      Res2 <- Res2[o1, ]
      rownames(Res2) <- NULL
      if (Verbose) {
        TimeEnd <- Sys.time()
        print(TimeEnd - TimeStart)
      }
      return(Res2)
    } else {
      if (Verbose) {
        TimeEnd <- Sys.time()
        print(TimeEnd - TimeStart)
      }
      return(Pairs)
    }
  }
  
}


