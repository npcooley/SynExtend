###### -- Expand blocks of predicted pairs ------------------------------------
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

ExpandDiagonal <- function(SynExtendObject,
                           DataBase,
                           InheritConfidence = TRUE,
                           GapTolerance = 100L,
                           DropSingletons = FALSE,
                           UserConfidence = list("PID" = 0.3),
                           Verbose = FALSE) {
  # start with timing
  if (Verbose) {
    pBar <- txtProgressBar(style = 1)
    TimeStart <- Sys.time()
  }
  # overhead checking
  if (!is(object = SynExtendObject,
          class2 = "PairSummaries")) {
    stop ("SynExtendObject must be an object of class 'PairSummaries'.")
  }
  # check DBPATH first
  if (is.character(DataBase)) {
    if (!requireNamespace(package = "RSQLite",
                          quietly = TRUE)) {
      stop("Package 'RSQLite' must be installed.")
    }
    if (!("package:RSQLite" %in% search())) {
      print("Eventually character vector access to DECIPHER DBs will be deprecated.")
      require(RSQLite, quietly = TRUE)
    }
    dbConn <- dbConnect(dbDriver("SQLite"), DataBase)
    on.exit(dbDisconnect(dbConn))
  } else {
    dbConn <- DataBase
    if (!dbIsValid(dbConn)) {
      stop("The connection has expired.")
    }
  }
  # check that user criteria makes sense
  # check the dimensions of the object

  AlignmentFun <- attr(x = SynExtendObject,
                       which = "AlignmentFunction")
  if (!(AlignmentFun %in% c("AlignPairs",
                            "AlignProfiles"))) {
    stop ("Unrecognized alignment function in the 'AlignmentFunction' attribute.")
  }
  if (AlignmentFun == "AlignProfiles") {
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
  }
  if (InheritConfidence) {
    UserConfidence <- attr(x = SynExtendObject,
                           which = "UserConfidence")
    if (is.null(UserConfidence)) {
      stop ("No 'UserConfidence' attribute present in supplied SynExtendObject.")
    }
  } else {
    if (is.null(UserConfidence)) {
      stop ("When InheritConfidence is FALSE, users must supply 'UserConfidence'.")
    } else {
      UserConfidence <- UserConfidence
    }
  }
  Criteria <- names(UserConfidence)
  Floor <- unname(UserConfidence)
  KmerSize <- attr(x = SynExtendObject,
                   which = "KVal")
  # break PairSummaries object down into a workable format
  # build overhead data in a way that makes sense
  GeneCalls <- attr(x = SynExtendObject,
                    which = "GeneCalls")
  # SuppliedFeatures <- as.integer(FeatureSeqs$IDs)
  GCIDs <- as.integer(names(GeneCalls))
  L <- length(GeneCalls)
  L2 <- (L * (L - 1L)) / 2L
  POIDs <- paste(SynExtendObject$p1,
                 SynExtendObject$p2,
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
  
  # change Res to lower case in all places!
  Res <- vector(mode = "list",
                length = nrow(GMat))
  
  # we have to scroll through each occupied cell of the upper triangle, like before
  # now instead of figuring out our edges based on hits, we're seeing whether any
  # 'adjacent' edges along the diagonals fit a user defined criteria for a positive edge
  prev_w1 <- 0L
  prev_w2 <- 0L
  for (m1 in seq(nrow(GMat))) {
    # what is the current comparison
    w1 <- GMat[m1, 1L]
    w2 <- GMat[m1, 2L]
    # wa1 <- match(x = as.character(w1),
    #              table = names(attr(x = Pairs,
    #                                 which = "GeneCalls")))
    # wa2 <- match(x = as.character(w2),
    #              table = names(attr(x = Pairs,
    #                                 which = "GeneCalls")))
    
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
    # i3 <- unique(GeneCalls[[i1]]$Index)
    # i4 <- unique(GeneCalls[[i2]]$Index)
    i3 <- seq_len(max(GeneCalls[[i1]]$Index))
    i4 <- seq_len(max(GeneCalls[[i2]]$Index))
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
                        if (length(x) > 0L) {
                          min(x)
                        } else {
                          NA # NA Placeholder - should never be accessed, but necessary to keep index matching
                        }
                      },
                      simplify = TRUE)
    i1upper <- sapply(i5,
                      function(x) {
                        if (length(x) > 0L) {
                          max(x)
                        } else {
                          NA
                        }
                      },
                      simplify = TRUE)
    i2lower <- sapply(i6,
                      function(x) {
                        if (length(x) > 0L) {
                          min(x)
                        } else {
                          NA # NA Placeholder - should never be accessed, but necessary to keep index matching
                        }
                      },
                      simplify = TRUE)
    i2upper <- sapply(i6,
                      function(x) {
                        if (length(x) > 0L) {
                          max(x)
                        } else {
                          NA
                        }
                      },
                      simplify = TRUE)
    
    # pull seqs, this needs to be more efficient in the future
    if (prev_w1 != w1) {
      out1 <- SearchDB(dbFile = dbConn,
                       tblName = "NTs",
                       identifier = as.character(w1),
                       verbose = FALSE,
                       nameBy = "description")
      out2 <- SearchDB(dbFile = dbConn,
                       tblName = "AAs",
                       identifier = as.character(w1),
                       verbose = FALSE,
                       nameBy = "description")
      
      # DBQUERY <- paste("select length, mod, code, cds from NTs where identifier is",
      #                  w1)
      # out3 <- dbGetQuery(conn = dbConn,
      #                    statement = DBQUERY)
      CurrentW1Seqs <- list("DNA" = out1,
                            "AA" = out2,
                            "Struct" = PredictHEC(myAAStringSet = out2,
                                                  type = "probabilities",
                                                  HEC_MI1 = MAT1,
                                                  HEC_MI2 = MAT2),
                            "NTCount" = width(out1),
                            "CodingVal1" = GeneCalls[[match(x = w1,
                                                            table = GCIDs)]]$Coding,
                            "CodingVal2" = width(out1) %% 3L == 0L,
                            "CDSCount" = lengths(GeneCalls[[match(x = w1,
                                                                  table = GCIDs)]]$Range))
    }
    
    if (prev_w2 != w2) {
      out1 <- SearchDB(dbFile = dbConn,
                       tblName = "NTs",
                       identifier = as.character(w2),
                       verbose = FALSE,
                       nameBy = "description")
      out2 <- SearchDB(dbFile = dbConn,
                       tblName = "AAs",
                       identifier = as.character(w2),
                       verbose = FALSE,
                       nameBy = "description")
      # DBQUERY <- paste("select length, mod, code, cds from NTs where identifier is",
      #                  w2)
      # out3 <- dbGetQuery(conn = dbConn,
      #                    statement = DBQUERY)
      CurrentW2Seqs <- list("DNA" = out1,
                            "AA" = out2,
                            "Struct" = PredictHEC(myAAStringSet = out2,
                                                  type = "probabilities",
                                                  HEC_MI1 = MAT1,
                                                  HEC_MI2 = MAT2),
                            "NTCount" = width(out1),
                            "CodingVal1" = GeneCalls[[match(x = w2,
                                                            table = GCIDs)]]$Coding,
                            "CodingVal2" = width(out1) %% 3L == 0L,
                            "CDSCount" = lengths(GeneCalls[[match(x = w2,
                                                                  table = GCIDs)]]$Range))
    }
    
    
    nuc1 <- oligonucleotideFrequency(x = CurrentW1Seqs$DNA,
                                     width = KmerSize,
                                     as.prob = TRUE)
    nuc2 <- oligonucleotideFrequency(x = CurrentW2Seqs$DNA,
                                     width = KmerSize,
                                     as.prob = TRUE)
    Features01Match <- match(x = names(CurrentW1Seqs$DNA),
                             table = names(CurrentW1Seqs$AA))
    Features02Match <- match(x = names(CurrentW2Seqs$DNA),
                             table = names(CurrentW2Seqs$AA))
    
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
      } # end m3 loop through dr3
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
      } # end m3 loop through dr4
      
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
      
      if (AlignmentFun == "AlignProfiles") {
        
        p1placeholder <- p2placeholder <- p1FeatureLength <- p2FeatureLength <- rep(NA_integer_,
                                                                                    times = VSize)
        PIDVector <- SCOREVector <- NucDist <- rep(NA_real_,
                                                   times = VSize)
        AType <- rep(NA_character_,
                     times = VSize)
        
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
            if (CurrentW1Seqs$CodingVal1[f1 - ci1lower + 1L] &
                CurrentW1Seqs$CodingVal2[f1 - ci1lower + 1L] &
                CurrentW2Seqs$CodingVal1[f2 - ci2lower + 1L] &
                CurrentW2Seqs$CodingVal2[f2 - ci2lower + 1L]) {
              # both are coding
              
              ali <- AlignProfiles(pattern = CurrentW1Seqs$AA[Features01Match[f1 - ci1lower + 1L]],
                                   subject = CurrentW2Seqs$AA[Features02Match[f2 - ci2lower + 1L]],
                                   p.struct = CurrentW1Seqs$Struct[Features01Match[f1 - ci1lower + 1L]],
                                   s.struct = CurrentW2Seqs$Struct[Features02Match[f2 - ci2lower + 1L]])
              PID <- 1 - DistanceMatrix(myXStringSet = ali,
                                        type = "matrix",
                                        includeTerminalGaps = TRUE,
                                        verbose = FALSE)[1L, 2L]
              SCORE <- ScoreAlignment(myXStringSet = ali,
                                      structures = PredictHEC(ali,
                                                              type = "probabilities",
                                                              HEC_MI1 = MAT1,
                                                              HEC_MI2 = MAT2),
                                      structureMatrix = structureMatrix)
              
              CType <- "AA"
            } else {
              # at least one is not coding
              ali <- AlignProfiles(pattern = CurrentW1Seqs$DNA[f1 - ci1lower + 1L],
                                   subject = CurrentW2Seqs$DNA[f2 - ci2lower + 1L])
              PID <- 1 - DistanceMatrix(myXStringSet = ali,
                                        type = "matrix",
                                        includeTerminalGaps = TRUE,
                                        verbose = FALSE)[1L, 2L]
              SCORE <- ScoreAlignment(myXStringSet = ali,
                                      substitutionMatrix = substitutionMatrix)
              CType <- "NT"
            } # end if else statement for coding / non
            CurrentDist <- sqrt(sum((nuc1[f1 - ci1lower + 1L, ] - nuc2[f2 - ci2lower + 1L, ])^2)) / ((sum(nuc1[f1 - ci1lower + 1L, ]) + sum(nuc2[f2 - ci2lower + 1L, ])) / 2)
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
              p1FeatureLength[Count] <- CurrentW1Seqs$NTCount[f1 - ci1lower + 1L]
              p2FeatureLength[Count] <- CurrentW2Seqs$NTCount[f2 - ci2lower + 1L]
              PIDVector[Count] <- PID
              SCOREVector[Count] <- SCORE
              AType[Count] <- CType
              NucDist[Count] <- CurrentDist
              
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
                SCOREVector <- c(SCOREVector,
                                 rep(NA_real_,
                                     times = VSize))
                AType <- c(AType,
                           rep(NA_character_,
                               times = VSize))
                NucDist <- c(NucDist,
                             rep(NA_real_,
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
          close(pBar)
        }
        
        if (any(!is.na(PIDVector))) {
          L2 <- max(which(!is.na(PIDVector)))
          
          p1placeholder <- p1placeholder[seq_len(L2)]
          p2placeholder <- p2placeholder[seq_len(L2)]
          p1FeatureLength <- p1FeatureLength[seq_len(L2)]
          p2FeatureLength <- p2FeatureLength[seq_len(L2)]
          PIDVector <- PIDVector[seq_len(L2)]
          SCOREVector <- SCOREVector[seq_len(L2)]
          AType <- AType[seq_len(L2)]
          NucDist <- NucDist[seq_len(L2)]
          if (is.null(attr(x = SynExtendObject,
                           which = "Retain"))) {
            NewClusterID <- -1L
          } else {
            NewClusterID <- max(as.integer(names(attr(x = SynExtendObject,
                                                      which = "Retain")))) + 1L
          }
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
                                        "Consensus" = rep(0,
                                                          times = L2),
                                        "p1featurelength" = p1FeatureLength,
                                        "p2featurelength" = p2FeatureLength,
                                        "blocksize" = rep(1L,
                                                          times = L2),
                                        "KDist" = NucDist,
                                        "TotalMatch" = rep(0L,
                                                           times = L2),
                                        "MaxMatch" = rep(0L,
                                                         times = L2),
                                        "UniqueMatches" = rep(0L,
                                                              times = L2),
                                        "PID" = PIDVector,
                                        "Score" = SCOREVector,
                                        "Alignment" = AType,
                                        "Block_UID" = rep(-1L,
                                                          times = L2),
                                        "ClusterID" = rep(NewClusterID,
                                                          times = L2),
                                        stringsAsFactors = FALSE)
          
        } else {
          # do nothing
        }
        
      } else if (AlignmentFun == "AlignPairs") {
        stop ("Expansion with AlignPairs() is not yet implemented.")
      }
      
    } # end m2 loop
  } # end m1 loop
  
  Res <- unlist(Res,
                recursive = FALSE)
  Res <- do.call(rbind,
                 Res)
  
  if (is.null(Res)) {
    # no new pairs discovered
    Res <- data.frame("p1" = character(0L),
                      "p2" = character(0L),
                      "Consensus" = numeric(0L),
                      "p1featurelength" = integer(0L),
                      "p2featurelength" = integer(0L),
                      "blocksize" = integer(0L),
                      "KDist" = integer(0L),
                      "TotalMatch" = integer(0L),
                      "MaxMatch" = integer(0L),
                      "UniqueMatches" = integer(0L),
                      "PID" = numeric(0L),
                      "SCORE" = numeric(0L),
                      "Alignment" = character(0L),
                      "Block_UID" = integer(0L),
                      "ClusterID" = integer(0L),
                      stringsAsFactors = FALSE)
  } else {
    # if Res is not NULL
    # remove duplicates
    if (nrow(Res) > 1L) {
      IDS <- paste(Res$p1,
                   Res$p2,
                   sep = "_")
      check1 <- !duplicated(IDS)
      Res <- Res[check1, , drop = FALSE]
      # remove IDs present in the original object IDs
      if (nrow(Res) > 0L) {
        IDS <- paste(Res$p1,
                     Res$p2,
                     sep = "_")
        check2 <- !(IDS %in% POIDs)
        Res <- Res[check2, , drop = FALSE]
      }
    }
  } # end logical check for whether Res is NULL or not
  # if an object that wasn't clustered beforehand was supplied, delete the column
  if (!("ClusterID" %in% colnames(SynExtendObject))) {
    Res <- Res[, -which(colnames(Res) == "ClusterID")]
  }
  # if new rows have been added, fold them in, recalculate blocksize
  if (nrow(Res) > 0L) {
    # return(list(SynExtendObject,
    #             Res,
    #             CurrentW1Seqs,
    #             CurrentW2Seqs))
    Res <- rbind(SynExtendObject,
                 Res)
    
    block_uid <- 1L
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
    o1 <- order(FeaturesMat2[, 4L],
                FeaturesMat2[, 1L],
                FeaturesMat2[, 5L],
                FeaturesMat2[, 2L],
                FeaturesMat2[, 6L],
                FeaturesMat2[, 3L],
                decreasing = FALSE)
    
    Res <- Res[o1, ]
    FeaturesMat2 <- FeaturesMat2[o1, ]
    rownames(Res) <- NULL
    # reset the blocksizes
    blockreplacement <- rep(1L,
                            nrow(Res))
    for (m2 in seq_len(nrow(GMat))) {
      wi1 <- which(FeaturesMat2[, 1L] == GMat[m2, 1L] &
                     FeaturesMat2[, 4L] == GMat[m2, 2L])
      names(wi1) <- rownames(FeaturesMat2[wi1, ])
      
      # block size determination
      if (length(wi1) > 1) {
        SubFeatures <- FeaturesMat2[wi1, c(2,3,5,6)]
        # only run block size checks if enough rows are present
        # FeaturesMat <- data.frame("i1" = IMatrix[, 1L],
        #                           "f1" = PMatrix[, 1L],
        #                           "i2" = IMatrix[, 2L],
        #                           "f2" = PMatrix[, 2L])
        dr1 <- SubFeatures[, 2L] + SubFeatures[, 4L]
        dr2 <- SubFeatures[, 2L] - SubFeatures[, 4L]
        InitialBlocks1 <- unname(split(x = SubFeatures,
                                       f = list(as.integer(SubFeatures[, 1L]),
                                                as.integer(SubFeatures[, 3L]),
                                                dr1),
                                       drop = TRUE))
        InitialBlocks2 <- unname(split(x = SubFeatures,
                                       f = list(as.integer(SubFeatures[, 1L]),
                                                as.integer(SubFeatures[, 3L]),
                                                dr2),
                                       drop = TRUE))
        Blocks <- c(InitialBlocks1[sapply(InitialBlocks1,
                                          function(x) nrow(x),
                                          simplify = TRUE) > 1],
                    InitialBlocks2[sapply(InitialBlocks2,
                                          function(x) nrow(x),
                                          simplify = TRUE) > 1])
        L01 <- length(Blocks)
        if (L01 > 0) {
          for (m3 in seq_along(Blocks)) {
            # blocks are guaranteed to contain more than 1 row
            
            sp1 <- vector(mode = "integer",
                          length = nrow(Blocks[[m3]]))
            # we need to check both columns here, this currently is not correct
            # in all cases
            sp2 <- Blocks[[m3]][, 4L]
            sp3 <- Blocks[[m3]][, 2L]
            
            it1 <- 1L
            it2 <- sp2[1L]
            it4 <- sp3[1L]
            # create a map vector on which to split the groups, if necessary
            for (m4 in seq_along(sp1)) {
              it3 <- sp2[m4]
              it5 <- sp3[m4]
              if ((it3 - it2 > 1L) |
                  (it5 - it4 > 1L)) {
                # if predicted pairs are not contiguous, update the iterator
                it1 <- it1 + 1L
              }
              sp1[m4] <- it1
              it2 <- it3
              it4 <- it5
            }
            
            # if the splitting iterator was updated at all, a gap was detected
            if (it1 > 1L) {
              Blocks[[m3]] <- unname(split(x = Blocks[[m3]],
                                           f = sp1))
            } else {
              Blocks[[m3]] <- Blocks[m3]
            }
            
          } # end m3 loop
          # Blocks is now a list where each position is a set of blocked pairs
          Blocks <- unlist(Blocks,
                           recursive = FALSE)
          # drop blocks of size 1, they do not need to be evaluated
          Blocks <- Blocks[sapply(X = Blocks,
                                  FUN = function(x) {
                                    nrow(x)
                                  },
                                  simplify = TRUE) > 1L]
          L01 <- length(Blocks)
          AbsBlockSize <- rep(1L,
                              nrow(SubFeatures))
          BlockID_Map <- rep(-1L,
                             nrow(SubFeatures))
          # only bother with this if there are blocks remaining
          # otherwise AbsBlockSize, which is initialized as a vector of 1s
          # will be left as a vector of 1s, all pairs are singleton pairs in this scenario
          if (L01 > 0L) {
            for (m3 in seq_along(Blocks)) {
              # rownames of the Blocks dfs relate to row positions in the original
              # matrix
              pos <- match(x = rownames(Blocks[[m3]]),
                           table = names(wi1))
              val <- rep(nrow(Blocks[[m3]]),
                         nrow(Blocks[[m3]]))
              # do not overwrite positions that are in larger blocks
              keep <- AbsBlockSize[pos] < val
              if (any(keep)) {
                AbsBlockSize[pos[keep]] <- val[keep]
                BlockID_Map[pos[keep]] <- rep(block_uid,
                                              sum(keep))
                block_uid <- block_uid + 1L
              }
            } # end m3 loop
          } # end logical check for block size
        } else {
          # no blocks observed, all pairs present are singleton pairs
          AbsBlockSize <- rep(1L,
                              nrow(SubFeatures))
          BlockID_Map <- rep(-1L,
                             nrow(SubFeatures))
        }
      } else {
        AbsBlockSize <- 1L
      }
      blockreplacement[wi1] <- AbsBlockSize
      
    } # end m2 loop through genome comparisons
    Res$blocksize <- blockreplacement
    attr(x = Res,
         which = "KVal") <- attr(x = SynExtendObject,
                                 which = "KVal")
    if (Verbose) {
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    return(Res)
  } else {
    if (Verbose) {
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    return(SynExtendObject)
  }
}

