###### -- summarize seqs ------------------------------------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu
# given a linked pairs object, and a FeatureSeqs object return a pairsummaries
# object
# this function will always align, unlike PairSummaries
# TODO:
# implement ShowPlot Processors and elipses
# It's not necessarily good to store the structure calls that align profiles needs
# so we're going to generate them on the fly
# we're switching the DataBase that we're pulling from because we shouldn't need
# the raw genome

SummarizePairs <- function(SynExtendObject,
                           DataBase01,
                           AlignmentFun = "AlignProfiles",
                           RetainAnchors = FALSE,
                           DefaultTranslationTable = "11",
                           KmerSize = 5,
                           IgnoreDefaultStringSet = FALSE,
                           Verbose = FALSE,
                           ShowPlot = FALSE,
                           Processors = 1,
                           Storage = 2,
                           ...) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  
  # overhead checking
  # object types
  if (!is(object = SynExtendObject,
          class2 = "LinkedPairs")) {
    stop ("'SynExtendObject' is not an object of class 'LinkedPairs'.")
  }
  Size <- nrow(SynExtendObject)
  # if (!is(object = FeatureSeqs,
  #         class2 = "FeatureSeqs")) {
  #   stop ("'FeatureSeqs' is not an object of class 'FeatureSeqs'.")
  # }
  # we should only need to talk to the DataBase IF FeatureSeqs is not the right length
  
  if (is.character(DataBase01)) {
    if (!requireNamespace(package = "RSQLite",
                          quietly = TRUE)) {
      stop("Package 'RSQLite' must be installed.")
    }
    if (!("package:RSQLite" %in% search())) {
      print("Eventually character vector access to DECIPHER DBs will be deprecated.")
      require(RSQLite, quietly = TRUE)
    }
    dbConn <- dbConnect(dbDriver("SQLite"), DataBase01)
    on.exit(dbDisconnect(dbConn))
  } else {
    dbConn <- DataBase01
    if (!dbIsValid(dbConn)) {
      stop("The connection has expired.")
    }
  }
  if (!is.character(AlignmentFun) |
      length(AlignmentFun) > 1) {
    stop("AlignmentFun must be either 'AlignPairs' or 'AlignProfiles'.")
  }
  if (!(AlignmentFun %in% c("AlignPairs",
                            "AlignProfiles"))) {
    stop("AlignmentFun must be either 'AlignPairs' or 'AlignProfiles'.")
  }
  if (!is.character(DefaultTranslationTable) |
      length(DefaultTranslationTable) > 1) {
    stop("DefaultTranslationTable must be a character of length 1.")
  }
  # check storage
  if (Storage < 0) {
    stop("Storage must be greater than zero.")
  } else {
    Storage <- Storage * 1e9 # conversion to gigabytes
  }
  # deal with Processors, this mimics Erik's error checking
  if (!is.null(Processors) && !is.numeric(Processors)) {
    stop("Processors must be a numeric.")
  }
  if (!is.null(Processors) && floor(Processors) != Processors) {
    stop("Processors must be a whole number.")
  }
  if (!is.null(Processors) && Processors < 1) {
    stop("Processors must be at least 1.")
  }
  if (is.null(Processors)) {
    Processors <- DECIPHER:::.detectCores()
  } else {
    Processors <- as.integer(Processors)
  }
  # deal with user arguments
  # ignore 'verbose'
  UserArgs <- list(...)
  # if (length(UserArgs) > 0) {
  #   UserArgNames <- names(UserArgs)
  #   AlignPairsArgs <- formals(AlignPairs)
  #   AlignPairsArgNames <- names(AlignPairsArgs)
  #   AlignProfilesArgs <- formals(AlignProfiles)
  #   AlignProfilesArgNames <- names(AlignProfilesArgs)
  #   DistanceMatrixArgs <- formals(DistanceMatrix)
  #   DistanceMatrixArgNames <- names(DistanceMatrixArgs)
  #   
  #   APaArgNames <- APrArgNames <- DMArgNames <- vector(mode = "character",
  #                                                    length = length(UserArgNames))
  #   for (m1 in seq_along(UserArgNames)) {
  #     APaArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                  choices = AlignPairsArgNames)
  #     APrArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                  choices = AlignProfilesArgNames)
  #     DMArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                 choices = DistanceMatrixArgNames)
  #   }
  #   # set the user args for AlignProfiles
  #   if (any(APaArgNames)) {
  #     
  #   }
  #   # set the user args for AlignPairs
  #   # set the user args for DistanceMatrix
  # }
  
  GeneCalls <- attr(x = SynExtendObject,
                    which = "GeneCalls")
  GeneCallIDs <- names(GeneCalls)
  ObjectIDs <- rownames(SynExtendObject)
  # when we subset a LinkedPairsObject it doesn't smartly handle the genecalls yet...
  if (!all(ObjectIDs %in% GeneCallIDs)) {
    stop("Function expects all IDs in the SynExtendObject to supplied GeneCalls.")
  }
  IDMatch <- match(x = ObjectIDs,
                   table = GeneCallIDs)
  GeneCalls <- GeneCalls[IDMatch]
  # we're only going to scroll through the cells that are supplied in the object
  DataPool <- vector(mode = "list",
                     length = length(ObjectIDs))
  
  MAT1 <- get(data("HEC_MI1",
                   package = "DECIPHER",
                   envir = environment()))
  MAT2 <- get(data("HEC_MI2",
                   package = "DECIPHER",
                   envir = environment()))
  # set initial progress bars and iterators
  # PH is going to be the container that 'res' eventually gets constructed from
  # while Total
  if (AlignmentFun == "AlignProfiles") {
    # progress bar ticks through each alignment
    # this is the slower non-default option
    Total <- (Size - (Size - 1L)) / 2
    PH <- vector(mode = "list",
                 length = Total)
    Total <- sum(sapply(SynExtendObject[upper.tri(SynExtendObject)],
                        function(x) nrow(x),
                        USE.NAMES = FALSE,
                        simplify = TRUE))
    
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
    
    # read this message after neighbors and k-mer dists ?
    if (Verbose) {
      cat("Aligning pairs.\n")
    }
  } else if (AlignmentFun == "AlignPairs") {
    # progress bar ticks through each cell
    # this is the faster default option
    # technically these alignments have the possibility of not being
    # as good as align profiles, but that's a hit we're willing to take
    Total <- (Size * (Size - 1L)) / 2L
    PH <- vector(mode = "list",
                 length = Total)
    Total <- Total * 2L
    if (Verbose) {
      cat("Collecting pairs.\n")
    }
  }
  # for testing:
  # NT_Anchors <- AA_Anchors <- vector(mode = "list",
  #                                    length = length(PH))
  
  Count <- 1L
  PBCount <- 0L
  # upper key!
  # QueryGene == 1
  # SubjectGene == 2
  # ExactOverlap == 3
  # QueryIndex == 4
  # SubjectIndex == 5
  # QLeft == 6
  # QRight == 7
  # SLeft == 8
  # SRight == 9
  # MaxKmer == 10
  # TotalKmer == 11
  
  # lower key!
  # same minus max and total, but for individual linking kmers
  # not all linking kmers
  block_uid <- 1L
  Prev_m1 <- 0L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      ###### -- only evaluate valid positions ---------------------------------
      if (nrow(SynExtendObject[[m1, m2]]) > 0L) {
        # links table is populated, do whatever
        PMatrix <- cbind(SynExtendObject[[m1, m2]][, 1L],
                         SynExtendObject[[m1, m2]][, 2L])
        
        IMatrix <- cbind(SynExtendObject[[m1, m2]][, 4L],
                         SynExtendObject[[m1, m2]][, 5L])
        # find the Consensus of each linking kmer hit
        # this reference to GeneCalls is still fine
        p1l <- GeneCalls[[m1]]$Stop[SynExtendObject[[m2, m1]][, 1L]] - GeneCalls[[m1]]$Start[SynExtendObject[[m2, m1]][, 1L]] + 1L
        p2l <- GeneCalls[[m2]]$Stop[SynExtendObject[[m2, m1]][, 2L]] - GeneCalls[[m2]]$Start[SynExtendObject[[m2, m1]][, 2L]] + 1L
        a1 <- SynExtendObject[[m2, m1]][, 6L]
        a2 <- SynExtendObject[[m2, m1]][, 7L]
        b1 <- SynExtendObject[[m2, m1]][, 8L]
        b2 <- SynExtendObject[[m2, m1]][, 9L]
        c1 <- GeneCalls[[m1]]$Start[SynExtendObject[[m2, m1]][, 1L]]
        c2 <- GeneCalls[[m1]]$Stop[SynExtendObject[[m2, m1]][, 1L]]
        d1 <- GeneCalls[[m2]]$Start[SynExtendObject[[m2, m1]][, 2L]]
        d2 <- GeneCalls[[m2]]$Stop[SynExtendObject[[m2, m1]][, 2L]]
        s1 <- GeneCalls[[m1]]$Strand[SynExtendObject[[m2, m1]][, 1L]] == 0L
        s2 <- GeneCalls[[m2]]$Strand[SynExtendObject[[m2, m1]][, 2L]] == 0L
        
        diff1 <- mapply(function(o, p, q, r, s, t, u, v, w, x, y, z) {
          if (o == p) {
            mean(c(abs((abs(w - s) / q) - (abs(y - u) / r)),
                   abs((abs(t - x) / q) - (abs(v - z) / r))))
          } else {
            mean(c(abs((abs(w - s) / q) - (abs(v - z) / r)),
                   abs((abs(t - x) / q) - (abs(y - u) / r))))
          }
        },
        o = s1,
        p = s2,
        q = p1l,
        r = p2l,
        s = a1,
        t = a2,
        u = b1,
        v = b2,
        w = c1,
        x = c2,
        y = d1,
        z = d2)
        
        # get the mean consensus
        diff2 <- vector(mode = "numeric",
                        length = nrow(SynExtendObject[[m1, m2]]))
        for (m3 in seq_along(diff2)) {
          diff2[m3] <- 1 - mean(diff1[SynExtendObject[[m2, m1]][, 1L] == SynExtendObject[[m1, m2]][m3, 1L] &
                                        SynExtendObject[[m2, m1]][, 2L] == SynExtendObject[[m1, m2]][m3, 2L]])
        }
        # max match size
        MatchMax <- SynExtendObject[[m1, m2]][, "MaxKmerSize"]
        # total unique matches
        UniqueMatches <- SynExtendObject[[m1, m2]][, "TotalKmerHits"]
        # total matches at all
        TotalMatch <- SynExtendObject[[m1, m2]][, "ExactOverlap"]
        
        # block size determination
        if (nrow(SynExtendObject[[m1, m2]]) > 1) {
          # only run block size checks if enough rows are present
          FeaturesMat <- data.frame("i1" = IMatrix[, 1L],
                                    "f1" = PMatrix[, 1L],
                                    "i2" = IMatrix[, 2L],
                                    "f2" = PMatrix[, 2L])
          dr1 <- FeaturesMat[, 2L] + FeaturesMat[, 4L]
          dr2 <- FeaturesMat[, 2L] - FeaturesMat[, 4L]
          InitialBlocks1 <- unname(split(x = FeaturesMat,
                                         f = list(as.integer(FeaturesMat[, 1L]),
                                                  as.integer(FeaturesMat[, 3L]),
                                                  dr1),
                                         drop = TRUE))
          InitialBlocks2 <- unname(split(x = FeaturesMat,
                                         f = list(as.integer(FeaturesMat[, 1L]),
                                                  as.integer(FeaturesMat[, 3L]),
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
                                nrow(FeaturesMat))
            BlockID_Map <- rep(-1L,
                               nrow(FeaturesMat))
            # only bother with this if there are blocks remaining
            # otherwise AbsBlockSize, which is initialized as a vector of 1s
            # will be left as a vector of 1s, all pairs are singleton pairs in this scenario
            if (L01 > 0L) {
              for (m3 in seq_along(Blocks)) {
                # rownames of the Blocks dfs relate to row positions in the original
                # matrix
                pos <- as.integer(rownames(Blocks[[m3]]))
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
                                nrow(FeaturesMat))
            BlockID_Map <- rep(-1L,
                               nrow(FeaturesMat))
          }
        } else {
          AbsBlockSize <- 1L
          BlockID_Map <- -1L
        }
        # BlockSize evaluation is complete
        
        # from here we need to get the kmer differences
        # the PIDs
        # the SCOREs
        
        # if both positions are present in the FeatureSeqs object, do nothing
        # if either or both is missing, they need to be pulled from the DB
        # this functionality needs to be expanded eventually to take in cases where the
        # we're overlaying something with gene calls on something that doesn't have gene calls
        # if (!all(ObjectIDs[c(m1, m2)] %in% FeatureSeqs$IDs)) {
        #   # an object ID does not have a seqs present, pull them
        #   # this is not a priority so we're leaving this blank for a second
        # } else {
        #   TMPSeqs01 <- FALSE
        #   TMPSeqs02 <- FALSE
        # }
        
        # regardless of how we align, we need to prepare and collect
        # the same basic statistics and look aheads
        # this needs to change a little,
        # I'll be pulling the seqs from a DB, and the integers and things as well
        if (is.null(DataPool[[m1]])) {
          # the pool position is empty, pull from the DB
          # and generate the AAStructures
          DataPool[[m1]]$DNA <- SearchDB(dbFile = dbConn,
                                         tblName = "NTs",
                                         identifier = ObjectIDs[m1],
                                         verbose = FALSE,
                                         nameBy = "description")
          DataPool[[m1]]$AA <- SearchDB(dbFile = dbConn,
                                        tblName = "AAs",
                                        identifier = ObjectIDs[m1],
                                        verbose = FALSE,
                                        nameBy = "description")
          DataPool[[m1]]$len <- width(DataPool[[m1]]$DNA)
          DataPool[[m1]]$mod <- DataPool[[m1]]$len %% 3L == 0
          DataPool[[m1]]$code <- GeneCalls[[m1]]$Coding
          DataPool[[m1]]$cds <- lengths(GeneCalls[[m1]]$Range)
          # DBQUERY <- paste("select len, mod, code, cds from NTs where identifier is",
          #                  ObjectIDs[m1])
          # DBOUT <- dbGetQuery(conn = dbConn,
          #                     statement = DBQUERY)
          # DataPool[[m1]]$len <- DBOUT$len
          # DataPool[[m1]]$mod <- DBOUT$mod
          # DataPool[[m1]]$code <- DBOUT$code
          # DataPool[[m1]]$cds <- DBOUT$cds
          DataPool[[m1]]$struct <- PredictHEC(myAAStringSet = DataPool[[m1]]$AA,
                                              type = "probabilities",
                                              HEC_MI1 = MAT1,
                                              HEC_MI2 = MAT2)
          
        } else {
          # the pool position is not empty, assume that it's populated with all the information
          # that it needs
          # DataPool[[m1]]$struct <- PredictHEC(myAAStringSet = DataPool[[m1]]$QueryAA,
          #                                     type = "probabilities",
          #                                     HEC_MI1 = MAT1,
          #                                     HEC_MI2 = MAT2)
        }
        
        if (is.null(DataPool[[m2]])) {
          # the pool position is empty, pull from DB
          # and generate the AAStructures
          DataPool[[m2]]$DNA <- SearchDB(dbFile = dbConn,
                                         tblName = "NTs",
                                         identifier = ObjectIDs[m2],
                                         verbose = FALSE,
                                         nameBy = "description")
          DataPool[[m2]]$AA <- SearchDB(dbFile = dbConn,
                                        tblName = "AAs",
                                        identifier = ObjectIDs[m2],
                                        verbose = FALSE,
                                        nameBy = "description")
          DataPool[[m2]]$len <- width(DataPool[[m2]]$DNA)
          DataPool[[m2]]$mod <- DataPool[[m2]]$len %% 3L == 0
          DataPool[[m2]]$code <- GeneCalls[[m2]]$Coding
          DataPool[[m2]]$cds <- lengths(GeneCalls[[m2]]$Range)
          # DBQUERY <- paste("select len, mod, code, cds from NTs where identifier is",
          #                  ObjectIDs[m2])
          # DBOUT <- dbGetQuery(conn = dbConn,
          #                     statement = DBQUERY)
          # DataPool[[m2]]$len <- DBOUT$len
          # DataPool[[m2]]$mod <- DBOUT$mod
          # DataPool[[m2]]$code <- DBOUT$code
          # DataPool[[m2]]$cds <- DBOUT$cds
          DataPool[[m2]]$struct <- PredictHEC(myAAStringSet = DataPool[[m2]]$AA,
                                              type = "probabilities",
                                              HEC_MI1 = MAT1,
                                              HEC_MI2 = MAT2)
        } else {
          # the pool position is not empty, assume that it's populated with all the information
          # that it needs
          # DataPool[[m2]]$struct <- PredictHEC(myAAStringSet = DataPool[[m2]]$QueryAA,
          #                                     type = "probabilities",
          #                                     HEC_MI1 = MAT1,
          #                                     HEC_MI2 = MAT2)
        }
        
        if (Prev_m1 != m1) {
          QueryDNA <- DataPool[[m1]]$DNA
          QueryAA <- DataPool[[m1]]$AA
          QNTCount <- DataPool[[m1]]$len
          QMod <- DataPool[[m1]]$mod
          QCode <- DataPool[[m1]]$code
          QCDSCount <- DataPool[[m1]]$cds
          QueryStruct <- DataPool[[m1]]$struct
        } else {
          # do something else?
        }
        
        SubjectDNA <- DataPool[[m2]]$DNA
        SubjectAA <- DataPool[[m2]]$AA
        SNTCount <- DataPool[[m2]]$len
        SMod <- DataPool[[m2]]$mod
        SCode <- DataPool[[m2]]$code
        SCDSCount <- DataPool[[m2]]$cds
        SubjectStruct <- DataPool[[m2]]$struct
        
        # align everyone as AAs who can be, i.e. modulo of 3, is coding, etc
        # then align everyone else as nucs
        # translate the hit locations to the anchor positions
        # every hit is an anchor
        # all hits are stored in the LinkedPairs object in nucleotide space
        # as left/right bounds, orientations will be dependant upon strandedness of the genes
        
        # prepare the kmer distance stuff:
        NucDist <- vector(mode = "numeric",
                          length = nrow(PMatrix))
        if (KmerSize < 10L) {
          nuc1 <- oligonucleotideFrequency(x = QueryDNA[PMatrix[, 1L]],
                                           width = KmerSize,
                                           as.prob = TRUE)
          nuc2 <- oligonucleotideFrequency(x = SubjectDNA[PMatrix[, 2L]],
                                           width = KmerSize,
                                           as.prob = TRUE)
        } else {
          stop("non-overlapping kmers not implemented yet")
        }
        QueryFeatureLength <- QNTCount[PMatrix[, 1L]]
        SubjectFeatureLength <- SNTCount[PMatrix[, 2L]]
        
        # Perform the alignments
        if (AlignmentFun == "AlignPairs") {
          
          # grab kmer distances ahead of time because AlignPairs doesn't need to loop
          # through anything
          for (m3 in seq_along(NucDist)) {
            NucDist[m3] <- sqrt(sum((nuc1[m3, ] - nuc2[m3, ])^2)) / ((sum(nuc1[m3, ]) + sum(nuc2[m3, ])) / 2)
          }
          
          # we need to check and default stringset logical and the retain anchor logical
          if (IgnoreDefaultStringSet) {
            # spit out something that makes it clear there's only a single call to
            # AlignPairs that is calling NTs
            # not much to do here, set the planning df as PMatrix and call it a day
            # df_plan <- data.frame("Pattern" = match(x = AASubSet[, 1L],
            #                                         table = aa_match1),
            #                       "Subject" = match(x = AASubSet[, 2L],
            #                                         table = aa_match2))
            df_nt <- data.frame("Pattern" = PMatrix[, 1L],
                                "Subject" = PMatrix[, 2L])
            df_aa <- data.frame("Pattern" = integer(),
                                "Subject" = integer())
            AASelect <- rep(FALSE, nrow(PMatrix))
            NTSelect <- rep(TRUE, nrow(PMatrix))
            AASubSet <- PMatrix[AASelect, , drop = FALSE]
            
          } else {
            # spit out the subset vectors and logicals to correctly call both AlignPairs calls
            # and both dfs
            AASelect <- QMod[PMatrix[, 1L]] & QCode[PMatrix[, 1L]] & SMod[PMatrix[, 2L]] & SCode[PMatrix[, 2L]]
            NTSelect <- !AASelect
            AASubSet <- PMatrix[AASelect, , drop = FALSE]
            
            aa_match1 <- strsplit(x = names(QueryAA),
                                  split = "_",
                                  fixed = TRUE)
            aa_match1 <- as.integer(sapply(X = aa_match1,
                                           FUN = function(x) {
                                             x[3]
                                           }))
            aa_match2 <- strsplit(x = names(SubjectAA),
                                  split = "_",
                                  fixed = TRUE)
            aa_match2 <- as.integer(sapply(X = aa_match2,
                                           FUN = function(x) {
                                             x[3]
                                           }))
            df_aa <- data.frame("Pattern" = match(x = AASubSet[, 1L],
                                                  table = aa_match1),
                                "Subject" = match(x = AASubSet[, 2L],
                                                  table = aa_match2))
            df_nt <- data.frame("Pattern" = PMatrix[NTSelect, 1L],
                                "Subject" = PMatrix[NTSelect, 2L])
          }
          if (RetainAnchors) {
            # set the anchors
            # start with nucleotide positions 
            hitsets <- SynExtendObject[[m2, m1]][, seq_len(2L), drop = FALSE]
            WithinQueryNucs <- mapply(USE.NAMES = FALSE,
                                      SIMPLIFY = FALSE,
                                      FUN = function(i1,
                                                     i2,
                                                     left1,
                                                     left2,
                                                     right1,
                                                     right2,
                                                     strand1,
                                                     strand2,
                                                     querygeneboundl,
                                                     querygeneboundr,
                                                     subgeneboundl,
                                                     subgeneboundr) {
                                        # a function
                                        if (strand1 == 0) {
                                          start1 <- left1 - querygeneboundl + 1L
                                          stop1 <- right1 - querygeneboundl + 1L
                                        } else {
                                          start1 <- querygeneboundr - right1 + 1L
                                          stop1 <- querygeneboundr - left1 + 1L
                                        }
                                        if (strand2 == 0) {
                                          start2 <- left2 - subgeneboundl + 1L
                                          stop2 <- right2 - subgeneboundl + 1L
                                        } else {
                                          start2 <- subgeneboundr - right2 + 1L
                                          stop2 <- subgeneboundr - left2 + 1L
                                        }
                                        return(matrix(data = c(start1,
                                                               stop1,
                                                               start2,
                                                               stop2),
                                                      ncol = 1L))
                                      },
                                      i1 = hitsets[, 1L],
                                      i2 = hitsets[, 2L],
                                      strand1 = GeneCalls[[IDMatch[m1]]]$Strand[hitsets[, 1L]],
                                      strand2 = GeneCalls[[IDMatch[m2]]]$Strand[hitsets[, 2L]],
                                      querygeneboundl = GeneCalls[[IDMatch[m1]]]$Start[hitsets[, 1L]],
                                      querygeneboundr = GeneCalls[[IDMatch[m1]]]$Stop[hitsets[, 1L]],
                                      subgeneboundl = GeneCalls[[IDMatch[m2]]]$Start[hitsets[, 2L]],
                                      subgeneboundr = GeneCalls[[IDMatch[m2]]]$Stop[hitsets[, 2L]],
                                      left1 = SynExtendObject[[m2, m1]][, "QLeftPos"],
                                      left2 = SynExtendObject[[m2, m1]][, "SLeftPos"],
                                      right1 = SynExtendObject[[m2, m1]][, "QRightPos"],
                                      right2 = SynExtendObject[[m2, m1]][, "SRightPos"])
            
            WithinQueryNucs <- unname(split(x = WithinQueryNucs,
                                            f = rep(x = seq(nrow(PMatrix)),
                                                    times = SynExtendObject[[m1, m2]][, "TotalKmerHits"])))
            WithinQueryNucs <- lapply(X = WithinQueryNucs,
                                      FUN = function(x) {
                                        do.call(cbind, x)
                                      })
            for (m3 in seq_along(WithinQueryNucs)) {
              o1 <- order(WithinQueryNucs[[m3]][1L, , drop = FALSE])
              o2 <- order(WithinQueryNucs[[m3]][3L, , drop = FALSE])
              if (length(o1) > 1L) {
                if (o1[1L] > o1[2L]) {
                  WithinQueryNucs[[m3]][seq_len(2L), ] <- WithinQueryNucs[[m3]][seq_len(2L), o1, drop = FALSE]
                }
                
                if (o2[1L] > o2[2L]) {
                  WithinQueryNucs[[m3]][3:4, ] <- WithinQueryNucs[[m3]][3:4, o2, drop = FALSE]
                }
              }
            }
            # if there are any AA alignments to do
            WithinQueryAAs <- WithinQueryNucs[AASelect]
            if (nrow(AASubSet) > 0) {
              CurrentQCDSs <- QCDSCount[AASubSet[, 1L]]
              CurrentSCDSs <- SCDSCount[AASubSet[, 2L]]
              # loop through the nucleotide anchor positions
              # drop all anchors in cases where there are more than one CDS
              # drop anchors that are out of frame, or that are in NT space
              for (m3 in seq_along(WithinQueryAAs)) {
                if (CurrentQCDSs[m3] > 1L |
                    CurrentSCDSs[m3] > 1L) {
                  # print(m3)
                  WithinQueryAAs[[m3]] <- integer()
                } else {
                  # check the anchors then drop or adjust
                  size_check <- apply(X = WithinQueryAAs[[m3]],
                                      MARGIN = 2L,
                                      FUN = function(x) {
                                        x[c(FALSE, TRUE)] - x[c(TRUE, FALSE)] + 1L
                                      }) %% 3L == 0L
                  frame_check <- apply(X = WithinQueryAAs[[m3]][c(TRUE,FALSE), , drop = FALSE],
                                       MARGIN = 2L,
                                       FUN = function(x) {
                                         ((x + 2L) %% 3L) == 0L
                                       })
                  # size check is a matrix of logicals
                  # frame check is a matrix of logicals
                  # drop hits that aren't coding
                  # and drop hits that aren't in the appropriate frame
                  WithinQueryAAs[[m3]] <- WithinQueryAAs[[m3]][, (colSums(size_check) == 2L) &
                                                                 (colSums(frame_check) == 2L),
                                                               drop = FALSE]
                  # we need to offset the start, but not the end
                  # this can return an integer of length zero when nothing passes
                  WithinQueryAAs[[m3]] <- matrix(data = as.integer((WithinQueryAAs[[m3]] + c(2L, 0L)) / 3L),
                                                 nrow = nrow(WithinQueryAAs[[m3]]),
                                                 ncol = ncol(WithinQueryAAs[[m3]]))
                }
                
              } # end of m3 loop
            } # WithinQueryAAs logical check
            df_aa$Position <- WithinQueryAAs
            df_nt$Position <- WithinQueryNucs[NTSelect]
            # return(list(WithinQueryAAs,
            #             WithinQueryNucs,
            #             QueryAA,
            #             SubjectAA))
            
          } else {
            # don't set anchors at all
            WithinQueryAAs <- WithinQueryNucs <- list()
          }
          
          if (sum(AASelect) > 0) {
            aapairs <- AlignPairs(pattern = QueryAA,
                                  subject = SubjectAA,
                                  pairs = df_aa,
                                  verbose = FALSE)
            current_aa_pids <- aapairs$Matches / aapairs$AlignmentLength
            current_aa_scores <- aapairs$Score
          } else {
            current_aa_pids <- numeric()
            current_aa_scores <- numeric()
          }
          if (Verbose) {
            PBCount <- PBCount + 1L
            setTxtProgressBar(pb = pBar,
                              value = PBCount / Total)
          }
          
          if (sum(NTSelect) > 0) {
            ntpairs <- AlignPairs(pattern = QueryDNA,
                                  subject = SubjectDNA,
                                  pairs = df_nt,
                                  verbose = FALSE)
            current_nt_pids <- ntpairs$Matches / ntpairs$AlignmentLength
            current_nt_scores <- ntpairs$Score
          } else {
            current_nt_pids <- numeric()
            current_nt_scores <- numeric()
          }
          if (Verbose) {
            PBCount <- PBCount + 1L
            setTxtProgressBar(pb = pBar,
                              value = PBCount / Total)
          }
          
          vec1 <- vec2 <- vector(mode = "numeric",
                                 length = nrow(PMatrix))
          vec1[AASelect] <- current_aa_pids
          vec1[NTSelect] <- current_nt_pids
          vec2[AASelect] <- current_aa_scores
          vec2[NTSelect] <- current_nt_scores
          # For Testing
          # AA_Anchors[[Count]] <- WithinQueryAAs
          # NT_Anchors[[Count]] <- WithinQueryNucs[NTSelect]
          
        } else if (AlignmentFun == "AlignProfiles") {
          
          # build out a map of who is being called where
          # if we're aligning nucleotides, just call the position in the DNA
          # stringsets,
          # if you not, use the matched positions
          ws1 <- match(table = names(QueryAA),
                       x = names(QueryDNA)[PMatrix[, 1L]])
          ws2 <- match(table = names(SubjectAA),
                       x = names(SubjectDNA)[PMatrix[, 2L]])
          vec1 <- vec2 <- vector(mode = "numeric",
                                 length = length(NucDist))
          
          # if we're aligning everything as NTs, just cheat and change the coding
          # logical
          if (IgnoreDefaultStringSet) {
            QCode <- rep(FALSE,
                         length(QCode))
          }
          AASelect <- QCode[PMatrix[, 1L]] & QMod[PMatrix[, 1L]] & SCode[PMatrix[, 2L]] & SMod[PMatrix[, 2L]]
          # return(list("AASelect" = AASelect,
          #             "QCode" = QCode[PMatrix[, 1L]],
          #             "QMod" = QMod[PMatrix[, 1L]],
          #             "SCode" = SCode[PMatrix[, 2L]],
          #             "SMod" = SMod[PMatrix[, 2L]],
          #             "PMatrix" = PMatrix,
          #             "DataPool" = DataPool,
          #             "m1" = m1,
          #             "m2" = m2,
          #             "ws1" = ws1,
          #             "ws2" = ws2))
          for (m3 in seq_along(NucDist)) {
            
            NucDist[m3] <- sqrt(sum((nuc1[m3, ] - nuc2[m3, ])^2)) / ((sum(nuc1[m3, ]) + sum(nuc2[m3, ])) / 2)
            
            if (AASelect[m3]) {
              # align as amino acids
              ph1 <- AlignProfiles(pattern = QueryAA[ws1[m3]],
                                   subject = SubjectAA[ws2[m3]],
                                   p.struct = QueryStruct[ws1[m3]],
                                   s.struct = SubjectStruct[ws2[m3]])
              ph2 <- DistanceMatrix(myXStringSet = ph1,
                                    includeTerminalGaps = TRUE,
                                    type = "matrix",
                                    verbose = FALSE)
              ph3 <- ScoreAlignment(myXStringSet = ph1,
                                    structures = PredictHEC(myAAStringSet = ph1,
                                                            type = "probabilities",
                                                            HEC_MI1 = MAT1,
                                                            HEC_MI2 = MAT2),
                                    structureMatrix = structureMatrix)
            } else {
              # align as nucleotides
              ph1 <- AlignProfiles(pattern = QueryDNA[PMatrix[m3, 1L]],
                                   subject = SubjectDNA[PMatrix[m3, 2L]])
              ph2 <- DistanceMatrix(myXStringSet = ph1,
                                    includeTerminalGaps = TRUE,
                                    type = "matrix",
                                    verbose = FALSE)
              ph3 <- ScoreAlignment(myXStringSet = ph1,
                                    substitutionMatrix = substitutionMatrix)
            }
            vec1[m3] <- 1 - ph2[1, 2]
            vec2[m3] <- ph3
            
            if (Verbose) {
              PBCount <- PBCount + 1L
              setTxtProgressBar(pb = pBar,
                                value = PBCount / Total)
            }
          } # end m3 loop
          
        } # end if else on alignment function
        PH[[Count]] <- data.frame("p1" = names(QueryDNA)[PMatrix[, 1]],
                                  "p2" = names(SubjectDNA)[PMatrix[, 2]],
                                  "Consensus" = diff2,
                                  "p1featurelength" = QueryFeatureLength,
                                  "p2featurelength" = SubjectFeatureLength,
                                  "blocksize" = AbsBlockSize,
                                  "KDist" = NucDist,
                                  "TotalMatch" = TotalMatch,
                                  "MaxMatch" = MatchMax,
                                  "UniqueMatches" = UniqueMatches,
                                  "PID" = vec1,
                                  "Score" = vec2,
                                  "Alignment" = ifelse(test = AASelect,
                                                       yes = "AA",
                                                       no = "NT"),
                                  "Block_UID" = BlockID_Map)
      } else {
        # link table is not populated
        PH[[Count]] <- data.frame("p1" = character(),
                                  "p2" = character(),
                                  "Consensus" = numeric(),
                                  "p1featurelength" = integer(),
                                  "p2featurelength" = integer(),
                                  "blocksize" = integer(),
                                  "KDist" = numeric(),
                                  "TotalMatch" = integer(),
                                  "MaxMatch" = integer(),
                                  "UniqueMatches" = integer(),
                                  "PID" = numeric(),
                                  "Score" = numeric(),
                                  "Alignment" = character(),
                                  "Block_UID" = integer())
      }
      # Count and PBCount are unlinked,
      # iterate through both separately and correctly
      Count <- Count + 1L
      
      if (object.size(DataPool) > Storage) {
        sw1 <- sapply(X = DataPool,
                      FUN = function(x) {
                        is.null(x)
                      })
        sw2 <- which(sw1)
        if (length(sw2) > 0) {
          # bonk the first one ... this might realistically need to happen in a while loop, but for now we can live with this
          DataPool[[sw2[1L]]] <- NULL
        } else {
          # i don't know if this case can happen, but we're putting a print statement here just in case
          print("Unexpected case occured while managing storage. Please contact maintainer.")
        }
      }
      Prev_m1 <- m1
    } # end m2
  } # end m1
  res <- do.call(rbind,
                 PH)
  # return(res)
  All_UIDs <- unique(res$Block_UID)
  res$Block_UID[res$Block_UID == -1L] <- seq(from = max(All_UIDs) + 1L,
                                             by = 1L,
                                             length.out = sum(res$Block_UID == -1L))
  attr(x = res,
       which = "GeneCalls") <- GeneCalls
  class(res) <- c("data.frame",
                  "PairSummaries")
  attr(x = res,
       which = "AlignmentFunction") <- AlignmentFun
  attr(x = res,
       which = "KVal") <- KmerSize
  # attr(x = res,
  #      which = "NT_Anchors") <- NT_Anchors
  
  # close pBar and return res
  if (Verbose) {
    TimeEnd <- Sys.time()
    close(pBar)
    print(TimeEnd - TimeStart)
  }
  return(res)
}


