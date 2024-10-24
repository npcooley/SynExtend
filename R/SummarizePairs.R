###### -- summarize seqs ------------------------------------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu
# given a linked pairs object, and a FeatureSeqs object return a pairsummaries
# object
# this function will always align, unlike PairSummaries
# TODO:
# implement:
# ShowPlot
# Processors -- partially implemented, at least for AlignPairs
# ellipses
# SearchIndex subsetting to only search non-occupied spaces -- erik thinks we can avoid this

SummarizePairs <- function(SynExtendObject,
                           DataBase01,
                           IncludeIndexSearch = TRUE,
                           AlignmentFun = "AlignPairs",
                           RetainAnchors = TRUE,
                           DefaultTranslationTable = "11",
                           KmerSize = 5,
                           IgnoreDefaultStringSet = FALSE,
                           Verbose = FALSE,
                           ShowPlot = FALSE,
                           Processors = 1,
                           Storage = 2,
                           IndexParams = list("K" = 6),
                           SearchParams = list("perPatternLimit" = 1),
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
      
      # build our data pool first
      # regardless of how we align or if we're including search index or not, we need to prepare and collect
      # the same basic statistics and look aheads
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
        
        DataPool[[m1]]$struct <- PredictHEC(myAAStringSet = DataPool[[m1]]$AA,
                                            type = "probabilities",
                                            HEC_MI1 = MAT1,
                                            HEC_MI2 = MAT2)
        if (IncludeIndexSearch & !IgnoreDefaultStringSet) {
          # return(list(DataPool[[m1]]$AA,
          #             DataPool[[m1]]$mod,
          #             DataPool[[m1]]$code))
          DataPool[[m1]]$index <- do.call(what = "IndexSeqs",
                                          args = c(list("subject" = DataPool[[m1]]$AA[DataPool[[m1]]$mod[DataPool[[m1]]$code]],
                                                        "verbose" = FALSE),
                                                   IndexParams))
        } else if (IncludeIndexSearch & IgnoreDefaultStringSet) {
          DataPool[[m1]]$index <- do.call(what = "IndexSeqs",
                                          args = c(list("subject" = DataPool[[m1]]$DNA,
                                                        "verbose" = FALSE),
                                                   IndexParams))
        }
        
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
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
        
        DataPool[[m2]]$struct <- PredictHEC(myAAStringSet = DataPool[[m2]]$AA,
                                            type = "probabilities",
                                            HEC_MI1 = MAT1,
                                            HEC_MI2 = MAT2)
        if (IncludeIndexSearch & !IgnoreDefaultStringSet) {
          DataPool[[m2]]$index <- do.call(what = "IndexSeqs",
                                          args = c(list("subject" = DataPool[[m2]]$AA[DataPool[[m2]]$mod[DataPool[[m2]]$code]],
                                                        "verbose" = FALSE),
                                                   IndexParams))
        } else if (IncludeIndexSearch & IgnoreDefaultStringSet) {
          DataPool[[m2]]$index <- do.call(what = "IndexSeqs",
                                          args = c(list("subject" = DataPool[[m2]]$DNA,
                                                        "verbose" = FALSE),
                                                   IndexParams))
        }
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
      }
      
      if (Prev_m1 != m1) {
        QueryDNA <- DataPool[[m1]]$DNA
        QueryAA <- DataPool[[m1]]$AA
        QNTCount <- DataPool[[m1]]$len
        QMod <- DataPool[[m1]]$mod
        QCode <- DataPool[[m1]]$code
        QCDSCount <- DataPool[[m1]]$cds
        QueryStruct <- DataPool[[m1]]$struct
        QueryIndex <- DataPool[[m1]]$index
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
      SubjectIndex <- DataPool[[m2]]$index
      
      if (IncludeIndexSearch) {
        # step 1: build indexes into the data pool if they don't exist already
        # # this is already accomplished
        
        # step 2: run the searches
        
        if (IgnoreDefaultStringSet) {
          search_df1 <- do.call(what = "SearchIndex",
                                args = c(list("pattern" = QueryDNA,
                                              "invertedIndex" = SubjectIndex,
                                              "verbose" = FALSE,
                                              "processors" = Processors),
                                         SearchParams))
          search_df2 <- do.call(what = "SearchIndex",
                                args = c(list("pattern" = SubjectDNA,
                                              "invertedIndex" = QueryIndex,
                                              "verbose" = FALSE,
                                              "processors" = Processors),
                                         SearchParams))
        } else {
          search_df1 <- do.call(what = "SearchIndex",
                                args = c(list("pattern" = QueryAA[QMod[QCode]],
                                              "invertedIndex" = SubjectIndex,
                                              "verbose" = FALSE,
                                              "processors" = Processors),
                                         SearchParams))
          search_df2 <- do.call(what = "SearchIndex",
                                args = c(list("pattern" = SubjectAA[SMod[SCode]],
                                              "invertedIndex" = QueryIndex,
                                              "verbose" = FALSE,
                                              "processors" = Processors),
                                         SearchParams))
          # return(list("seq" = QueryAA,
          #             "code" = QCode,
          #             "mod" = QMod))
          #direction 1 is pattern -> subject
          #direction 2 is subject -> pattern
          search_df1 <- search_df1[order(search_df1$Pattern,
                                         search_df1$Subject), ]
          search_df2 <- search_df2[order(search_df2$Subject,
                                         search_df2$Pattern), ]
          # return(list(search_df1,
          #             search_df2))
          search_i1 <- paste(search_df1$Pattern,
                             search_df1$Subject,
                             sep = "_")
          search_i2 <- paste(search_df2$Subject,
                             search_df2$Pattern,
                             sep = "_")
          search_pairs <- data.frame("p1" = names(QueryAA[QMod[QCode]])[search_df1$Pattern[search_i1 %in% search_i2]],
                                     "p2" = names(SubjectAA[SMod[SCode]])[search_df1$Subject[search_i1 %in% search_i2]])
          place_holder1 <- do.call(rbind,
                                   strsplit(x = search_pairs$p1,
                                            split = "_",
                                            fixed = TRUE))
          search_pairs$i1 <- as.integer(place_holder1[, 2L])
          search_pairs$f1 <- as.integer(place_holder1[, 3L])
          place_holder2 <- do.call(rbind,
                                   strsplit(x = search_pairs$p2,
                                            split = "_",
                                            fixed = TRUE))
          search_pairs$i2 <- as.integer(place_holder2[, 2L])
          search_pairs$f2 <- as.integer(place_holder2[, 3L])
          search_pairs$f_hits <- search_df1$Position[search_i1 %in% search_i2]
          search_pairs$s1 <- GeneCalls[[m1]]$Strand[search_pairs$f1]
          search_pairs$s2 <- GeneCalls[[m2]]$Strand[search_pairs$f2]
          search_pairs$start1 <- GeneCalls[[m1]]$Start[search_pairs$f1]
          search_pairs$start2 <- GeneCalls[[m2]]$Start[search_pairs$f2]
          search_pairs$stop1 <- GeneCalls[[m1]]$Stop[search_pairs$f1]
          search_pairs$stop2 <- GeneCalls[[m2]]$Stop[search_pairs$f2]
          
          # slam everything together -- i.e. build out rows that need to be
          # added to the linked pairs object
          hit_adjust_start <- do.call(cbind,
                                      search_pairs$f_hits)
          hit_key <- vapply(X = search_pairs$f_hits,
                            FUN = function(x) {
                              ncol(x)
                            },
                            FUN.VALUE = vector(mode = "integer",
                                               length = 1L))
          hit_q_partner <- rep(search_pairs$f1,
                               times = hit_key)
          hit_s_partner <- rep(search_pairs$f2,
                               times = hit_key)
          # should be a list in the same shape as the SearchIndex Positions
          # but with the hits mapped to (mostly) the right spots
          # !!! IMPORTANT !!! This is limited to sequences without introns
          # for this to work correctly with introns, I need to be able to divy
          # these offsets up across CDSs, which though possible now will need
          # some significant infrastructure changes to make work cleanly
          
          # return(list("a" = search_pairs$f_hits,
          #             "d" = search_pairs$s1,
          #             "e" = search_pairs$s2,
          #             "f" = search_pairs$start1,
          #             "g" = search_pairs$start2,
          #             "h" = search_pairs$stop1,
          #             "i" = search_pairs$stop2))
          
          hit_arrangement <- mapply(SIMPLIFY = FALSE,
                                    FUN = function(a, d, e, f, g, h, i) {
                                      # 0 == FALSE
                                      # 1 == TRUE
                                      # strand is 0 or 1, 1 == negative strand
                                      if (d) {
                                        # negative strand, flip, offset, then assign
                                        # feature starts at 100
                                        # 1 == 100
                                        # 2 == 97
                                        # 3 == 94
                                        # etc
                                        # right - (x - 1) * 3 = position
                                        h1 <- h - (a[c(1, 2), , drop = FALSE] - 1) * 3
                                        h1 <- apply(X = h1,
                                                    MARGIN = 2,
                                                    FUN = function(x) {
                                                      sort(x)
                                                    },
                                                    simplify = TRUE)
                                      } else {
                                        # positive strand, just offset and assign
                                        # feature starts at 10,
                                        # 1 == 10
                                        # 2 == 13
                                        # 3 == 14 
                                        # etc
                                        # left + (x - 1) * 3 = position
                                        h1 <- f + (a[c(1, 2), , drop = FALSE] - 1) * 3
                                        
                                      }
                                      # repeat for feature 2
                                      if (e) {
                                        h2 <- i - (a[c(3,4), , drop = FALSE] - 1) * 3
                                        h2 <- apply(X = h2,
                                                    MARGIN = 2,
                                                    FUN = function(x) {
                                                      sort(x)
                                                    },
                                                    simplify = TRUE)
                                      } else {
                                        h2 <- g + (a[c(3,4), , drop = FALSE] - 1) * 3
                                      }
                                      return(rbind(h1, h2))
                                    },
                                    a = search_pairs$f_hits,
                                    d = search_pairs$s1,
                                    e = search_pairs$s2,
                                    f = search_pairs$start1,
                                    g = search_pairs$start2,
                                    h = search_pairs$stop1,
                                    i = search_pairs$stop2)
          
          hit_rearrangement <- t(do.call(cbind,
                                         hit_arrangement))
          block_bounds <- lapply(X = hit_arrangement,
                                 FUN = function(x) {
                                   c(min(x[1, ]),
                                     max(x[2, ]),
                                     min(x[3, ]),
                                     max(x[4, ]))
                                 })
          block_bounds <- do.call(rbind,
                                  block_bounds)
          
          hit_widths <- lapply(X = search_pairs$f_hits,
                               FUN = function(x) {
                                 x[2, ] - x[1, ] + 1L
                               })
          hit_totals <- vapply(X = hit_widths,
                               FUN = function(x) {
                                 sum(x)
                               },
                               FUN.VALUE = vector(mode = "integer",
                                                  length = 1))
          hit_max <- vapply(X = hit_widths,
                            FUN = function(x) {
                              max(x)
                            },
                            FUN.VALUE = vector(mode = "integer",
                                               length = 1))
          # return(list("QueryGene" = hit_q_partner,
          #             "SubjectGene" = hit_s_partner,
          #             "ExactOverlap" = unlist(hit_widths),
          #             "QueryIndex" = rep(search_pairs$i1,
          #                                times = hit_key),
          #             "SubjectIndex" = rep(search_pairs$i2,
          #                                  times = hit_key),
          #             "QLeftPos" = hit_rearrangement[, 1],
          #             "QRightPos" = hit_rearrangement[, 2],
          #             "SLeftPos" = hit_rearrangement[, 3],
          #             "SRightPos" = hit_rearrangement[, 4]))
          add_by_hit <- data.frame("QueryGene" = hit_q_partner,
                                   "SubjectGene" = hit_s_partner,
                                   "ExactOverlap" = unlist(hit_widths),
                                   "QueryIndex" = rep(search_pairs$i1,
                                                      times = hit_key),
                                   "SubjectIndex" = rep(search_pairs$i2,
                                                        times = hit_key),
                                   "QLeftPos" = hit_rearrangement[, 1],
                                   "QRightPos" = hit_rearrangement[, 2],
                                   "SLeftPos" = hit_rearrangement[, 3],
                                   "SRightPos" = hit_rearrangement[, 4])
          add_by_block <- data.frame("QueryGene" = search_pairs$f1,
                                     "SubjectGene" = search_pairs$f2,
                                     "ExactOverlap" = hit_totals,
                                     "QueryIndex" = search_pairs$i1,
                                     "SubjectIndex" = search_pairs$i2,
                                     "QLeftPos" = block_bounds[, 1],
                                     "QRightPos" = block_bounds[, 2],
                                     "SLeftPos" = block_bounds[, 3],
                                     "SRightPos" = block_bounds[, 4],
                                     "MaxKmerSize" = hit_max,
                                     "TotalKmerHits" = hit_key)
          
          # using forward hits from here:
          # append pairs that don't appear in the current linked pairs object
          # onto the current linked pairs object
        }
        
        # hits are now transposed into the context of the whole sequence
        # as opposed to the feature
        # build out rows to add, and then add them
        if (nrow(SynExtendObject[[m1, m2]]) > 0) {
          select_row1 <- paste(SynExtendObject[[m1, m2]][, 1L],
                               SynExtendObject[[m1, m2]][, 4L],
                               SynExtendObject[[m1, m2]][, 2L],
                               SynExtendObject[[m1, m2]][, 5L],
                               sep = "_")
          select_row2 <- paste(SynExtendObject[[m2, m1]][, 1L],
                               SynExtendObject[[m2, m1]][, 4L],
                               SynExtendObject[[m2, m1]][, 2L],
                               SynExtendObject[[m2, m1]][, 5L],
                               sep = "_")
        } else {
          select_row1 <- select_row2 <- vector(mode = "character",
                                               length = 0)
        }
        # upper diagonal is blocks,
        # lower diagonal is hits
        select_row3 <- paste(add_by_block$QueryGene,
                             add_by_block$QueryIndex,
                             add_by_block$SubjectGene,
                             add_by_block$SubjectIndex,
                             sep = "_")
        select_row4 <- paste(add_by_hit$QueryGene,
                             add_by_hit$QueryIndex,
                             add_by_hit$SubjectGene,
                             add_by_hit$SubjectIndex,
                             sep = "_")
        sr_upper <- !(select_row3 %in% select_row1)
        sr_lower <- !(select_row4 %in% select_row2)
        
        # return(list(select_row1,
        #             select_row2,
        #             select_row3,
        #             select_row4,
        #             SynExtendObject[[m1, m2]],
        #             SynExtendObject[[m2, m1]],
        #             add_by_hit,
        #             add_by_block))
        if (any(sr_upper)) {
          SynExtendObject[[m1, m2]] <- rbind(SynExtendObject[[m1, m2]],
                                             as.matrix(add_by_block[sr_upper, ]))
          SynExtendObject[[m2, m1]] <- rbind(SynExtendObject[[m2, m1]],
                                             as.matrix(add_by_hit[sr_lower, ]))
        }
        
        
        
        # step 3: morph the searches into the SynExtendObject so other things
        # go smoothly
        
        # alignments happen later and profile vs pairs is chosen by the user
        
      } # end include index search logical
      
      # return(list(add_by_hit,
      #             add_by_block))
      
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
        # we're eventually moving blocksize stuff down to happen post everything else
        
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
            hitsets <- SynExtendObject[[m2, m1]][, c(1, 2), drop = FALSE]
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
                  WithinQueryNucs[[m3]][c(1,2), ] <- WithinQueryNucs[[m3]][c(1,2), o1, drop = FALSE]
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
                  # ensure that anchors are in ascending order
                }
                
              } # end of m3 loop
            } # WithinQueryAAs logical check
            # df_aa$Position <- WithinQueryAAs
            # ph_obj <- WithinQueryAAs
            if (IncludeIndexSearch) {
              full_aa_set <- paste(df_aa$Pattern,
                                   df_aa$Subject,
                                   sep = "_")
              # this isn't right because these need to be indexed to the
              # available coding sequences
              aa_w_index <- paste(match(x = search_pairs$f1,
                                        table = aa_match1),
                                  match(x = search_pairs$f2,
                                        table = aa_match2),
                                  sep = "_")
              # return(list("a" = full_aa_set,
              #             "b" = aa_w_index))
              WithinQueryAAs[match(x = aa_w_index,
                                   table = full_aa_set)] <- search_pairs$f_hits
            }
            df_aa$Position <- WithinQueryAAs
            df_aa$Position <- mapply(SIMPLIFY = FALSE,
                                     FUN = function(x, y, z) {
                                       if (length(x) > 4) {
                                         cbind(matrix(0L,
                                                      4),
                                               x,
                                               matrix(c(y, y, z, z),
                                                      4))
                                       } else {
                                         cbind(matrix(0L,
                                                      4),
                                               matrix(data = c(y,y,z,z),
                                                      4))
                                       }
                                     },
                                     x = df_aa$Position,
                                     y = width(QueryAA[df_aa$Pattern]) + 1L,
                                     z = width(SubjectAA[df_aa$Subject]) + 1L)
            
            df_nt$Position <- WithinQueryNucs[NTSelect]
            df_nt$Position <- mapply(SIMPLIFY = FALSE,
                                     FUN = function(x, y, z) {
                                       if (length(x) > 4) {
                                         cbind(matrix(0L,
                                                      4),
                                               x,
                                               matrix(c(y, y, z, z),
                                                      4))
                                       } else {
                                         cbind(matrix(0L,
                                                      4),
                                               matrix(data = c(y,y,z,z),
                                                      4))
                                       }
                                       
                                     },
                                     x = df_nt$Position,
                                     y = QNTCount[df_nt$Pattern] + 1L,
                                     z = SNTCount[df_nt$Subject] + 1L)
            # return(list(WithinQueryAAs,
            #             WithinQueryNucs,
            #             QueryAA,
            #             SubjectAA))
            # one last check to reject hits that alignpairs views as overlapping
            for (m3 in seq_along(df_aa$Position)) {
              while (is.unsorted(df_aa$Position[[m3]][1, ]) |
                     is.unsorted(df_aa$Position[[m3]][3, ])) {
                hit_sizes <- df_aa$Position[[m3]][2, ] - df_aa$Position[[m3]][1, ] + 1L
                # drop the smallest hit that is not an anchor
                z1 <- which.min(hit_sizes[hit_sizes > 1])
                # i'm not sure if this is safe, but this will be refactored to
                # occur before anchoring anyway
                if (length(z1) < 1) {
                  stop("an unexpected condition occurred, please contact maintainer")
                } else {
                  df_aa$Position[[m3]] <- df_aa$Position[[m3]][, -(z1 + 1L), drop = FALSE]
                }
              }
            }
          } else {
            # don't set anchors at all
            WithinQueryAAs <- WithinQueryNucs <- list()
          }
          
          if (sum(AASelect) > 0) {
            # return(list("q" = QueryAA,
            #             "s" = SubjectAA,
            #             "df" = df_aa,
            #             "ph" = ph_obj))
            aapairs <- AlignPairs(pattern = QueryAA,
                                  subject = SubjectAA,
                                  pairs = df_aa,
                                  verbose = FALSE,
                                  processors = Processors)
            current_aa_pids <- aapairs$Matches / aapairs$AlignmentLength
            current_aa_scores <- aapairs$Score / aapairs$AlignmentLength
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
                                  verbose = FALSE,
                                  processors = Processors)
            current_nt_pids <- ntpairs$Matches / ntpairs$AlignmentLength
            current_nt_scores <- ntpairs$Score / ntpairs$AlignmentLength
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
          # return(list(aapairs,
          #             ntpairs))
          
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


