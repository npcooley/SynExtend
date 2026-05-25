# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: Nicholas.Cooley@ul.ie

NucleotideOverlap <- function(SyntenyObject,
                              GeneCalls,
                              LimitIndex = FALSE,
                              AcceptContigNames = TRUE,
                              Verbose = FALSE) {
  ######
  # Error Checking
  # Require names for synteny object, DECIPHER be loaded, object sizes, some index checking
  # Parameters Default Filter, FilterOverlap, and Filter Coverage have been depreciated and
  # a filtering function is under construction
  ######
  
  L <- nrow(SyntenyObject)
  if (!("DECIPHER" %in% .packages())) {
    stop ("Required package DECIPHER is not loaded.")
  }
  if (!is(SyntenyObject, "Synteny")) {
    stop ("Object is not a synteny object.")
  }
  # this needs to be dropped if i'm going to work with self-synteny objects
  # if (L != length(GeneCalls)) {
  #   stop ("Synteny Object and Gene Calls have differing number of genomes.")
  # }
  if (is.null(names(GeneCalls))) {
    stop ("Gene Predictions must be named.")
  }
  if (is.null(rownames(SyntenyObject))) {
    stop ("Synteny Object must have named identifiers.")
  }
  if (any(names(GeneCalls) != rownames(SyntenyObject))) {
    stop ("Names between Synteny Object and Gene Predictions do not match.")
  }
  if (!all(rownames(SyntenyObject) %in% names(GeneCalls))) {
    stop ("All distinct identifiers in the Synteny object require matching gene calls.")
  }
  if (L <= 1L) {
    stop ("SyntenyObject is too small.")
  }
  
  if (Verbose) {
    TotalTimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  ResultMatrix <- matrix(data = list(),
                         nrow = L,
                         ncol = L)
  # The formula for the count of the upper and lower triangles minus the diagonal...
  # TotalLength <- L^2 - L
  # we're only ticking the pb once per comparison now, used to be twice ...
  TotalLength <- (L * (L - 1)) / 2
  TotalCounter <- 1L
  ######
  # scroll through every hit table in the synteny object
  ######
  
  ###### -- Deal with Different GeneCall types --------------------------------
  
  SynNames <- unname(sapply(diag(SyntenyObject),
                            function(x) gsub(x = names(x),
                                             pattern = " .*",
                                             replacement = ""),
                            simplify = FALSE,
                            USE.NAMES = FALSE))
  # this can no longer be length L and needs to just be the length of GeneCalls
  FeatureRepresentations <- ContigNames <- vector(mode = "list",
                                                  length = length(GeneCalls))
  names(FeatureRepresentations) <- names(GeneCalls)
  
  GRangeCheck <- vapply(X = GeneCalls,
                        FUN = function(x) {
                          is(object = x,
                             class2 = "GRanges")
                        },
                        FUN.VALUE = vector(mode = "logical",
                                           length = 1))
  if (any(GRangeCheck)) {
    warning("GRanges support is tenuous, exercise caution.")
  }
  IndexMatching <- vector("integer",
                          length = length(GeneCalls))
  
  if (Verbose) {
    cat("\nReconciling genecalls.\n")
  }
  for (m1 in seq_along(GeneCalls)) {
    if (!is(GeneCalls[[m1]],
            "GRanges")) {
      IndexMatching[m1] <- length(unique(GeneCalls[[m1]][, "Index"]))
    } else {
      IndexMatching[m1] <- 1L
    }
  }
  SyntenyIndices <- unname(lengths(diag(SyntenyObject)))
  if (any(IndexMatching > SyntenyIndices)) {
    warning ("Indices do not match. Setting LimitIndex to TRUE.")
    LimitIndex <- TRUE
  }
  
  for (m1 in seq_along(GeneCalls)) {
    if (is(GeneCalls[[m1]],
           "GRanges")) {
      # if GRanges, force contig name matching
      # stop if a contig name exists in the synteny object that
      # does not exist in the GRanges object
      
      ContigNames[[m1]] <- seq(length(GeneCalls[[m1]]@seqnames@lengths))
      names(ContigNames[[m1]]) <- unique(GeneCalls[[m1]]@seqnames@values)[ContigNames[[m1]]]
      
      if (any(is.na(match(x = names(ContigNames[[m1]]),
                          table = SynNames[[m1]])))) {
        stop ("Contig names imply inorrectly matched objects at diag position ",
                     names(GeneCalls)[m1])
      }
      
      TypePlaceHolder <- as.character(GeneCalls[[m1]]$type)
      GeneCalls[[m1]] <- GeneCalls[[m1]][TypePlaceHolder %in% c("gene",
                                                                "pseudogene"), ]
      StrandConversion <- ifelse(test = as.character(GeneCalls[[m1]]@strand == "+"),
                                 yes = 0L,
                                 no = 1L)
      StartConversion <- GeneCalls[[m1]]@ranges@start
      StopConversion <- GeneCalls[[m1]]@ranges@start + GeneCalls[[m1]]@ranges@width - 1L
      # assign integer positions to contings
      IndexConversion <- rep(unname(ContigNames[[m1]]),
                             GeneCalls[[m1]]@seqnames@lengths)
      
      LengthsConversion <- StopConversion - StartConversion + 1L
      # reorder and reindex lines
      C.Index <- IndexConversion
      ph1 <- unname(ContigNames[[m1]])
      ph2 <- match(x = names(ContigNames[[m1]]),
                   table = SynNames[[m1]])
      ph3 <- vector(mode = "list",
                    length = length(ph1))
      for (m3 in seq_along(ph1)) {
        ph3[[m3]] <- which(C.Index == ph1[m3])
      }
      for (m3 in seq_along(ph1)) {
        if (length(ph3[[m3]]) > 0L) {
          C.Index <- replace(x = C.Index,
                             list = ph3[[m3]],
                             values = ph2[m3])
        }
      }
      rm(list = c("ph1",
                  "ph2",
                  "ph3"))
      o <- order(C.Index,
                 StartConversion)
      
      StrandConversion <- StrandConversion[o]
      IndexConversion <- C.Index[o]
      StartConversion <- StartConversion[o]
      StopConversion <- StopConversion[o]
      LengthsConversion <- LengthsConversion[o]
      
      StrandMax <- rep(unname(SyntenyObject[[m1, m1]]),
                       table(IndexConversion))
      R <- mapply(function(x, y) IRanges(start = x,
                                         end = y),
                  x = StartConversion,
                  y = StopConversion,
                  SIMPLIFY = FALSE)
      FeatureRepresentations[[m1]] <- DataFrame("Index" = IndexConversion,
                                                "Strand" = StrandConversion,
                                                "Start" = StartConversion,
                                                "Stop" = StopConversion,
                                                # "Lengths" = LengthsConversion,
                                                "Type" = rep("gene",
                                                             length(o)),
                                                "Translation_Table" = rep(NA_character_,
                                                                          length(o)),
                                                "Coding" = rep(FALSE,
                                                               length(o)),
                                                "Range" = IRangesList(R))
      
      FeatureRepresentations[[m1]] <- FeatureRepresentations[[m1]][StopConversion <= StrandMax, ]
      rm(list = c("IndexConversion",
                  "StrandConversion",
                  "StartConversion",
                  "StopConversion",
                  "LengthsConversion"))
      
      ph <- seq(length(ContigNames[[m1]]))
      # ph <- unname(SyntenyObject[[m1, m1]])
      names(ph) <- SynNames[[m1]]
      ContigNames[[m1]] <- ph
      
    } else if (is(GeneCalls[[m1]],
                  "DFrame")) {
      
      ph <- unique(GeneCalls[[m1]][, c("Index", "Contig")])
      ContigNames[[m1]] <- ph$Index
      names(ContigNames[[m1]]) <- ph$Contig
      rm(list = c("ph"))
      
      # if this is to work with the sub-features as opposed to the entire feature,
      # i need to:
      # slam the IRanges list together,
      # rep the strand, and indices identifiers to expand to the new length that the need
      # generate a key 
      CurrentIndices <- as.integer(GeneCalls[[m1]]$Index)
      CurrentStarts <- as.integer(GeneCalls[[m1]]$Start)
      if (AcceptContigNames) {
        C.Index <- CurrentIndices
        ph1 <- unname(ContigNames[[m1]])
        ph2 <- match(x = names(ContigNames[[m1]]),
                     table = SynNames[[m1]])
        ph3 <- vector(mode = "list",
                      length = length(ph1))
        for (m3 in seq_along(ph1)) {
          ph3[[m3]] <- which(C.Index == ph1[m3])
        }
        for (m3 in seq_along(ph1)) {
          if (length(ph3) > 0L) {
            C.Index <- replace(x = C.Index,
                               list = ph3[[m3]],
                               values = ph2[m3])
          }
        }
        rm(list = c("ph1",
                    "ph2",
                    "ph3"))
        o <- order(C.Index,
                   CurrentStarts)
      } else {
        o <- order(CurrentIndices,
                   CurrentStarts)
      }
      # CurrentStarts <- CurrentStarts[o]
      # CurrentStops <- as.integer(GeneCalls[[m1]]$Stop)[o]
      # CurrentStrands <- as.integer(GeneCalls[[m1]]$Strand)[o]
      if (AcceptContigNames) {
        # CurrentIndices <- C.Index[o]
        GeneCalls[[m1]][, "Index"] <- C.Index
      } else {
        # CurrentIndices <- CurrentIndices[o]
        GeneCalls[[m1]][, "Index"] <- CurrentIndices
      }
      # this will return a single IRanges object
      ans <- GeneCalls[[m1]][o, ]
      
      SubFeatureBlocks <- FrameDownward(genecalls = ans)
      FeatureRepresentations[[m1]] <- SubFeatureBlocks
      
    } else if (is(GeneCalls[[m1]],
                  "Genes")) {
      # convert Erik's gene calls to a temporary DataFrame
      # the column "Gene" assigns whether or not that particular line is the gene
      # that the caller actually picked, calls must be subset to where Gene == 1
      ans <- GeneCalls[[m1]]
      
      CurrentIndices <- as.integer(ans[, "Index"])
      CurrentStarts <- as.integer(ans[, "Begin"])
      CurrentStrand <- as.integer(ans[, "Strand"])
      CurrentStops <- as.integer(ans[, "End"])
      CurrentType <- rep("gene",
                         nrow(ans))
      CurrentGene <- as.integer(ans[, "Gene"])
      CurrentCoding <- ifelse(test = CurrentGene > 0L,
                              yes = TRUE,
                              no = FALSE)
      R <- mapply(function(x, y) IRanges(start = x,
                                         end = y),
                  x = ans[, "Begin"],
                  y = ans[, "End"],
                  SIMPLIFY = FALSE)
      ContigNames[[m1]] <- unique(CurrentIndices)
      names(ContigNames[[m1]]) <- gsub(x = names(attr(ans, "widths")),
                                       pattern = " .*",
                                       replacement = "")
      if (AcceptContigNames) {
        C.Index <- CurrentIndices
        ph1 <- unname(ContigNames[[m1]])
        ph2 <- match(x = names(ContigNames[[m1]]),
                     table = SynNames[[m1]])
        ph3 <- vector(mode = "list",
                      length = length(ph1))
        for (m3 in seq_along(ph1)) {
          ph3[[m3]] <- which(C.Index == ph1[m3])
        }
        for (m3 in seq_along(ph1)) {
          if (length(ph3) > 0L) {
            C.Index <- replace(x = C.Index,
                               list = ph3[[m3]],
                               values = ph2[m3])
          }
        }
        ContigNames[[m1]] <- ContigNames[[m1]][ph2]
        rm(list = c("ph1",
                    "ph2",
                    "ph3"))
        o <- order(C.Index,
                   CurrentStarts)
      } else {
        # there should be no need to generically re-order a GeneCalls object
        # unless it has been fiddled with by the user, but we'll include this anyway
        o <- order(CurrentIndices,
                   CurrentStarts)
      }
      D <- DataFrame("Index" = CurrentIndices,
                     "Strand" = CurrentStrand,
                     "Start" = CurrentStarts,
                     "Stop" = CurrentStops,
                     "Type" = CurrentType,
                     "Range" = IRangesList(R),
                     "Gene" = CurrentGene,
                     "Coding" = CurrentCoding,
                     "Translation_Table" = rep(NA_character_,
                                               length(o)),
                     "Contig" = rep(names(ContigNames[[m1]]),
                                    times = table(CurrentIndices)))
      D <- D[o, ]
      D <- D[as.vector(ans[, "Gene"]) != 0L, ]
      rownames(D) <- NULL
      FeatureRepresentations[[m1]] <- D
      # ContigNames[[m1]] <- unique(D$Index)
      # names(ContigNames[[m1]]) <- ContigNames[[m1]]
      
      
      rm(list = c("D",
                  "ans",
                  "R",
                  "CurrentIndices",
                  "CurrentStarts",
                  "CurrentStrand",
                  "CurrentStops",
                  "CurrentType",
                  "CurrentGene",
                  "CurrentCoding"))
    }
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / length(GeneCalls))
    }
  }
  
  if (Verbose) {
    cat("\nFinding connected features.\n")
  }
  
  if (AcceptContigNames) {
    for (m1 in seq_along(ContigNames)) {
      if (any(is.na(match(x = names(ContigNames[[m1]]),
                          table = SynNames[[m1]])))) {
        stop ("Contig names imply incorrectly matched objects at diag position ",
                     names(GeneCalls)[m1])
      } else {
        ph <- unname(SyntenyObject[[m1, m1]])
        # ph <- seq(length(ContigNames[[m1]]))
        names(ph) <- SynNames[[m1]]
        ContigNames[[m1]] <- ph
      }
    }
  } else {
    for (m1 in seq_along(ContigNames)) {
      ph <- unname(SyntenyObject[[m1, m1]])
      names(ph) <- SynNames[[m1]]
      ContigNames[[m1]] <- ph
    }
  }
  diag(ResultMatrix) <- ContigNames
  feature_match <- match(x = rownames(SyntenyObject),
                         table = names(FeatureRepresentations))
  
  # temp_res <<- FeatureRepresentations
  ###### -- End Gene call stuff -----------------------------------------------
  
  pBar <- txtProgressBar(style = 1L)
  for (m1 in seq_len(L - 1L)) {
    for (m2 in (m1 +1L):L) {
      
      CurrentHitTable <- SyntenyObject[m1, m2, drop = FALSE][[1]]
      # hit table
      hit_table <- data.frame("idx1" = CurrentHitTable[, "index1"],
                              "idx2" = CurrentHitTable[, "index2"],
                              "inv" = CurrentHitTable[, "strand"],
                              "width" = CurrentHitTable[, "width"],
                              "doc1_l" = CurrentHitTable[, "start1"],
                              "doc1_r" = CurrentHitTable[, "start1"] + CurrentHitTable[, "width"] - 1L,
                              "doc2_l" = rep(-1L,
                                             nrow(CurrentHitTable)),
                              "doc2_r" = rep(-1L,
                                             nrow(CurrentHitTable)))
      inv_w <- hit_table$inv == 0
      hit_table$doc2_l[inv_w] <- CurrentHitTable[inv_w, "start2"]
      hit_table$doc2_r[inv_w] <- CurrentHitTable[inv_w, "start2"] + CurrentHitTable[inv_w, "width"] - 1L
      hit_table$doc2_l[!inv_w] <- CurrentHitTable[!inv_w, "start2"] - CurrentHitTable[!inv_w, "width"] + 1L
      hit_table$doc2_r[!inv_w] <- CurrentHitTable[!inv_w, "start2"]
      # genecalls in query
      d1_genecalls_source <- FeatureRepresentations[[feature_match[m1]]]
      d1_genecalls <- split(x = d1_genecalls_source,
                            f = d1_genecalls_source$Index)
      d1_ranges <- lapply(X = d1_genecalls,
                          FUN = function(x) {
                            IRanges(start = x$Left,
                                    end = x$Right)
                          })
      # genecalls in subject
      d2_genecalls_source <- FeatureRepresentations[[feature_match[m2]]]
      d2_genecalls <- split(x = d2_genecalls_source,
                            f = d2_genecalls_source$Index)
      d2_ranges <- lapply(X = d2_genecalls,
                          FUN = function(x) {
                            IRanges(start = x$Left,
                                    end = x$Right)
                          })
      # list names should correspond to the indices of the matched orders of 
      # the contigs
      
      doc1_ranges <- IRanges(start = hit_table$doc1_l,
                             end = hit_table$doc1_r)
      # split by the index1 value from the synteny object
      doc1_ranges <- split(x = doc1_ranges,
                           f = hit_table$idx1)
      # also split the hit data
      index_splits <- split(x = hit_table,
                            f = hit_table$idx1)
      sub_tbl <- vector(mode = "list",
                        length = length(doc1_ranges))
      
      for (d3 in seq_along(doc1_ranges)) {
        # using names here should be fine, because it should be forced to link
        # correctly to the contig names
        # ctg_w <- which(d1_genecalls$Index == as.integer(names(doc1_ranges)[d3]))
        # print(d3)
        # return(list(d1_genecalls,
        #             doc1_ranges,
        #             d3))
        w3 <- which(names(d1_genecalls) == names(doc1_ranges)[d3])
        if (length(w3) < 1) {
          warning("internal name alignment appears to be failing")
          next
        }
        w4 <- which(names(d1_ranges) == names(doc1_ranges)[d3])
        if (length(w4) < 1) {
          warning("internal name alignment appears to be failing")
          next
        }
        w5 <- which(names(index_splits) == names(doc1_ranges)[d3])
        if (length(w5) < 1) {
          warning("internal name alignment appears to be failing")
          next
        }
        curr_d1_genecalls <- d1_genecalls[[w3]]
        curr_d1_ranges <- d1_ranges[[w4]]
        curr_hit_tbl <- index_splits[[w5]]
        ph1 <- findOverlaps(query = doc1_ranges[[d3]],
                            subject = curr_d1_ranges)
        if (length(ph1) < 1) {
          # next in the d3 loop, not the m2 or m1
          next
        }
        ph2 <- overlapsRanges(query = doc1_ranges[[d3]],
                              subject = curr_d1_ranges,
                              hits = ph1)
        curr_from <- from(ph1)
        curr_to <- to(ph1)
        # find the range adjustment
        # ctg_w is an integer, and shouldn't be abnormal as long as the overhead
        # checking has done its job
        # FOR data.frames: the rownames are split as well, so we can reliably
        # access them for sorting, ordering or chaining back to the original object if we need
        # THIS IS NOT THE CASE FOR DataFrames
        adj_vals <- data.frame("l" = start(ph2) - start(doc1_ranges[[d3]][curr_from]),
                               "r" = end(doc1_ranges[[d3]][curr_from]) - end(ph2),
                               "hit_idx" = as.integer(rownames(curr_hit_tbl)[curr_from]),
                               "feature" = as.integer(rownames(curr_d1_genecalls)[curr_to]))
        # return(list(d1_genecalls_source,
        #             d1_genecalls,
        #             d1_ranges,
        #             curr_d1_ranges,
        #             ph1,
        #             curr_hit_tbl,
        #             adj_vals))
        # the row name of the current gene calls != the feature key
        # they're being left here because we need the index, but later that
        # index needs to be converted for both the d1 and d2 genecalls
        adj_tbl <- data.frame("idx1" = curr_hit_tbl$idx1[curr_from],
                              "idx2" = curr_hit_tbl$idx2[curr_from],
                              "inv" = curr_hit_tbl$inv[curr_from],
                              "width" = width(ph2),
                              "doc1_l" = curr_hit_tbl$doc1_l[curr_from] + adj_vals$l,
                              "doc1_r" = curr_hit_tbl$doc1_r[curr_from] - adj_vals$r,
                              "doc2_l" = rep(x = -1,
                                             times = nrow(adj_vals)),
                              "doc2_r" = rep(x = -1,
                                             times = nrow(adj_vals)),
                              "source_hit" = adj_vals$hit_idx,
                              "range1_hit" = adj_vals$feature)
        # forward-forward search space vs forward-reverse search space
        
        adj_tbl$doc2_l[adj_tbl$inv == 0] <- curr_hit_tbl$doc2_l[curr_from][adj_tbl$inv == 0] + adj_vals$l[adj_tbl$inv == 0]
        adj_tbl$doc2_l[adj_tbl$inv != 0] <- curr_hit_tbl$doc2_l[curr_from][adj_tbl$inv != 0] + adj_vals$r[adj_tbl$inv != 0]
        adj_tbl$doc2_r[adj_tbl$inv == 0] <- curr_hit_tbl$doc2_r[curr_from][adj_tbl$inv == 0] - adj_vals$r[adj_tbl$inv == 0]
        adj_tbl$doc2_r[adj_tbl$inv != 0] <- curr_hit_tbl$doc2_r[curr_from][adj_tbl$inv != 0] - adj_vals$l[adj_tbl$inv != 0]
        
        sub_tbl[[d3]] <- adj_tbl
        # return(list(hit_table,
        #             doc1_ranges,
        #             d1_genecalls,
        #             d2_genecalls,
        #             ph1,
        #             ph2,
        #             adj_tbl))
        
      }
      # drop NULL positions
      # do.call sets its behavior based on either the class or typeof
      # the first element, if it is NULL, that's a problem for you
      check_this <- vapply(X = sub_tbl,
                           FUN = function(x) {
                             is.null(x)
                           },
                           FUN.VALUE = vector(mode = "logical",
                                              length = 1L))
      sub_tbl <- sub_tbl[!check_this]
      # if there's only 1 tbl left, that's it,
      # if there are none, fill out a dummy table to move forward with
      # if there are many, slam them together
      if (length(sub_tbl) == 1) {
        hit_table <- sub_tbl[[1]]
      } else if (length(sub_tbl) > 1) {
        hit_table <- do.call(rbind,
                             sub_tbl)
        # need to reorder this
        # hits that have the same left bound will be ordered from shortest
        # to longest
        o_hits <- order(hit_table$idx1,
                        hit_table$idx2,
                        hit_table$doc1_l,
                        hit_table$doc1_r,
                        hit_table$doc2_l,
                        hit_table$doc2_r,
                        decreasing = FALSE)
        hit_table <- hit_table[o_hits,]
        rownames(hit_table) <- NULL
      } else {
        # create dummy table
        hit_table <- data.frame("idx1" = integer(),
                                "idx2" = integer(),
                                "inv" = integer(),
                                "width" = integer(),
                                "doc1_l" = integer(),
                                "doc1_r" = integer(),
                                "doc2_l" = integer(),
                                "doc2_r" = integer(),
                                "source_hit" = integer(),
                                "range1_hit" = integer())
      }
      # at this point using the range overlap functions feels logically obtuse
      # and it should be simpler to scroll through this object than previous
      # iterations, and hopefully i'm less of a dumbass than i was when i wrote
      # the original versions of this function
      
      if (nrow(hit_table) >= 1) {
        # set the hits table for the upper triangle
        # try to create the pairs table
        
        # at this point, hits are framed into the space of the doc1 genecalls
        # meaning that they could have been:
        # truncated
        # split apart and truncated
        # or left intact
        # because they've been split apart, I don't need to re-register the query ever
        # but i do need to reregister the subject search
        subject_search <- split(x = hit_table,
                                f = list(hit_table$idx1,
                                         hit_table$idx2),
                                drop = TRUE)
        # return(subject_search)
        prev_idx1 <- 0L
        prev_idx2 <- 0L
        modified_subject <- vector(mode = "list",
                                   length = length(subject_search))
        for (d3 in seq_along(subject_search)) {
          # overhead setting
          curr_idx_1 <- subject_search[[d3]]$idx1[1]
          curr_idx_2 <- subject_search[[d3]]$idx2[1]
          curr_subject <- subject_search[[d3]]
          # do i need this?
          # if (curr_idx_1 != prev_idx1) {
          #   
          # }
          if (curr_idx_2 != prev_idx2) {
            d2_w <- which(as.integer(names(d2_genecalls)) == curr_idx_2)
            curr_d2_genecalls <- d2_genecalls[[d2_w]]
            curr_d2_ranges <- d2_ranges[[d2_w]]
          }
          curr_focal_ranges <- IRanges(start = curr_subject$doc2_l,
                                       end = curr_subject$doc2_r)
          ph1 <- findOverlaps(query = curr_focal_ranges,
                              subject = curr_d2_ranges)
          if (length(ph1) < 1) {
            # next in the d3 loop, not the m2 or m1
            next
          }
          ph2 <- overlapsRanges(query = curr_focal_ranges,
                                subject = curr_d2_ranges,
                                hits = ph1)
          curr_from <- from(ph1)
          curr_to <- to(ph1)
          # hit idx and features are no long row based here, they need to be based
          # on the the pass through from the query search
          adj_vals <- data.frame("l" = start(ph2) - start(curr_focal_ranges[curr_from]),
                                 "r" = end(curr_focal_ranges[curr_from]) - end(ph2),
                                 "hit_idx" = curr_subject$source_hit[curr_from],
                                 "d1_feature" = curr_subject$range1_hit[curr_from])
          # adjustments here are framed forward forward in the doc2 space 
          # if inv == 1, then
          # in the doc1 space they get flipped (i think?)
          adj_tbl <- data.frame("idx1" = curr_subject$idx1[curr_from],
                                "idx2" = curr_subject$idx2[curr_from],
                                "inv" = curr_subject$inv[curr_from],
                                "width" = width(ph2),
                                "doc1_l" = rep(-1,
                                               length(curr_from)),
                                "doc1_r" = rep(-1,
                                               length(curr_from)),
                                "doc2_l" = curr_subject$doc2_l[curr_from] + adj_vals$l,
                                "doc2_r" = curr_subject$doc2_r[curr_from] - adj_vals$r,
                                "source_hit" = adj_vals$hit_idx,
                                "d1_hit" = adj_vals$d1_feature,
                                "d2_hit" = as.integer(rownames(curr_d2_genecalls)[curr_to]))
          # return(list(adj_tbl,
          #             d1_genecalls_source,
          #             d2_genecalls_source))
          # adj_tbl$d1_hit <- d1_genecalls_source$Key[adj_tbl$d1_hit]
          # adj_tbl$d2_hit <- d2_genecalls_source$Key[adj_tbl$d2_hit]
          inv_log1 <- adj_tbl$inv == 0
          inv_log2 <- adj_tbl$inv != 0
          adj_tbl$doc1_l[inv_log1] <- curr_subject$doc1_l[curr_from][inv_log1] + adj_vals$l[inv_log1]
          adj_tbl$doc1_l[inv_log2] <- curr_subject$doc1_l[curr_from][inv_log2] + adj_vals$r[inv_log2]
          adj_tbl$doc1_r[inv_log1] <- curr_subject$doc1_r[curr_from][inv_log1] - adj_vals$r[inv_log1]
          adj_tbl$doc1_r[inv_log2] <- curr_subject$doc1_r[curr_from][inv_log2] - adj_vals$l[inv_log2]
          
          modified_subject[[d3]] <- adj_tbl
          
          prev_idx1 <- curr_idx_1
          prev_idx2 <- curr_idx_2
        }
        # return(modified_subject)
        check_this <- vapply(X = modified_subject,
                             FUN = function(x) {
                               is.null(x)
                             },
                             FUN.VALUE = vector(mode = "logical",
                                                length = 1))
        modified_subject <- modified_subject[!check_this]
        
        # return(modified_subject)
        # manage our table possibilities
        # in any case where there is more than 0 rows in the final table, we need
        # to re-register the genecall ids from the sub-features back to the features
        if (length(modified_subject) == 1) {
          final_hit_table <- modified_subject[[1]]
          rownames(final_hit_table) <- NULL
        } else if (length(modified_subject) > 1) {
          final_hit_table <- do.call(rbind,
                                     modified_subject)
          rownames(final_hit_table) <- NULL
        } else {
          # fill in a dummy table
          final_hit_table <- data.frame("idx1" = integer(),
                                        "idx2" = integer(),
                                        "inv" = integer(),
                                        "width" = integer(),
                                        "doc1_l" = integer(),
                                        "doc1_r" = integer(),
                                        "doc2_l" = integer(),
                                        "doc2_r" = integer(),
                                        "source_hit" = integer(),
                                        "d1_hit" = integer(),
                                        "d2_hit" = integer(),
                                        "d1_subkey" = integer(),
                                        "d2_subkey" = integer())
        }
        
        # return(final_hit_table)
        if (nrow(final_hit_table) > 0) {
          # return(list(final_hit_table,
          #             d1_genecalls_source,
          #             d2_genecalls_source))
          # assign the subkey first because the hit vals are still the row indices
          # then assign the key value
          final_hit_table$d1_subkey <- d1_genecalls_source$SubKey[final_hit_table$d1_hit]
          final_hit_table$d1_hit <- d1_genecalls_source$Key[final_hit_table$d1_hit]
          final_hit_table$d2_subkey <- d2_genecalls_source$SubKey[final_hit_table$d2_hit]
          final_hit_table$d2_hit <- d2_genecalls_source$Key[final_hit_table$d2_hit]
          o_final <- order(final_hit_table$d1_hit,
                           final_hit_table$d2_hit,
                           final_hit_table$d1_subkey,
                           final_hit_table$d2_subkey)
          final_hit_table <- final_hit_table[o_final, ]
          
          # create the feature pairwise table
          pairwise_table <- split(x = final_hit_table,
                                  f = list(final_hit_table$d1_hit,
                                           final_hit_table$d2_hit),
                                  drop = TRUE)
          pairwise_table <- lapply(X = pairwise_table,
                                   FUN = function(x) {
                                     if (nrow(x) == 1) {
                                       # res <- data.frame("QueryGene" = x$d1_hit,
                                       #                   "SubjectGene" = x$d2_hit,
                                       #                   "ExactOverlap" = x$width,
                                       #                   "QueryIndex" = x$idx1,
                                       #                   "SubjectIndex" = x$idx2,
                                       #                   "QLeftPos" = x$doc1_l,
                                       #                   "QRightPos" = x$doc1_r,
                                       #                   "SLeftPos" = x$doc2_l,
                                       #                   "SRightPos" = x$doc2_r,
                                       #                   "MaxKmerSize" = x$width,
                                       #                   "TotalKmerHits" = 1)
                                       res <- c(x$d1_hit,
                                                x$d2_hit,
                                                x$width,
                                                x$idx1,
                                                x$idx2,
                                                x$doc1_l,
                                                x$doc1_r,
                                                x$doc2_l,
                                                x$doc2_r,
                                                x$width,
                                                1L)
                                     } else {
                                       # res <- data.frame("QueryGene" = x$d1_hit[1],
                                       #                   "SubjectGene" = x$d2_hit[1],
                                       #                   "ExactOverlap" = sum(x$width),
                                       #                   "QueryIndex" = x$idx1[1],
                                       #                   "SubjectIndex" = x$idx2[1],
                                       #                   "QLeftPos" = min(x$doc1_l),
                                       #                   "QRightPos" = max(x$doc1_r),
                                       #                   "SLeftPos" = min(x$doc2_l[1]),
                                       #                   "SRightPos" = max(x$doc2_r),
                                       #                   "MaxKmerSize" = sum(x$width),
                                       #                   "TotalKmerHits" = nrow(x))
                                       
                                       # not yet implemented, but need to figure out how
                                       # to correctly place bounds that recycle to the 
                                       # beginning of a contig
                                       # o_doc1 <- order(x$d1_hit,
                                       #                 x$d1_subkey)
                                       # o_doc2 <- order(x$d2_hit,
                                       #                 x$d2_subkey)
                                       res <- c(x$d1_hit[1],
                                                x$d2_hit[1],
                                                sum(x$width),
                                                x$idx1[1],
                                                x$idx2[1],
                                                min(x$doc1_l),
                                                max(x$doc1_r),
                                                min(x$doc2_l[1]),
                                                max(x$doc2_r),
                                                max(x$width),
                                                nrow(x))
                                     }
                                     return(res)
                                   })
          pairwise_table <- do.call(rbind,
                                    pairwise_table)
          rownames(pairwise_table) <- NULL
          o_pair <- order(pairwise_table[, 1],
                          pairwise_table[, 2],
                          pairwise_table[, 10],
                          pairwise_table[, 11])
          pairwise_table <- pairwise_table[o_pair, ]
        } else {
          # create dummy pairs table
          pairwise_table <- matrix(data = integer(),
                                   ncol = 10L)
        } # end if else on the final hit table
        
        # there needs to be a hits table for the lower triangle
        # and a pairs table for the upper triangle
        
      } else {
        # create a dummy hits table and a dummy pairs table
        final_hit_table <- data.frame("idx1" = integer(),
                                      "idx2" = integer(),
                                      "inv" = integer(),
                                      "width" = integer(),
                                      "doc1_l" = integer(),
                                      "doc1_r" = integer(),
                                      "doc2_l" = integer(),
                                      "doc2_r" = integer(),
                                      "source_hit" = integer(),
                                      "d1_hit" = integer(),
                                      "d2_hit" = integer(),
                                      "d1_subkey" = integer(),
                                      "d2_subkey" = integer())
        pairwise_table <- matrix(data = integer(),
                                 ncol = 10L)
      } # end if else on hit table rows
      
      colnames(pairwise_table) <- c("QueryGene",
                                    "SubjectGene",
                                    "ExactOverlap",
                                    "QueryIndex",
                                    "SubjectIndex",
                                    "QLeftPos",
                                    "QRightPos",
                                    "SLeftPos",
                                    "SRightPos",
                                    "MaxKmerSize",
                                    "TotalKmerHits")
      final_hit_table <- cbind("QueryGene" = final_hit_table$d1_hit,
                               "SubjectGene" = final_hit_table$d2_hit,
                               "ExactOverlap" = final_hit_table$width,
                               "QueryIndex" = final_hit_table$idx1,
                               "SubjectIndex" = final_hit_table$idx2,
                               "QLeftPos" = final_hit_table$doc1_l,
                               "QRightPos" = final_hit_table$doc1_r,
                               "SLeftPos" = final_hit_table$doc2_l,
                               "SRightPos" = final_hit_table$doc2_r,
                               "QuerySubKey" = final_hit_table$d1_subkey,
                               "SubjectSubKey" = final_hit_table$d2_subkey,
                               "SyntenicOrigin" = final_hit_table$source_hit)
      
      
      ResultMatrix[[m1, m2]] <- pairwise_table
      ResultMatrix[[m2, m1]] <- final_hit_table
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = TotalCounter / TotalLength)
        TotalCounter <- TotalCounter + 1L
      }
    } # end of columns loop
  } # end of rows loop
  if (Verbose) {
    TotalTimeStop <- Sys.time()
    close(pBar)
    print(TotalTimeStop - TotalTimeStart)
  }
  dimnames(ResultMatrix) <- dimnames(SyntenyObject)
  class(ResultMatrix) <- "LinkedPairs"
  ph1 <- lapply(X = FeatureRepresentations,
                FUN = function(x) {
                  attr(x = x,
                       which = "origin")
                })
  names(ph1) <- names(FeatureRepresentations)
  attr(ResultMatrix, "GeneCalls") <- ph1
  return(ResultMatrix)
}


