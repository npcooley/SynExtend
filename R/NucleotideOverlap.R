# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu

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
  if (L != length(GeneCalls)) {
    stop ("Synteny Object and Gene Calls have differing number of genomes.")
  }
  if (is.null(names(GeneCalls))) {
    stop ("Gene Predictions must be named.")
  }
  if (is.null(rownames(SyntenyObject))) {
    stop ("Synteny Object must have named identifiers.")
  }
  if (any(names(GeneCalls) != rownames(SyntenyObject))) {
    stop ("Names between Synteny Object and Gene Predictions do not match.")
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
  TotalLength <- L^2 - L
  TotalCounter <- 0L
  ######
  # scroll through every hit table in the synteny object
  ######
  
  ###### -- Deal with Different GeneCall types --------------------------------
  
  SynNames <- unname(sapply(diag(SyntenyObject),
                            function(x) gsub(x = names(x),
                                             pattern = " .+",
                                             replacement = ""),
                            simplify = FALSE,
                            USE.NAMES = FALSE))
  FeatureRepresentations <- vector(mode = "list",
                                   length = L)
  names(FeatureRepresentations) <- names(GeneCalls)
  ContigNames <- vector(mode = "list",
                        length = L)
  
  GCallClasses <- sapply(GeneCalls,
                         function(x) class(x),
                         USE.NAMES = FALSE,
                         simplify = TRUE)
  if (any(GCallClasses == "GRanges")) {
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
                 StartConversion)
      
      StrandConversion <- StrandConversion[o]
      IndexConversion <- C.Index[o]
      StartConversion <- StartConversion[o]
      StopConversion <- StopConversion[o]
      LengthsConversion <- LengthsConversion[o]
      
      # FeatureRepresentations[[m1]] <- matrix(data = c(IndexConversion,
      #                                                 StrandConversion,
      #                                                 StartConversion,
      #                                                 StopConversion,
      #                                                 LengthsConversion),
      #                                        nrow = length(o),
      #                                        ncol = 5L)
      
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
      # CurrentLengths <- CurrentStops - CurrentStarts + 1L
      
      # FeatureRepresentations[[m1]] <- matrix(data = c(CurrentIndices,
      #                                                 CurrentStrands,
      #                                                 CurrentStarts,
      #                                                 CurrentStops,
      #                                                 CurrentLengths),
      #                                        nrow = length(o),
      #                                        ncol = 5L)
      FeatureRepresentations[[m1]] <- GeneCalls[[m1]][o, ]
      # rm(list = c("CurrentIndices",
      #             "CurrentStrands",
      #             "CurrentStarts",
      #             "CurrentStops",
      #             "CurrentLengths"))
      
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
                                       pattern = " .+",
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
                                               length(o)))
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
  
  ###### -- End Gene call stuff -----------------------------------------------
  
  pBar <- txtProgressBar(style = 1L)
  for (m1 in seq_len(L - 1L)) {
    for (m2 in (m1 +1L):L) {
      ######
      # Collect Index Start Stop and Strand from current subject and query gene calls
      # Collect Index Start Stop Strand from hit table
      ######
      
      Q.Index <- FeatureRepresentations[[m1]][, 1L]
      Q.Start <- FeatureRepresentations[[m1]][, 3L]
      Q.Stop <- FeatureRepresentations[[m1]][, 4L]
      QG.Strand <- FeatureRepresentations[[m1]][, 2L]
      
      S.Index <- FeatureRepresentations[[m2]][, 1L]
      S.Start <- FeatureRepresentations[[m2]][, 3L]
      S.Stop <- FeatureRepresentations[[m2]][, 4L]
      SG.Strand <- FeatureRepresentations[[m2]][, 2L]
      
      # if (AcceptContigNames) {
      #   Q.Index <- FeatureRepresentations[[m1]][, 1L]
      #   ph1 <- unname(ContigNames[[m1]])
      #   ph2 <- match(x = names(ContigNames[[m1]]),
      #                table = SynNames[[m1]])
      #   ph3 <- vector(mode = "list",
      #                 length = length(ph1))
      #   for (m3 in seq_along(ph1)) {
      #     ph3[[m3]] <- which(Q.Index == ph1[m3])
      #   }
      #   for (m3 in seq_along(ph1)) {
      #     if (length(ph3) > 0L) {
      #       Q.Index <- replace(x = Q.Index,
      #                          list = ph3[[m3]],
      #                          values = ph2[m3])
      #     }
      #   }
      #   rm(list = c("ph1",
      #               "ph2",
      #               "ph3"))
      #   S.Index <- FeatureRepresentations[[m2]][, 1L]
      #   ph1 <- unname(ContigNames[[m2]])
      #   ph2 <- match(x = names(ContigNames[[m2]]),
      #                table = SynNames[[m2]])
      #   ph3 <- vector(mode = "list",
      #                 length = length(ph1))
      #   for (m3 in seq_along(ph1)) {
      #     ph3[[m3]] <- which(S.Index == ph1[m3])
      #   }
      #   for (m3 in seq_along(ph2)) {
      #     if (length(ph3[[m3]]) > 0L) {
      #       S.Index <- replace(x = S.Index,
      #                          list = ph3[[m3]],
      #                          values = ph2[m3])
      #     }
      #   }
      #   rm(list = c("ph1",
      #               "ph2",
      #               "ph3"))
      #   o <- order(S.Index,
      #              S.Start)
      #   S.Start <- S.Start[o]
      #   S.Stop <- S.Stop[o]
      #   SG.Strand <- SG.Strand[o]
      #   S.Index <- S.Index[o]
      #   o <- order(Q.Index,
      #              Q.Start)
      #   Q.Start <- Q.Start[o]
      #   Q.Stop <- Q.Stop[o]
      #   QG.Strand <- QG.Strand[o]
      #   Q.Index <- Q.Index[o]
      #   rm(list = c("o"))
      # } else {
      #   Q.Index <- FeatureRepresentations[[m1]][, 1L]
      #   S.Index <- FeatureRepresentations[[m2]][, 1L]
      # }
      
      CurrentHitTable <- SyntenyObject[m1, m2, drop = FALSE][[1]]
      # return(list(unique(CurrentHitTable[, c(1,2)]),
      #             FeatureRepresentations[[m1]],
      #             FeatureRepresentations[[m2]],
      #             LimitIndex,
      #             diag(ResultMatrix),
      #             SynNames,
      #             FeatureRepresentations,
      #             Q.Index,
      #             FeatureRepresentations[[m1]][, 1L],
      #             S.Index,
      #             FeatureRepresentations[[m2]][, 1L]))
      if (LimitIndex) {
        ######
        # Select current hit table, subset to only single index
        ######
        CurrentHitTable <- CurrentHitTable[which(CurrentHitTable[, "index1"] == 1L), ]
        CurrentHitTable <- CurrentHitTable[which(CurrentHitTable[, "index2"] == 1L), ]
        CurrentHitTable <- CurrentHitTable[order(CurrentHitTable[,
                                                                 "start1",
                                                                 drop = FALSE]),
                                           ,
                                           drop = FALSE]
      } else {
        ######
        # If not limiting by index, order hit table by the index1 first, then the query starts
        ######
        CurrentHitTable <- CurrentHitTable[order(CurrentHitTable[,
                                                                 "index1",
                                                                 drop = FALSE],
                                                 CurrentHitTable[,
                                                                 "start1",
                                                                 drop = FALSE]),
                                           ,
                                           drop = FALSE]
      }
      ######
      # Collect the starts and stops for all hits correct for strandedness
      # Collect indices as well
      ######
      HitWidths <- CurrentHitTable[, "width"]
      Q.HitStarts <- CurrentHitTable[, "start1"]
      S.HitStarts <- CurrentHitTable[, "start2"]
      Q.HitEnds <- Q.HitStarts + HitWidths - 1L
      S.HitEnds <- S.HitStarts + HitWidths - 1L
      Strand <- CurrentHitTable[, "strand"]
      QHI <- CurrentHitTable[, "index1"]
      SHI <- CurrentHitTable[, "index2"]
      for (i in seq_along(HitWidths)) {
        if (Strand[i] == 0L) {
          S.HitStarts[i] <- CurrentHitTable[i, "start2"]
          S.HitEnds[i] <- CurrentHitTable[i, "start2"] + CurrentHitTable[i, "width"] - 1L
        } else if (Strand[i] == 1L) {
          S.HitStarts[i] <- CurrentHitTable[i, "start2"] - CurrentHitTable[i, "width"] + 1L
          S.HitEnds[i] <- CurrentHitTable[i, "start2"]
        }
      }
      ######
      # Record a query matrix for every hit that lands inside a gene in our query geneset/genome
      # contains:
      # the gene index and hit index in question
      # modifies the hit boundaries IN THE SUBJECT
      # records the distances between the hit boundaries and the gene boundaries
      # the strand for the gene in the query that is being recorded
      ######
      QueryMatrix <- matrix(NA_integer_,
                            nrow = 8L,
                            ncol = nrow(CurrentHitTable))
      # return(list(CurrentHitTable,
      #             FeatureRepresentations,
      #             Q.Index,
      #             S.Index))
      HitCounter <- 1L
      AddCounter <- 1L
      DimLimit <- nrow(CurrentHitTable) / 2
      DimAdjust <- 2L
      
      # return(list(Q.Start,
      #             Q.Stop,
      #             QG.Strand,
      #             Q.Index,
      #             S.Start,
      #             S.Stop,
      #             SG.Strand,
      #             S.Index,
      #             SHI))
      
      for (z1 in seq_along(Q.Start)) {
        ######
        # loop through the starts of the query genes as they have been ordered
        ######
        CurrentGene <- HitCounter
        Q.NucOverLapL <- NA_integer_
        Q.NucOverLapR <- NA_integer_
        S.NucPositionL <- NA_integer_
        S.NucPositionR <- NA_integer_
        S.Strand <- NA_integer_
        
        ######
        # re-register the hit search if a new gene has a start position that
        # occurs before the end of the previous gene
        # re-registers can only occur when
        ######
        if (z1 > 1L &
            HitCounter > 1L &
            HitCounter <= length(HitWidths)) {
          # not the first gene
          # not the first hit
          if (Q.Start[z1] < Q.Stop[z1 - 1L]) {
            # start of the current gene occurs before end of previous gene
            # print(paste(QHI[HitCounter],
            #             Q.Index[z1],
            #             HitCounter,
            #             z1))
            if (QHI[HitCounter] == Q.Index[z1] &
                Q.Index[z1] == Q.Index[z1 - 1L]) {
              # gene and hit are in the same index
              # gene and previous gene are in the same index
              AvailableReRegisters <- QHI == Q.Index[z1] & Q.HitStarts < Q.Start[z1]
              if (sum(AvailableReRegisters) > 0L) {
                # more than zero re-register positions available
                TempHit <- max(which(AvailableReRegisters))
                if (TempHit < HitCounter) {
                  HitCounter <- TempHit
                }
              } # else do nothing, no hits that start before the current start and occupy the correct index
            } # else do nothing - indices are not matched
          } # else do nothing - start of current gene is after stop of previous
        } # else do nothing - either first gene or hit
        
        while (HitCounter <= length(HitWidths)) {
          ######
          # interrogation loop, depending upon where the current hit is
          # and where the current gene is
          # either go to the next hit (HitCounter + 1L)
          # or
          # go to the next gene (break)
          # AddCounter tells the loop which row to add information to in
          # the initialized QueryMatrix, when there is information to add
          ######
          
          if (Q.HitEnds[HitCounter] < Q.Start[z1] &
              QHI[HitCounter] == Q.Index[z1]) {
            # Hit ends before current query gene begins
            # AND the indices are the same
            # record the hit counter every time you iterate to the next hit
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] < Q.Start[z1] &
                     Q.HitEnds[HitCounter] >= Q.Start[z1] &
                     Q.HitEnds[HitCounter] <= Q.Stop[z1] &
                     QHI[HitCounter] == Q.Index[z1]) {
            # Hit overlaps left bound of current query gene
            # Stay on the current gene
            # And go to the next hit
            CurrentGene <- z1
            Q.NucOverLapL <- Q.Start[z1]
            Q.NucOverLapR <- Q.HitEnds[HitCounter]
            S.Strand <- Strand[HitCounter]
            NewWidth <- Q.NucOverLapR - Q.NucOverLapL + 1L
            Current.QHI <- QHI[HitCounter]
            Current.SHI <- SHI[HitCounter]
            ######
            # Trim hit in subject, overlaps that are in the reverse strand require different trims
            ######
            if (S.Strand == 0L) {
              S.NucPositionR <- S.HitEnds[HitCounter]
              S.NucPositionL <- S.HitEnds[HitCounter] - NewWidth + 1L
            } else {
              S.NucPositionL <- S.HitStarts[HitCounter]
              S.NucPositionR <- S.HitStarts[HitCounter] + NewWidth - 1L
            }
            
            # Add To Vector!
            QueryMatrix[, AddCounter] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR,
                                           S.Strand,
                                           Current.QHI,
                                           Current.SHI)
            if (AddCounter >= DimLimit) {
              QueryMatrix <- cbind(QueryMatrix,
                                   matrix(data = NA_integer_,
                                          nrow = nrow(QueryMatrix),
                                          ncol = ncol(QueryMatrix) * DimAdjust))
              DimLimit <- ncol(QueryMatrix) / 2
              DimAdjust <- DimAdjust * 2L
            }
            # record the add counter every time a new row is added to the query matrix
            AddCounter <- AddCounter + 1L
            # record the hit counter every time you iterate to the next hit
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] >= Q.Start[z1] &
                     Q.HitEnds[HitCounter] <= Q.Stop[z1] &
                     QHI[HitCounter] == Q.Index[z1]) {
            # Hit occurs entirely within current query gene
            # Stay on the current gene
            # And go to the next hit
            CurrentGene <- z1
            Q.NucOverLapL <- Q.HitStarts[HitCounter]
            Q.NucOverLapR <- Q.HitEnds[HitCounter]
            S.Strand <- Strand[HitCounter]
            S.NucPositionL <- S.HitStarts[HitCounter]
            S.NucPositionR <- S.HitEnds[HitCounter]
            Current.QHI <- QHI[HitCounter]
            Current.SHI <- SHI[HitCounter]
            # Add To Vector!
            QueryMatrix[, AddCounter] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR,
                                           S.Strand,
                                           Current.QHI,
                                           Current.SHI)
            if (AddCounter >= DimLimit) {
              QueryMatrix <- cbind(QueryMatrix,
                                   matrix(data = NA_integer_,
                                          nrow = nrow(QueryMatrix),
                                          ncol = ncol(QueryMatrix) * DimAdjust))
              DimLimit <- ncol(QueryMatrix) / 2
              DimAdjust <- DimAdjust * 2L
            }
            # record the add counter every time a new row is added to the query matrix
            AddCounter <- AddCounter + 1L
            # record the hit counter every time you iterate to the next hit
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] >= Q.Start[z1] &
                     Q.HitStarts[HitCounter] <= Q.Stop[z1] &
                     Q.HitEnds[HitCounter] > Q.Stop[z1] &
                     QHI[HitCounter] == Q.Index[z1]) {
            # Hit overlaps right bound of current query gene
            # Stay on the current hit and go to the next gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.HitStarts[HitCounter]
            Q.NucOverLapR <- Q.Stop[z1]
            S.Strand <- Strand[HitCounter]
            NewWidth <- Q.NucOverLapR - Q.NucOverLapL + 1L
            Current.QHI <- QHI[HitCounter]
            Current.SHI <- SHI[HitCounter]
            if (S.Strand == 0L) {
              S.NucPositionL <- S.HitStarts[HitCounter]
              S.NucPositionR <- S.HitStarts[HitCounter] + NewWidth - 1L
            } else {
              S.NucPositionR <- S.HitEnds[HitCounter]
              S.NucPositionL <- S.HitEnds[HitCounter] - NewWidth + 1L
            }
            # Add To Vector!
            QueryMatrix[, AddCounter] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR,
                                           S.Strand,
                                           Current.QHI,
                                           Current.SHI)
            if (AddCounter >= DimLimit) {
              QueryMatrix <- cbind(QueryMatrix,
                                   matrix(data = NA_integer_,
                                          nrow = nrow(QueryMatrix),
                                          ncol = ncol(QueryMatrix) * DimAdjust))
              DimLimit <- ncol(QueryMatrix) / 2
              DimAdjust <- DimAdjust * 2L
            }
            # record the add counter every time a new row is added to the query matrix
            AddCounter <- AddCounter + 1L
            break
          } else if (Q.HitStarts[HitCounter] < Q.Start[z1] &
                     Q.HitEnds[HitCounter] > Q.Stop[z1] &
                     QHI[HitCounter] == Q.Index[z1]) {
            # Hit eclipses current query gene
            # Stay on the current hit, and go to the next gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.Start[z1]
            Q.NucOverLapR <- Q.Stop[z1]
            TrimLeft <- Q.Start[z1] - Q.HitStarts[HitCounter]
            TrimRight <- Q.HitEnds[HitCounter] - Q.Stop[z1]
            S.Strand <- Strand[HitCounter]
            NewWidth <- Q.NucOverLapR - Q.NucOverLapL + 1L
            Current.QHI <- QHI[HitCounter]
            Current.SHI <- SHI[HitCounter]
            if (S.Strand == 0L) {
              S.NucPositionL <- S.HitStarts[HitCounter] + TrimLeft
              S.NucPositionR <- S.NucPositionL + NewWidth - 1L
            } else {
              S.NucPositionR <- S.HitEnds[HitCounter] - TrimLeft
              S.NucPositionL <- S.NucPositionR - NewWidth + 1L
            }
            # Add To Vector!
            QueryMatrix[, AddCounter] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR,
                                           S.Strand,
                                           Current.QHI,
                                           Current.SHI)
            if (AddCounter >= DimLimit) {
              QueryMatrix <- cbind(QueryMatrix,
                                   matrix(data = NA_integer_,
                                          nrow = nrow(QueryMatrix),
                                          ncol = ncol(QueryMatrix) * DimAdjust))
              DimLimit <- ncol(QueryMatrix) / 2
              DimAdjust <- DimAdjust * 2L
            }
            # record the add counter every time a new row is added to the query matrix
            AddCounter <- AddCounter + 1L
            break
          } else if (Q.HitStarts[HitCounter] > Q.Stop[z1] &
                     QHI[HitCounter] == Q.Index[z1]) {
            # Hit occurs after current query gene
            # And the indices are matched
            # Go to next gene
            break
          } else if (QHI[HitCounter] > Q.Index[z1]) {
            # If the index of the hits in the query has outpaced the index of the genes
            break
          } else if (QHI[HitCounter] < Q.Index[z1]) {
            # If the index of the genes has outpaced the index of the hits
            # record the hit counter every time you iterate to the next hit
            HitCounter <- HitCounter + 1L
          } # end of else if conditionals
        } # end while loop through hits
      } # end for loop through genes
      if (Verbose) {
        TotalCounter <- TotalCounter + 1L
        setTxtProgressBar(pb = pBar,
                          value = TotalCounter/TotalLength)
      }
      QueryMatrix <- t(QueryMatrix)
      colnames(QueryMatrix) <- c("CurrentGene", # The gene that was recorded with a hit present
                                 "QueryNucleotideOverLapLeft", # leftmost nucleotide position of a hit that is within CurrentGene
                                 "QueryNucleotideOverLapRight", # rightmost nucleotide position of a hit that is within CurrentGene
                                 "SubjectNucleotidePositionLeft", # corresponding leftmost position in the subject
                                 "SubjectNucleotidePositionRight", # corresponding rightmost position in the subject
                                 "SubjectStrand", # strand of the hit in the subject
                                 "QueryIndex", # index that was matched between the hit and the gene in the query
                                 "SubjectIndex" # index (chromosome, plasmid, whatever) of the hit in the subject
                                 )
      # return(QueryMatrix)
      ######
      # Remove unfilled extra rows
      ######
      QueryMatrix <- QueryMatrix[apply(QueryMatrix,
                                       1L,
                                       function(x) !all(is.na(x))),
                                 ,
                                 drop = FALSE]
      if (dim(QueryMatrix)[1] == 0L) {
        OutPutMatrix <- matrix(NA_integer_,
                               nrow = 1L,
                               ncol = 7L)
        OverLapMatrix <- matrix(NA_integer_,
                                nrow = 1L,
                                ncol = 7L)
      } else {
        HitCounter <- 1L
        AddCounter <- 1L
        OverLapMatrix <- matrix(NA_integer_,
                                ncol = nrow(QueryMatrix),
                                nrow = 9L)
        QueryMatrix <- QueryMatrix[order(QueryMatrix[,
                                                     "SubjectIndex",
                                                     drop = FALSE],
                                         QueryMatrix[,
                                                     "SubjectNucleotidePositionLeft",
                                                     drop = FALSE]),
                                   ,
                                   drop = FALSE]
        DimAdjust <- 2L
        DimLimit <- ncol(OverLapMatrix) / 2
        ######
        # Part 2!
        # Hits that were recorded as being within a gene in the query
        # are now tested again the genes in the subject
        ######
        QueryMap <- QueryMatrix[, "CurrentGene"]
        S.HitStarts <- QueryMatrix[, "SubjectNucleotidePositionLeft"]
        S.HitEnds <- QueryMatrix[, "SubjectNucleotidePositionRight"]
        Q.HitStarts <- QueryMatrix[, "QueryNucleotideOverLapLeft"]
        Q.HitEnds <- QueryMatrix[, "QueryNucleotideOverLapRight"]
        QueryIndices <- QueryMatrix[, "QueryIndex"]
        SubjectHitIndex <- QueryMatrix[, "SubjectIndex"]
        
        # return(list(QueryMatrix,
        #             S.Start,
        #             S.Stop,
        #             S.Index,
        #             FeatureRepresentations[[m1]],
        #             FeatureRepresentations[[m2]]))
        for (z2 in seq_along(S.Start)) {
          ######
          # Loop through the subject genes
          ######
          while (HitCounter <= length(QueryMap)) {
            if (S.HitEnds[HitCounter] < S.Start[z2] &
                SubjectHitIndex[HitCounter] == S.Index[z2]) {
              # Hit ends before current subject gene begins
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] < S.Start[z2] &
                       S.HitEnds[HitCounter] >= S.Start[z2] &
                       S.HitEnds[HitCounter] <= S.Stop[z2] &
                       SubjectHitIndex[HitCounter] == S.Index[z2]) {
              # Hit overlaps left bound of current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              CurrentQueryIndex <- QueryIndices[HitCounter]
              CurrentSubjectIndex <- SubjectHitIndex[HitCounter]
              TrimLeft <- S.Start[z2] - S.HitStarts[HitCounter]
              ExactOverLap <- S.HitEnds[HitCounter] - S.Start[z2] + 1L
              SubjectHitLeft <- S.HitStarts[HitCounter] + TrimLeft
              SubjectHitRight <- S.HitEnds[HitCounter]
              QueryHitLeft <- Q.HitStarts[HitCounter] + TrimLeft
              QueryHitRight <- Q.HitEnds[HitCounter]
              # Add to vector !
              OverLapMatrix[, AddCounter] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap,
                                               CurrentQueryIndex,
                                               CurrentSubjectIndex,
                                               QueryHitLeft,
                                               QueryHitRight,
                                               SubjectHitLeft,
                                               SubjectHitRight)
              if (AddCounter >= DimLimit) {
                OverLapMatrix <- cbind(OverLapMatrix,
                                       matrix(data = NA_integer_,
                                              nrow = nrow(OverLapMatrix),
                                              ncol = ncol(OverLapMatrix) * DimAdjust))
                DimLimit <- ncol(OverLapMatrix) / 2
                DimAdjust <- DimAdjust * 2L
              }
              AddCounter <- AddCounter + 1L
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] >= S.Start[z2] &
                       S.HitEnds[HitCounter] <= S.Stop[z2] &
                       SubjectHitIndex[HitCounter] == S.Index[z2]) {
              # Hit occurs entirely within current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.HitEnds[HitCounter] - S.HitStarts[HitCounter] + 1L
              CurrentQueryIndex <- QueryIndices[HitCounter]
              CurrentSubjectIndex <- SubjectHitIndex[HitCounter]
              SubjectHitLeft <- S.HitStarts[HitCounter]
              SubjectHitRight <- S.HitEnds[HitCounter]
              QueryHitLeft <- Q.HitStarts[HitCounter]
              QueryHitRight <- Q.HitEnds[HitCounter]
              # Add to vector !
              OverLapMatrix[, AddCounter] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap,
                                               CurrentQueryIndex,
                                               CurrentSubjectIndex,
                                               QueryHitLeft,
                                               QueryHitRight,
                                               SubjectHitLeft,
                                               SubjectHitRight)
              if (AddCounter >= DimLimit) {
                OverLapMatrix <- cbind(OverLapMatrix,
                                       matrix(data = NA_integer_,
                                              nrow = nrow(OverLapMatrix),
                                              ncol = ncol(OverLapMatrix) * DimAdjust))
                DimLimit <- ncol(OverLapMatrix) / 2
                DimAdjust <- DimAdjust * 2L
              }
              AddCounter <- AddCounter + 1L
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] >= S.Start[z2] &
                       S.HitStarts[HitCounter] <= S.Stop[z2] &
                       S.HitEnds[HitCounter] > S.Stop[z2] &
                       SubjectHitIndex[HitCounter] == S.Index[z2]) {
              # Hit overlaps right bound of current subject gene
              # Stay on the current hit and go to the next gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.Stop[z2] - S.HitStarts[HitCounter] + 1L
              CurrentQueryIndex <- QueryIndices[HitCounter]
              CurrentSubjectIndex <- SubjectHitIndex[HitCounter]
              TrimRight <- S.HitEnds[HitCounter] - S.Stop[z2]
              SubjectHitLeft <- S.HitStarts[HitCounter]
              SubjectHitRight <- S.HitEnds[HitCounter] - TrimRight
              QueryHitLeft <- Q.HitStarts[HitCounter]
              QueryHitRight <- Q.HitEnds[HitCounter] - TrimRight
              # Add to vector !
              OverLapMatrix[, AddCounter] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap,
                                               CurrentQueryIndex,
                                               CurrentSubjectIndex,
                                               QueryHitLeft,
                                               QueryHitRight,
                                               SubjectHitLeft,
                                               SubjectHitRight)
              if (AddCounter >= DimLimit) {
                OverLapMatrix <- cbind(OverLapMatrix,
                                       matrix(data = NA_integer_,
                                              nrow = nrow(OverLapMatrix),
                                              ncol = ncol(OverLapMatrix) * DimAdjust))
                DimLimit <- ncol(OverLapMatrix) / 2
                DimAdjust <- DimAdjust * 2L
              }
              AddCounter <- AddCounter + 1L
              break
            } else if (S.HitStarts[HitCounter] <= S.Start[z2] &
                       S.HitEnds[HitCounter] >= S.Stop[z2] &
                       SubjectHitIndex[HitCounter] == S.Index[z2]) {
              # Hit eclipses current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.Stop[z2] - S.Start[z2] + 1L
              CurrentQueryIndex <- QueryIndices[HitCounter]
              CurrentSubjectIndex <- SubjectHitIndex[HitCounter]
              TrimLeft <- S.Start[z2] - S.HitStarts[HitCounter]
              TrimRight <- S.HitEnds[HitCounter] - S.Stop[z2]
              SubjectHitLeft <- S.HitStarts[HitCounter] + TrimLeft
              SubjectHitRight <- S.HitEnds[HitCounter] - TrimRight
              QueryHitLeft <- Q.HitStarts[HitCounter] + TrimLeft
              QueryHitRight <- Q.HitEnds[HitCounter] - TrimRight
              # Add to vector !
              OverLapMatrix[, AddCounter] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap,
                                               CurrentQueryIndex,
                                               CurrentSubjectIndex,
                                               QueryHitLeft,
                                               QueryHitRight,
                                               SubjectHitLeft,
                                               SubjectHitRight)
              if (AddCounter >= DimLimit) {
                OverLapMatrix <- cbind(OverLapMatrix,
                                       matrix(data = NA_integer_,
                                              nrow = nrow(OverLapMatrix),
                                              ncol = ncol(OverLapMatrix) * DimAdjust))
                DimLimit <- ncol(OverLapMatrix) / 2
                DimAdjust <- DimAdjust * 2L
              }
              AddCounter <- AddCounter + 1L
              break
            } else if (S.HitStarts[HitCounter] > S.Stop[z2]) {
              # Hit occurs after current subject gene
              break
            } else if (SubjectHitIndex[HitCounter] < S.Index[z2]) {
              # The indices do not match the only direction to go is forward
              HitCounter <- HitCounter + 1L
            } else if (SubjectHitIndex[HitCounter] > S.Index[z2]) {
              # if you have outpaced the gene index with the hits go to the next gene
              break
            }
          }
        }
      }
      
      OverLapMatrix <- t(OverLapMatrix)
      ######
      # Remove empty extra rows
      ######
      OverLapMatrix <- OverLapMatrix[apply(OverLapMatrix,
                                           1L,
                                           function(x) !all(is.na(x))),
                                     ,
                                     drop = FALSE]
      ######
      # If the overlap matrix is empty, assign an empty single row matrix
      # If it is not, sum up the overlap in nucleotide space and condense to a single row
      # Select the smallest distances for all four displacements
      ######
      if (dim(OverLapMatrix)[1] == 0L) {
        OverLapMatrix <- matrix(NA_integer_,
                                nrow = 1L,
                                ncol = 11L)
        OutPutMatrix <- matrix(NA_integer_,
                               nrow = 1L,
                               ncol = 11L)
      } else if (dim(OverLapMatrix)[1] >= 1L) {
        OverLapMatrix <- OverLapMatrix[order(OverLapMatrix[, 1L, drop = FALSE],
                                             OverLapMatrix[, 2L, drop = FALSE],
                                             OverLapMatrix[, 4L, drop = FALSE],
                                             OverLapMatrix[, 5L, drop = FALSE]),
                                       ,
                                       drop = FALSE]
        OutPutMatrix <- matrix(NA_integer_,
                               ncol = 11L,
                               nrow = nrow(OverLapMatrix))
        ######
        # All recordings so far are by hit
        # condense hits as appropriate, by the genes that they link
        ######
        RowCount <- 1L
        CondenseCount <- 1L
        Row <- 2L
        while (CondenseCount <= nrow(OverLapMatrix)) {
          while (Row <= nrow(OverLapMatrix)) {
            if (OverLapMatrix[Row, 5L] != OverLapMatrix[CondenseCount, 5L]) {
              break
            }
            if (OverLapMatrix[Row, 4L] != OverLapMatrix[CondenseCount, 4L]) {
              break
            }
            if (OverLapMatrix[Row, 2L] != OverLapMatrix[CondenseCount, 2L]) {
              break
            }
            if (OverLapMatrix[Row, 1L] != OverLapMatrix[CondenseCount, 1L]) {
              break
            }
            Row <- Row + 1L
          }
          
          z5 <- CondenseCount:(Row - 1L)
          OutPutMatrix[RowCount, ] <- c(OverLapMatrix[CondenseCount, 1L],
                                        OverLapMatrix[CondenseCount, 2L],
                                        sum(OverLapMatrix[z5, 3L]),
                                        OverLapMatrix[CondenseCount, 4L],
                                        OverLapMatrix[CondenseCount, 5L],
                                        min(OverLapMatrix[z5, 6L]),
                                        max(OverLapMatrix[z5, 7L]),
                                        min(OverLapMatrix[z5, 8L]),
                                        max(OverLapMatrix[z5, 9L]),
                                        max(OverLapMatrix[z5, 3L]),
                                        nrow(OverLapMatrix[z5, , drop = FALSE]))
          RowCount <- RowCount + 1L
          CondenseCount <- CondenseCount + length(z5)
        }
        OutPutMatrix <- OutPutMatrix[apply(OutPutMatrix,
                                           1L,
                                           function(x) !all(is.na(x))),
                                     ,
                                     drop = FALSE]
      }
      OutPutMatrix <- OutPutMatrix[order(OutPutMatrix[,
                                                      1L,
                                                      drop = FALSE]),
                                   ,
                                   drop = FALSE]
      
      colnames(OutPutMatrix) <- c("QueryGene",
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
      colnames(OverLapMatrix) <- c("QueryGene",
                                   "SubjectGene",
                                   "ExactOverlap",
                                   "QueryIndex",
                                   "SubjectIndex",
                                   "QLeftPos",
                                   "QRightPos",
                                   "SLeftPos",
                                   "SRightPos")
      # return(list(OverLapMatrix,
      #             OutPutMatrix))
      # QueryStartDisplacement <- ifelse(test = QG.Strand[OutPutMatrix[, "QueryGene"]] == 1L,
      #                                  yes = abs(OutPutMatrix[, 7L] - Q.Stop[OutPutMatrix[, 1L]]),
      #                                  no = abs(OutPutMatrix[, 6L] - Q.Start[OutPutMatrix[, 1L]]))
      # QueryStopDisplacement <- ifelse(test = QG.Strand[OutPutMatrix[, "QueryGene"]] == 1L,
      #                                 yes = abs(OutPutMatrix[, 6L] - Q.Start[OutPutMatrix[, 1L]]),
      #                                 no = abs(OutPutMatrix[, 7L] - Q.Stop[OutPutMatrix[, 1L]]))
      # SubjectStartDisplacement <- ifelse(test = SG.Strand[OutPutMatrix[, "SubjectGene"]] == 1L,
      #                                    yes = abs(OutPutMatrix[, 9L] - S.Stop[OutPutMatrix[, 2L]]),
      #                                    no = abs(OutPutMatrix[, 8L] - S.Start[OutPutMatrix[, 2L]]))
      # SubjectStopDisplacement <- ifelse(test = SG.Strand[OutPutMatrix[, "SubjectGene"]] == 1L,
      #                                   yes = abs(OutPutMatrix[, 8L] - S.Start[OutPutMatrix[, 2L]]),
      #                                   no = abs(OutPutMatrix[, 9L] - S.Stop[OutPutMatrix[, 2L]]))
      # DisplacementMatrix <- cbind(QueryStartDisplacement,
      #                             QueryStopDisplacement,
      #                             SubjectStartDisplacement,
      #                             SubjectStopDisplacement)
      
      if (Verbose) {
        TotalCounter <- TotalCounter + 1L
        setTxtProgressBar(pb = pBar,
                          value = TotalCounter/TotalLength)
      }
      if (nrow(OutPutMatrix) == 1L &
          all(is.na(OutPutMatrix))) {
        OutPutMatrix <- OutPutMatrix[-1L, ]
        # DisplacementMatrix[1, ] <- rep(0L, ncol(DisplacementMatrix))
        OverLapMatrix[1, ] <- rep(0L, ncol(OverLapMatrix))
      }
      ResultMatrix[m1, m2] <- list(OutPutMatrix)
      ResultMatrix[m2, m1] <- list(OverLapMatrix)
    } # end of columns loop
  } # end of rows loop
  if (Verbose) {
    TotalTimeStop <- Sys.time()
    
    cat("\n")
    print(TotalTimeStop - TotalTimeStart)
  }
  dimnames(ResultMatrix) <- dimnames(SyntenyObject)
  class(ResultMatrix) <- "LinkedPairs"
  attr(ResultMatrix, "GeneCalls") <- FeatureRepresentations
  return(ResultMatrix)
}










