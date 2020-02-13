# Authors: Adelle Fernando, Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu

######
# TODO
# 1: deal with other features - mobile genetic elements, CRISPER arrays, repeat regions?
# 2: Very small genes?
######


gffToDataFrame <- function(GFF,
                           AdditionalAttrs = NULL,
                           AdditionalTypes = NULL,
                           RawTableOnly = FALSE,
                           Verbose = FALSE) {
  
  if (Verbose) {
    FunStartTime <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  if (grepl(pattern = "www.|http:|https:|ftp:|ftps:",
            x = GFF)) {
    CONN <- gzcon(url(GFF))
    Z01 <- readLines(CONN)
    close.connection(con = CONN)
  } else {
    Z01 <- readLines(GFF)
  }
  
  T01 <- strsplit(Z01,
                  split = "\t",
                  fixed = TRUE)
  T01 <- T01[which(lengths(T01) == 9L)]
  T01 <- do.call(rbind,
                 T01)
  T02 <- strsplit(T01[, 9L],
                  split = ";",
                  fixed = TRUE)
  T01 <- T01[, -9L]
  
  # create a DF with the correct raw formats
  
  T03 <- as.data.frame(T01,
                       stringsAsFactors = FALSE)
  
  T03[, 4L] <- as.integer(T03[, 4L])
  T03[, 5L] <- as.integer(T03[, 5L])
  T03[, 7L] <- ifelse(test = T03[, 7L] == "+",
                      yes = 0L,
                      no = 1L)
  
  CollectAttr <- c("ID",
                   "Parent",
                   "Name",
                   "gbkey",
                   "gene",
                   "product",
                   "protein_id",
                   "gene_biotype",
                   "Note")
  
  if (!is.null(AdditionalAttrs)) {
    CollectAttr <- c(CollectAttr,
                     AdditionalAttrs)
    CollectAttr <- unique(CollectAttr)
  }
  
  # rectangularize attribute fields
  T04 <- sapply(T02,
                function(y) sapply(CollectAttr,
                                   function(z) if (length(grep(pattern = paste(z,
                                                                               "=",
                                                                               sep = ""),
                                                               x = y)) == 1L) {
                                     gsub(pattern = paste("(",
                                                          z,
                                                          "=)",
                                                          "(.*)",
                                                          sep = ""),
                                          replacement = "\\2",
                                          x = y[grepl(pattern = paste(z,
                                                                      "=",
                                                                      sep = ""),
                                                      x = y)])
                                   } else if (length(grep(pattern = paste(z,
                                                                          "=",
                                                                          sep = ""),
                                                          x = y)) == 0L) {
                                     NA
                                   } else if (length(grep(pattern = paste(z,
                                                                          "=",
                                                                          sep = ""),
                                                          x = y)) > 1L) {
                                     gsub(pattern = paste("(",
                                                          z,
                                                          "=)",
                                                          "(.*)",
                                                          sep = ""),
                                          replacement = "\\2",
                                          x = y[grepl(pattern = paste(z,
                                                                      "=",
                                                                      sep = ""),
                                                      x = y)][1L])
                                     warning("An attribute field is repeated, only the first is collected.")
                                   }))
  
  T04 <- t(T04)
  T04 <- as.data.frame(T04,
                       stringsAsFactors = FALSE)
  
  colnames(T03) <- c("Contig",
                     "Source",
                     "Type",
                     "Start",
                     "Stop",
                     "Score",
                     "Strand",
                     "Phase")
  
  CompleteTable <- cbind(T03,
                         T04,
                         stringsAsFactors = FALSE)
  
  Contigs <- CompleteTable[!duplicated(CompleteTable$Contig), ]
  ContigMaxes <- Contigs$Stop[order(Contigs$Stop,
                                    decreasing = TRUE)]
  Contigs <- Contigs[order(Contigs$Stop,
                           decreasing = TRUE), ]
  Contigs <- Contigs$Contig
  
  Index <- sapply(CompleteTable$Contig,
                  function(x) which(Contigs == x),
                  USE.NAMES = FALSE,
                  simplify = TRUE)
  
  CompleteTable <- cbind(CompleteTable,
                         "Index" = Index)
  
  CompleteTable <- CompleteTable[order(CompleteTable$Index,
                                       CompleteTable$Start), ]
  
  if (RawTableOnly) {
    
    if (Verbose) {
      cat("\n")
      FunEndTime <- Sys.time()
      print(FunEndTime - FunStartTime)
    }
    
    return(CompleteTable)
  }
  
  if (Verbose) {
    setTxtProgressBar(pb = pBar,
                      value = 1 / 3)
  }
  
  MatchTypes <- c("gene",
                  "pseudogene")
  
  if (!is.null(AdditionalTypes)) {
    MatchTypes <- c(MatchTypes,
                    AdditionalTypes)
    MatchTypes <- unique(MatchTypes)
  }
  
  MatchTable <- CompleteTable[CompleteTable$Type %in% MatchTypes, ]
  MatchTable <- MatchTable[order(MatchTable$Index,
                                 MatchTable$Start), ]
  
  # Match by Parent/ID matching
  
  AllChildrenList <- vector(mode = "list",
                            length = nrow(MatchTable))
  AllCDSChildrenList <- vector(mode = "list",
                               length = nrow(MatchTable))
  
  for (m1 in seq_len(nrow(MatchTable))) {
    ExpandLines <- TRUE
    AllParents <- MatchTable$ID[m1]
    while (ExpandLines) {
      PrevIter <- length(AllParents)
      NewParents <- CompleteTable$ID[CompleteTable$Parent %in% AllParents]
      NewParents <- NewParents[!is.na(NewParents)]
      AllParents <- c(AllParents,
                      NewParents)
      AllParents <- unique(AllParents)
      CurrentIter <- length(AllParents)
      
      if (PrevIter == CurrentIter) {
        ExpandLines <- FALSE
        AllChildren <- CompleteTable[CompleteTable$Parent %in% AllParents, ]
        CDSChildren <- AllChildren[AllChildren$Type == "CDS", ]
      }
      # while loop will repeat here until all child lines are found
    }
    AllCDSChildrenList[[m1]] <- CDSChildren
    AllChildrenList[[m1]] <- AllChildren
    
    if (nrow(AllCDSChildrenList[[m1]]) >= 1L) {
      # Correct for phase:
      AllCDSChildrenList[[m1]][, "Start"] <- mapply(function(x, y) ifelse(test = !is.na(as.integer(x)),
                                                                          yes = y + as.integer(x),
                                                                          no = y),
                                                    x = AllCDSChildrenList[[m1]][, "Phase"],
                                                    y = AllCDSChildrenList[[m1]][, "Start"])
      # remove aberant lines # maybe flip, depending on NCBI response
      if (any(AllCDSChildrenList[[m1]]$Start >= AllCDSChildrenList[[m1]]$Stop)) {
        AllCDSChildrenList[[m1]] <- AllCDSChildrenList[[m1]][!(AllCDSChildrenList[[m1]]$Start >= AllCDSChildrenList[[m1]]$Stop), ]
      }
    }
    
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = (nrow(MatchTable) + m1) / (nrow(MatchTable) * 3L))
    }
  }
  
  # Generate associated character and logical columns ...
  CodingSelect <- vector(mode = "logical",
                         length = nrow(MatchTable))
  MatchLine <- vector(mode = "character",
                      length = nrow(MatchTable))
  ProductLine <- vector(mode = "character",
                        length = nrow(MatchTable))
  NoteLine <- vector(mode = "character",
                     length = nrow(MatchTable))
  ParseNote <- vector(mode = "character",
                      length = nrow(MatchTable))
  
  for (m1 in seq_len(nrow(MatchTable))) {
    if (nrow(AllChildrenList[[m1]]) == 0L &
        nrow(AllCDSChildrenList[[m1]]) == 0L) {
      # Case 1
      # Generate CodingSelect, MatchLine, ProductLine and NoteLine
      CodingSelect[m1] <- FALSE
      MatchLine[m1] <- paste(MatchTable[m1, "Start"],
                             MatchTable[m1, "Stop"],
                             sep = "X")
      if (is.na(MatchTable[m1, "product"])) {
        ProductLine[m1] <- ""
      } else {
        ProductLine[m1] <- MatchTable[m1, "product"]
      }
      if (is.na(MatchTable[m1, "Note"])) {
        NoteLine[m1] <- ""
      } else {
        NoteLine[m1] <- MatchTable[m1, "Note"]
      }
    } else if (nrow(AllCDSChildrenList[[m1]]) == 0L &
               nrow(AllChildrenList[[m1]]) >= 1L) {
      # Case 2
      CodingSelect[m1] <- FALSE
      CurrentTable <- rbind(MatchTable[m1, ],
                            AllChildrenList[[m1]])
      # Parse Product Line
      if (all(is.na(CurrentTable$product))) {
        ProductLine[m1] <- ""
      } else if (length(unique(CurrentTable$product[!is.na(CurrentTable$product)])) == 1L) {
        ProductLine[m1] <- unique(CurrentTable$product[!is.na(CurrentTable$product)])
      } else if (length(unique(CurrentTable$product[!is.na(CurrentTable$product)])) > 1L) {
        ProductLine[m1] <- unique(CurrentTable$product[!is.na(CurrentTable$product)])[1L]
        warning(paste("Multiple Products exist for Row ",
                      m1,
                      ", taking only the first."))
      }
      # Parse Note Line
      if (all(is.na(CurrentTable$Note))) {
        NoteLine[m1] <- ""
      } else if (length(unique(CurrentTable$Note[!is.na(CurrentTable$Note)])) == 1L) {
        NoteLine[m1] <- unique(CurrentTable$Note[!is.na(CurrentTable$Note)])
      } else if (length(unique(CurrentTable$Note[!is.na(CurrentTable$Note)])) > 1L) {
        NoteLine[m1] <- unique(CurrentTable$Note[!is.na(CurrentTable$Note)])[1L]
        warning(paste("Multiple Notes exist for Row ",
                      m1,
                      ", taking only the first."))
      }
      # split into 2A & 2B & 2C
      if (length(unique(CurrentTable$Start)) == 1L &
          length(unique(CurrentTable$Stop)) == 1L) {
        # 2A no disagreements in bounds!
        MatchLine[m1] <- paste(MatchTable[m1, "Start"],
                               MatchTable[m1, "Stop"],
                               sep = "X")
      } else if (xor(length(unique(CurrentTable$Start)) > 1L,
                     length(unique(CurrentTable$Stop)) > 1L)) {
        # 2B if only one bound has a deviation, take the most expansive!
        MatchLine[m1] <- paste(min(CurrentTable$Start),
                               max(CurrentTable$Stop),
                               sep = "X")
      } else if (length(unique(CurrentTable$Start)) > 1L &
                 length(unique(CurrentTable$Stop)) > 1L) {
        # 2C this is an anticipated case with no clear a priori guidance
        MatchLine[m1] <- paste(min(CurrentTable$Start),
                               max(CurrentTable$Stop),
                               sep = "X")
        ParseNote[m1] <- "A discrepancy exists in feature bounds here."
      }
    } else if (nrow(AllCDSChildrenList[[m1]]) >= 1L) {
      # Case 3
      CodingSelect[m1] <- TRUE
      CurrentTable <- rbind(MatchTable[m1, ],
                            AllCDSChildrenList[[m1]])
      # Parse Product Line
      if (all(is.na(CurrentTable$product))) {
        ProductLine[m1] <- ""
      } else if (length(unique(CurrentTable$product[!is.na(CurrentTable$product)])) == 1L) {
        ProductLine[m1] <- unique(CurrentTable$product[!is.na(CurrentTable$product)])
      } else if (length(unique(CurrentTable$product[!is.na(CurrentTable$product)])) > 1L) {
        ProductLine[m1] <- unique(CurrentTable$product[!is.na(CurrentTable$product)])[1L]
        warning(paste("Multiple Products exist for Row ",
                      m1,
                      ", taking only the first."))
      }
      # Parse Note Line
      if (all(is.na(CurrentTable$Note))) {
        NoteLine[m1] <- ""
      } else if (length(unique(CurrentTable$Note[!is.na(CurrentTable$Note)])) == 1L) {
        NoteLine[m1] <- unique(CurrentTable$Note[!is.na(CurrentTable$Note)])
      } else if (length(unique(CurrentTable$Note[!is.na(CurrentTable$Note)])) > 1L) {
        NoteLine[m1] <- unique(CurrentTable$Note[!is.na(CurrentTable$Note)])[1L]
        warning(paste("Multiple Notes exist for Row ",
                      m1,
                      ", taking only the first."))
      }
      
      # split into 3A, 3B, and 3C
      if (length(unique(AllCDSChildrenList[[m1]]$Start)) == 1L &
          length(unique(AllCDSChildrenList[[m1]]$Stop)) == 1L) {
        # 3A no disagreements in bounds!
        # Take the original gene line?!
        # MatchLine[m1] <- paste(MatchTable[m1, "Start"],
        #                        MatchTable[m1, "Stop"],
        #                        sep = "X")
        # vs take the CDS lines
        MatchLine[m1] <- paste(unique(AllCDSChildrenList[[m1]][, "Start"]),
                               unique(AllCDSChildrenList[[m1]][, "Stop"]),
                               sep = "X")
      } else if (length(unique(AllCDSChildrenList[[m1]]$Start)) > 1L &
                 length(unique(AllCDSChildrenList[[m1]]$Stop)) > 1L) {
        # 3B multiple CDSs
        MatchLine[m1] <- paste(apply(X = AllCDSChildrenList[[m1]][, c("Start", "Stop")],
                                     MARGIN = 1,
                                     FUN = function(x) paste(x[1],
                                                             x[2],
                                                             sep = "X")),
                               collapse = "Y")
      } else if (xor(length(unique(AllCDSChildrenList[[m1]]$Start)) > 1L,
                     length(unique(AllCDSChildrenList[[m1]]$Stop)) > 1L)) {
        # only one CDS bound is has multiple positions, take the most expansive CDS bound
        MatchLine[m1] <- paste(min(AllChildrenList[[m1]][, "Start"]),
                               max(AllChildrenList[[m1]][, "Stop"]),
                               sep = "X")
        ParseNote[m1] <- "A discrepancy exists in CDS Bounds here."
      }
    } # end line reconcilation cases
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = ((nrow(MatchTable) * 2) + m1) / (nrow(MatchTable) * 3L))
    }
  }
  
  MatchTable <- cbind(MatchTable,
                      "Match" = MatchLine,
                      "Coding" = CodingSelect,
                      "Product" = ProductLine,
                      "AnnotationNote" = NoteLine,
                      "ParseNotes" = ParseNote,
                      stringsAsFactors = FALSE)
  MatchTable <- MatchTable[, c("Index",
                               "Strand",
                               "Start",
                               "Stop",
                               "Type",
                               "Match",
                               "Product",
                               "AnnotationNote",
                               "gene_biotype",
                               "Coding",
                               "Contig",
                               "ParseNotes")]
  rownames(MatchTable) <- NULL
  
  # rewrite any matches where genes extend over the end of an index
  for (m1 in seq_len(nrow(MatchTable))) {
    if (MatchTable$Stop[m1] > ContigMaxes[MatchTable$Index[m1]]) {
      CurrentBounds <- sapply(strsplit(strsplit(MatchTable$Match[m1],
                                                split = "Y",
                                                fixed = TRUE)[[1]],
                                       split = "X",
                                       fixed = TRUE),
                              function(x) as.integer(x))
      BreakPoint <- which(CurrentBounds[1, ] <= ContigMaxes[MatchTable$Index[m1]] &
                            CurrentBounds[2, ] > ContigMaxes[MatchTable$Index[m1]])
      Section <- CurrentBounds[, BreakPoint]
      NewSection <- matrix(data = c(Section[1],
                                    ContigMaxes[MatchTable$Index[m1]],
                                    1L,
                                    Section[2] - ContigMaxes[MatchTable$Index[m1]]),
                           ncol = 2L)
      if (BreakPoint < ncol(CurrentBounds) &
          BreakPoint != 1L) {
        # if breakpoint is not either the first or last bound set
        AdjustedBounds <- CurrentBounds[,
                                        (BreakPoint + 1L):ncol(CurrentBounds),
                                        drop = FALSE] - ContigMaxes[MatchTable$Index[m1]]
        NewBoundSet <- cbind(CurrentBounds[, 1L:(BreakPoint - 1L)],
                             NewSection,
                             AdjustedBounds)
      } else if (BreakPoint == 1L &
                 ncol(CurrentBounds) > 1L) {
        # if breakboint is the first bound set
        CurrentBounds <- CurrentBounds - ContigMaxes[MatchTable$Index[m1]]
        NewBoundSet <- cbind(NewSection,
                             CurrentBounds[, (BreakPoint + 1L):ncol(CurrentBounds)])
      } else if (BreakPoint == ncol(CurrentBounds) &
                 ncol(CurrentBounds) > 1L) {
        # if breakpoint is the last bound set
        NewBoundSet <- cbind(CurrentBounds[, 1:(BreakPoint - 1L)],
                             NewSection)
      } else {
        # breakpoint is the only bound set
        NewBoundSet <- NewSection
      }
      NewMatchSet <- paste(apply(X = NewBoundSet,
                                 MARGIN = 2L,
                                 FUN = function(x) paste(x,
                                                         collapse = "X")),
                           collapse = "Y")
      MatchTable$Match[m1] <- NewMatchSet
    }
  }
  
  if (Verbose) {
    FunEndTime <- Sys.time()
    cat("\n")
    print(FunEndTime - FunStartTime)
  }
  return(MatchTable)
}






