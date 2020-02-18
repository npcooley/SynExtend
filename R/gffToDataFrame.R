# Authors: Adelle Fernando, Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu
# TODO
# 1: deal with other features - mobile genetic elements, CRISPER arrays, repeat regions?
# 2: Very small genes?


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
    InitLines <- readLines(CONN)
    close.connection(con = CONN)
  } else if (grepl(pattern = ".gz$",
                   x = GFF)) {
    CONN <- gzfile(GFF)
    InitLines <- readLines(CONN)
    close.connection(con = CONN)
  } else if (grepl(pattern = ".gff$",
                   x = GFF)) {
    InitLines <- readLines(GFF)
  } else {
    stop("File does not have a recognizable extension.")
  }
  
  LinesAsList <- strsplit(InitLines,
                          split = "\t",
                          fixed = TRUE)
  LinesAsList <- LinesAsList[which(lengths(LinesAsList) == 9L)]
  LinesAsList <- do.call(rbind,
                         LinesAsList)
  RawAttrField <- strsplit(LinesAsList[, 9L],
                           split = ";",
                           fixed = TRUE)
  LinesAsList <- LinesAsList[, -9L]
  
  # create a DF with the correct raw formats
  
  StandardizedColumns <- as.data.frame(LinesAsList,
                                       stringsAsFactors = FALSE)
  
  StandardizedColumns[, 4L] <- as.integer(StandardizedColumns[, 4L])
  StandardizedColumns[, 5L] <- as.integer(StandardizedColumns[, 5L])
  StandardizedColumns[, 7L] <- ifelse(test = StandardizedColumns[, 7L] == "+",
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
  ParsedAttrs <- sapply(RawAttrField,
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
  
  ParsedAttrs <- t(ParsedAttrs)
  ParsedAttrs <- as.data.frame(ParsedAttrs,
                               stringsAsFactors = FALSE)
  
  colnames(StandardizedColumns) <- c("Contig",
                     "Source",
                     "Type",
                     "Start",
                     "Stop",
                     "Score",
                     "Strand",
                     "Phase")
  
  CompleteTable <- cbind(StandardizedColumns,
                         ParsedAttrs,
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
      # This is not what phase is for ...
      # AllCDSChildrenList[[m1]][, "Start"] <- mapply(function(x, y) ifelse(test = !is.na(as.integer(x)),
      #                                                                     yes = y + as.integer(x),
      #                                                                     no = y),
      #                                               x = AllCDSChildrenList[[m1]][, "Phase"],
      #                                               y = AllCDSChildrenList[[m1]][, "Start"])
      AllCDSChildrenList[[m1]] <- AllCDSChildrenList[[m1]][order(AllCDSChildrenList[[m1]]$Start,
                                                                 decreasing = FALSE), ]
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
  
  # Generate associated character and logical columns and lists
  CodingSelect <- vector(mode = "logical",
                         length = nrow(MatchTable))
  MatchLine <- vector(mode = "list",
                      length = nrow(MatchTable))
  ProductLine <- vector(mode = "character",
                        length = nrow(MatchTable))
  NoteLine <- vector(mode = "character",
                     length = nrow(MatchTable))
  ParseNote <- vector(mode = "list",
                      length = nrow(MatchTable))
  
  for (m1 in seq_len(nrow(MatchTable))) {
    if (nrow(AllChildrenList[[m1]]) == 0L &
        nrow(AllCDSChildrenList[[m1]]) == 0L) {
      # Case 1
      # Generate CodingSelect, MatchLine, ProductLine and NoteLine
      CodingSelect[m1] <- FALSE
      MatchLine[[m1]] <- IRanges(start = MatchTable[m1, "Start"],
                                 end = MatchTable[m1, "Stop"])
      if (is.na(MatchTable[m1, "product"])) {
        ProductLine[m1] <- ""
      } else {
        ProductLine[m1] <- MatchTable[m1, "product"]
      }
      if (is.na(MatchTable[m1, "Note"])) {
        NoteLine[m1] <- ""
      } else {
        NoteLine[[m1]] <- MatchTable[m1, "Note"]
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
        ParseNote[[m1]] <- c(NoteLine[[m1]],
                             paste("Multiple Products exist for Row ",
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
        ParseNote[[m1]] <- c(NoteLine[[m1]],
                             paste("Multiple Notes exist for Row ",
                                   m1,
                                   ", taking only the first."))
      }
      # split into 2A & 2B & 2C
      if (length(unique(CurrentTable$Start)) == 1L &
          length(unique(CurrentTable$Stop)) == 1L) {
        # 2A no disagreements in bounds!
        MatchLine[[m1]] <- IRanges(start = MatchTable[m1, "Start"],
                                   end = MatchTable[m1, "Stop"])
      } else if (xor(length(unique(CurrentTable$Start)) > 1L,
                     length(unique(CurrentTable$Stop)) > 1L)) {
        # 2B if only one bound has a deviation, take the most expansive!
        MatchLine[[m1]] <- IRanges(start = min(CurrentTable$Start),
                                   end = max(CurrentTable$Stop))
      } else if (length(unique(CurrentTable$Start)) > 1L &
                 length(unique(CurrentTable$Stop)) > 1L) {
        # 2C this is an anticipated case with no clear a priori guidance
        MatchLine[[m1]] <- IRanges(start = min(CurrentTable$Start),
                                   end = max(CurrentTable$Stop))
        ParseNote[[m1]] <- c(ParseNote[[m1]],
                             "A discrepancy exists in feature bounds here.")
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
        ParseNote[[m1]] <- c(NoteLine[[m1]],
                             paste("Multiple Products exist for Row ",
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
        ParseNote[[m1]] <- c(NoteLine[[m1]],
                             paste("Multiple Notes exist for Row ",
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
        MatchLine[[m1]] <- IRanges(start = unique(AllCDSChildrenList[[m1]][, "Start"]),
                                   end = unique(AllCDSChildrenList[[m1]][, "Stop"]))
        # case where CDSs are differing by ends?starts?
      } else if (length(unique(AllCDSChildrenList[[m1]]$Start)) > 1L &
                 length(unique(AllCDSChildrenList[[m1]]$Stop)) > 1L) {
        # 3B multiple CDSs
        # MatchLine[m1] <- paste(apply(X = AllCDSChildrenList[[m1]][, c("Start", "Stop")],
        #                              MARGIN = 1,
        #                              FUN = function(x) paste(x[1],
        #                                                      x[2],
        #                                                      sep = "X")),
        #                        collapse = "Y")
        MatchLine[[m1]] <- IRanges(start = AllCDSChildrenList[[m1]]$Start,
                                   end = AllCDSChildrenList[[m1]]$Stop)
      } else if (xor(length(unique(AllCDSChildrenList[[m1]]$Start)) > 1L,
                     length(unique(AllCDSChildrenList[[m1]]$Stop)) > 1L)) {
        # only one CDS bound is has multiple positions, take the most expansive CDS bound
        # MatchLine[m1] <- paste(min(AllChildrenList[[m1]][, "Start"]),
        #                        max(AllChildrenList[[m1]][, "Stop"]),
        #                        sep = "X")
        MatchLine[[m1]] <- IRanges(start = min(AllChildrenList[[m1]][, "Start"]),
                                   end = max(AllChildrenList[[m1]][, "Stop"]))
        ParseNote[[m1]] <- c(ParseNote[[m1]],
                             "A discrepancy exists in CDS Bounds here.")
      }
    } # end line reconcilation cases
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = ((nrow(MatchTable) * 2) + m1) / (nrow(MatchTable) * 3L))
    }
  }
  
  MatchTable <- DataFrame(MatchTable,
                          "Range" = IRangesList(MatchLine),
                          "Coding" = CodingSelect,
                          # "AnnotationNote" = NoteLine,
                          # "ParseNotes" = ParseNote,
                          # stringsAsFactors = FALSE,
                          "Product" = ProductLine)
  
  MatchTable <- MatchTable[, c("Index",
                               "Strand",
                               "Start",
                               "Stop",
                               "Type",
                               "ID",
                               "Range",
                               "Product",
                               # "AnnotationNote",
                               # "gene_biotype",
                               "Coding",
                               "Contig")]
  rownames(MatchTable) <- NULL
  
  # rewrite any ranges where bounds extend over the end of an index
  for (m1 in seq_len(nrow(MatchTable))) {
    if (any((MatchLine[[m1]]@start + MatchLine[[m1]]@width - 1L) > ContigMaxes[MatchTable$Index[m1]])) {
      # at least one bound extends over the end of the index
      B_Starts <- MatchLine[[m1]]@start # Bound starts
      B_Ends <- MatchLine[[m1]]@start + MatchLine[[m1]]@width - 1L # Bound ends
      F_Break <- ContigMaxes[MatchTable$Index[m1]] # FASTA break
      # do this with a while loop
      RemainingBounds <- TRUE
      C_Counts <- 1L
      N_Counts <- 1L
      N_Starts <- rep(NA_integer_,
                      length(B_Starts) * 2L)
      N_Ends <- rep(NA_integer_,
                    length(B_Ends) * 2L)
      while(RemainingBounds) {
        if (F_Break >= B_Ends[C_Counts] &
            C_Counts == length(B_Ends)) {
          # print("break after, no further evaluations")
          # Case 1:
          # Break occurs after this bound set
          # This bound set is the only remaining set
          # add bounds to new vectors
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- B_Ends[C_Counts]
          RemainingBounds <- FALSE
        } else if (F_Break < B_Starts[C_Counts] &
                   C_Counts == length(B_Starts)) {
          # print("break prior, no further evaluations")
          # Case 2:
          # Break occurs before this bound set
          # This bound set is the only remaining set
          # adjust current bound set
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts] - F_Break
          N_Ends[N_Counts] <- B_Ends[C_Counts] - F_Break
          RemainingBounds <- FALSE
        } else if (F_Break >= B_Starts[C_Counts] &
                   F_Break < B_Ends[C_Counts] &
                   C_Counts == length(B_Starts)) {
          # print("break interior, no further evaluations")
          # Case 3:
          # Break occurs inside current bound set
          # this is the only remaining bound set
          # add set with adjusted current end
          # add new set
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- F_Break
          N_Starts[N_Counts + 1L] <- 1L
          N_Ends[N_Counts + 1L] <- B_Ends[C_Counts] - F_Break
          RemainingBounds <- FALSE
        } else if (F_Break >= B_Ends[C_Counts] &
                   C_Counts < length(B_Ends)) {
          # print("break after, further evaluations exist")
          # Case 4:
          # Break occurs after this bound set
          # Bound sets remain to be evaluated
          # add current set whole cloth
          # update counters
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- B_Ends[C_Counts]
          N_Counts <- N_Counts + 1L
          C_Counts <- C_Counts + 1L
        } else if (F_Break < B_Starts[C_Counts] &
                   C_Counts < length(B_Starts)) {
          # print("break prior, further evaluations exist")
          # Case 5:
          # Break occurs before this bound set
          # Bound sets remain to be evaluated
          # Adjust all following bound sets, but not current
          # update counters
          # remaining sets will be visited in future
          B_Starts[(C_Counts + 1L):length(B_Starts)] <- B_Starts[(C_Counts + 1L):length(B_Starts)] - F_Break
          B_Ends[(C_Counts + 1L):length(B_Ends)] <- B_Ends[(C_Counts + 1L):length(B_Ends)] - F_Break
          N_Starts[N_Counts] <- B_Starts[C_Counts] - F_Break
          N_Ends[N_Counts] <- B_Ends[C_Counts] - F_Break
          N_Counts <- N_Counts + 1L
          C_Counts <- C_Counts + 1L
        } else if (F_Break >= B_Starts[C_Counts] &
                   F_Break < B_Ends[C_Counts]) {
          # print("break interior, further evaluations exist")
          # Case 6:
          # Break occurs inside current bound set
          # Bound sets remain to be evaluated
          # Adjust following bound sets
          # add set with adjusted current end
          # add new set
          # update following bounds
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- F_Break
          N_Starts[N_Counts + 1L] <- 1L
          N_Ends[N_Counts + 1L] <- B_Ends[C_Counts] - F_Break # check off by 1!
          B_Starts[(C_Counts + 1L):length(B_Starts)] <- B_Starts[(C_Counts + 1L):length(B_Starts)] - F_Break
          B_Ends[(C_Counts + 1L):length(B_Ends)] <- B_Ends[(C_Counts + 1L):length(B_Ends)] - F_Break
          C_Counts <- C_Counts + 1L
          N_Counts <- N_Counts + 2L
        } else {
          print("an unmet condition exists")
          stop("An unmet condition was discovered while adjusting bounds, please contact maintainer!")
        }
      } # while loop
      N_Starts <- N_Starts[!is.na(N_Starts)]
      N_Ends <- N_Ends[!is.na(N_Ends)]
      
      MatchTable$Range[[m1]] <- IRanges(start = N_Starts,
                                        end = N_Ends)
    } # conditional to check
  } # end loop
  
  if (Verbose) {
    FunEndTime <- Sys.time()
    cat("\n")
    print(FunEndTime - FunStartTime)
  }
  return(MatchTable)
}






