# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu

PairSummaries <- function(SyntenyLinks,
                          GeneCalls,
                          DBPATH,
                          PIDs = TRUE,
                          IgnoreDefaultStringSet = FALSE,
                          Verbose = TRUE,
                          GapPenalty = TRUE,
                          TerminalPenalty = TRUE,
                          Model = "Global",
                          Correction = "none") {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ###### -- Overhead checking -------------------------------------------------
  
  if (any(lengths(SyntenyLinks[lower.tri(SyntenyLinks)]) == 0L)) {
    stop ("LinkedPairs object must not be in 'Sparse format'")
  }
  if (length(GeneCalls) != ncol(SyntenyLinks)) {
    stop ("LinkedPairs object and gene predictions are not compatible")
  }
  if (!("DECIPHER" %in% .packages())) {
    stop ("Required package DECIPHER is not loaded")
  }
  if (!is(SyntenyLinks, "LinkedPairs")) {
    stop ("Object is not an LinkedPairs object.")
  }
  GCallClasses <- sapply(GeneCalls,
                         function(x) class(x),
                         simplify = TRUE,
                         USE.NAMES = FALSE)
  if (any(GCallClasses == "GRanges")) {
    warning("GRanges objects only support Nucleotide Alignments.")
  }
  
  if (PIDs) {
    Genomes <- vector("list",
                      length = length(GeneCalls))
    for (i in seq_along(Genomes)) {
      Genomes[[i]] <- SearchDB(dbFile = DBPATH,
                               identifier = names(GeneCalls[i]),
                               nameBy = "identifier",
                               verbose = FALSE)
    }
    Genomes <- DNAStringSetList(Genomes)
  }
  
  ###### -- Deal with GRanges objects -----------------------------------------
  
  for (m1 in seq_along(GeneCalls)) {
    if (GCallClasses[m1] == "GRanges") {
      if (length(levels(GeneCalls[[m1]]@seqnames)) > 1L) {
        
        warning("GRange object ",
                m1,
                " contains more than 1 index and has been truncated.")
        IndexPlaceHolder <- as.character(GeneCalls[[m1]]@seqnames)[1L]
        GeneCalls[[m1]] <- GeneCalls[[m1]][as.character(GeneCalls[[m1]]@seqnames == IndexPlaceHolder), ]
        
      }
      TypePlaceHolder <- as.character(GeneCalls[[m1]]$type)
      GeneCalls[[m1]] <- GeneCalls[[m1]][TypePlaceHolder %in% c("gene",
                                                                "pseudogene"), ]
      StrandConversion <- ifelse(test = as.character(GeneCalls[[m1]]@strand == "+"),
                                 yes = 0L,
                                 no = 1L)
      StartConversion <- GeneCalls[[m1]]@ranges@start
      StopConversion <- GeneCalls[[m1]]@ranges@start + GeneCalls[[m1]]@ranges@width - 1L
      IndexConversion <- rep(1L,
                             length(StartConversion))
      LengthsConversion <- StopConversion - StartConversion + 1L
      
      o <- order(StartConversion)
      StrandConversion <- StrandConversion[o]
      IndexConversion <- IndexConversion[o]
      StartConversion <- StartConversion[o]
      StopConversion <- StopConversion[o]
      LengthsConversion <- LengthsConversion[o]
      
      NewRanges <- vector(mode = "list",
                          length = length(o))
      for (m2 in seq_along(NewRanges)) {
        NewRanges[[m2]] <- IRanges(start = StartConversion[m2],
                                   end = StopConversion[m2])
      }
      NewRanges <- IRangesList(NewRanges)
      CodingSelect <- rep(FALSE,
                          length(o))
      GRangeIndices <- rep(1L,
                           length(o))
      GeneCalls[[m1]] <- DataFrame("Index" = GRangeIndices,
                                   "Strand" = StrandConversion,
                                   "Start" = StartConversion,
                                   "Stop" = StopConversion,
                                   "Range" = NewRanges,
                                   "Coding" = CodingSelect)
      
      
      
      if (any(StopConversion > length(Genomes[[m1]][[1]]))) {
        GeneCalls[[m1]] <- GeneCalls[[m1]][-which(StopConversion > length(Genomes[[m1]][[1]])), ]
      }
    } else if (GCallClasses[m1] == "Genes") {
      # convert Erik's gene calls to a temporary DataFrame
      # the column "Gene" assigns whether or not that particular line is the gene
      # that the caller actually picked, calls must be subset to where Gene == 1
      ans <- GeneCalls[[m1]]
      R <- mapply(function(x, y) IRanges(start = x,
                                         end = y),
                  x = ans[, "Begin"],
                  y = ans[, "End"],
                  SIMPLIFY = FALSE)
      D <- DataFrame("Index" = as.integer(ans[, "Index"]),
                     "Strand" = as.integer(ans[, "Strand"]),
                     "Start" = as.integer(ans[, "Begin"]),
                     "Stop" = as.integer(ans[, "End"]),
                     "Type" = rep("gene",
                                  nrow(ans)),
                     "Range" = IRangesList(R),
                     "Coding" = rep(TRUE,
                                    nrow(ans)))
      D <- D[ans[, "Gene"] == 1L, ]
      rownames(D) <- NULL
      GeneCalls[[m1]] <- D
      rm(c("D", "ans", "R"))
      
    }
  }
  
  Size <- dim(SyntenyLinks)[1]
  Total <- (Size^2 - Size) / 2
  
  TotalCoverage <- MinCoverage <- MaxCoverage <- vector("list",
                                                        length = Total)
  PairMatrix <- IndexMatrix <- vector("list",
                                      length = Total)
  QueryGeneLength <- SubjectGeneLength <- CombinedGeneLength <- vector("list",
                                                                       length = Total)
  NormGeneDiff <- AbsStartDelta <- AbsStopDelta <- NormDeltaStart <- NormDeltaStop <- vector("list",
                                                                                             length = Total)
  Scores <- QueryCharacter <- SubjectCharacter <- LabelNames <- vector("list",
                                                                       length = Total)
  IDType <- ExactMatch <- vector(mode = "list",
                                 length = Total)
  ExactLeftQ <- ExactRightQ <- ExactLeftS <- ExactRightS <- vector(mode = "list",
                                                                   length = Total)
  Count <- 1L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      o <- order(GeneCalls[[m1]][, "Index"],
                 GeneCalls[[m1]][, "Start"],
                 decreasing = FALSE)
      GeneCalls[[m1]] <- GeneCalls[[m1]][o, ]
      o <- order(GeneCalls[[m2]][, "Index"],
                 GeneCalls[[m2]][, "Start"],
                 decreasing = FALSE)
      GeneCalls[[m2]] <- GeneCalls[[m2]][o, ]
      ######
      # Find the distance from the closest hit to the start of the gene
      # in nucleotide space
      # for both the subject and the query
      # take the delta of these differences, and convert to a percentage of the combined gene length
      # do the same for stops
      ######
      AbsStartDelta[[Count]] <- abs(SyntenyLinks[m2, m1][[1]][, "QueryStartDisplacement"] - SyntenyLinks[m2, m1][[1]][, "SubjectStartDisplacement"])
      AbsStopDelta[[Count]] <- abs(SyntenyLinks[m2, m1][[1]][, "QueryStopDisplacement"] - SyntenyLinks[m2, m1][[1]][, "SubjectStopDisplacement"])
      PairMatrix[[Count]] <- cbind(SyntenyLinks[m1, m2][[1]][, "QueryGene"],
                                   SyntenyLinks[m1, m2][[1]][, "SubjectGene"])
      IndexMatrix[[Count]] <- cbind(SyntenyLinks[m1, m2][[1]][, "QueryIndex"],
                                    SyntenyLinks[m1, m2][[1]][, "SubjectIndex"])
      QueryGeneLength[[Count]] <- GeneCalls[[m1]][PairMatrix[[Count]][, 1L], "Stop"] - GeneCalls[[m1]][PairMatrix[[Count]][, 1L], "Start"] + 1L
      SubjectGeneLength[[Count]] <- GeneCalls[[m2]][PairMatrix[[Count]][, 2L], "Stop"] - GeneCalls[[m2]][PairMatrix[[Count]][, 2L], "Start"] + 1L
      CombinedGeneLength[[Count]] <- QueryGeneLength[[Count]] + SubjectGeneLength[[Count]]
      NormGeneDiff[[Count]] <- abs(QueryGeneLength[[Count]] - SubjectGeneLength[[Count]]) / CombinedGeneLength[[Count]]
      NormDeltaStart[[Count]] <- AbsStartDelta[[Count]] / CombinedGeneLength[[Count]]
      NormDeltaStop[[Count]] <- AbsStopDelta[[Count]] / CombinedGeneLength[[Count]]
      TotalCoverage[[Count]] <- (SyntenyLinks[m1, m2][[1]][, "ExactOverlap"] * 2L) / CombinedGeneLength[[Count]]
      ExactMatch[[Count]] <- SyntenyLinks[m1, m2][[1]][, "ExactOverlap"]
      ExactLeftQ[[Count]] <- SyntenyLinks[m1, m2][[1]][, "QLeftPos"]
      ExactRightQ[[Count]] <- SyntenyLinks[m1, m2][[1]][, "QRightPos"]
      ExactLeftS[[Count]] <- SyntenyLinks[m1, m2][[1]][, "SLeftPos"]
      ExactRightS[[Count]] <- SyntenyLinks[m1, m2][[1]][, "SRightPos"]
      MinCoverage[[Count]] <- SyntenyLinks[m1, m2][[1]][, "ExactOverlap"] / apply(cbind(QueryGeneLength[[Count]],
                                                                                           SubjectGeneLength[[Count]]),
                                                                                     MARGIN = 1L,
                                                                                     FUN = function(x) min(x))
      MaxCoverage[[Count]] <- SyntenyLinks[m1, m2][[1]][, "ExactOverlap"] / apply(cbind(QueryGeneLength[[Count]],
                                                                                           SubjectGeneLength[[Count]]),
                                                                                     MARGIN = 1L,
                                                                                     FUN = function(x) max(x))
      QueryCharacter[[Count]] <- vector("character",
                                        length = nrow(SyntenyLinks[[m1, m2]]))
      SubjectCharacter[[Count]] <- vector("character",
                                          length = nrow(SyntenyLinks[[m1, m2]]))
      LabelNames[[Count]] <- vector("character",
                                    length = nrow(SyntenyLinks[[m1, m2]]))
      if (PIDs) {
        Scores[[Count]] <- vector("list",
                                  length = nrow(SyntenyLinks[[m1, m2]]))
        IDType[[Count]] <- vector(mode = "character",
                                  length = nrow(SyntenyLinks[[m1, m2]]))
      }
      
      for (i in seq_len(nrow(SyntenyLinks[m1, m2][[1]]))) {
        QueryCharacter[[Count]][i] <- paste(m1,
                                            IndexMatrix[[Count]][i, 1L],
                                            PairMatrix[[Count]][i, 1L],
                                            sep = "_")
        SubjectCharacter[[Count]][i] <- paste(m2,
                                              IndexMatrix[[Count]][i, 2L],
                                              PairMatrix[[Count]][i, 2L],
                                              sep = "_")
        LabelNames[[Count]][i] <- paste(QueryCharacter[[Count]][i],
                                        SubjectCharacter[[Count]][i],
                                        sep = " ")
        if (PIDs) {
          
          QuerySeq <- DNAStringSet(unlist(extractAt(x = Genomes[[m1]][[SyntenyLinks[[m1, m2]][i, "QueryIndex"]]],
                                                    at = GeneCalls[[m1]]$Range[[SyntenyLinks[[m1, m2]][i, "QueryGene"]]])))
          SubjectSeq <- DNAStringSet(unlist(extractAt(x = Genomes[[m2]][[SyntenyLinks[[m1, m2]][i, "SubjectIndex"]]],
                                                      at = GeneCalls[[m2]]$Range[[SyntenyLinks[[m1, m2]][i, "SubjectGene"]]])))
          Scores[[Count]][[i]] <- c(QuerySeq,
                                    SubjectSeq)
          if (GeneCalls[[m1]][SyntenyLinks[m1, m2][[1]][i, "QueryGene"], "Strand"] == 1L) {
            Scores[[Count]][[i]][1] <- reverseComplement(Scores[[Count]][[i]][1])
          }
          if (GeneCalls[[m2]][SyntenyLinks[m1, m2][[1]][i, "SubjectGene"], "Strand"] == 1L) {
            Scores[[Count]][[i]][2] <- reverseComplement(Scores[[Count]][[i]][2])
          }
          
          if (!(GeneCalls[[m1]][SyntenyLinks[[m1, m2]][i, "QueryGene"], "Coding"]) |
              !(GeneCalls[[m2]][SyntenyLinks[[m1, m2]][i, "SubjectGene"], "Coding"]) |
              IgnoreDefaultStringSet) {
            # align as nucleotides
            PlaceHolder01 <- AlignSeqs(myXStringSet = Scores[[Count]][[i]],
                                       verbose =	FALSE)
            PlaceHolder02 <- "DNA"
          } else {
            # align as AAs
            PlaceHolder01 <- AlignTranslation(myXStringSet = Scores[[Count]][[i]],
                                              sense = "+",
                                              direction = "5' to 3'",
                                              type = "DNAStringSet",
                                              readingFrame = 1L,
                                              verbose = FALSE)
            PlaceHolder02 <- "AA"
          } # end select AA vs DNA Alignment
          IDType[[Count]][i] <- PlaceHolder02
          Scores[[Count]][[i]] <- 1 - DistanceMatrix(myXStringSet = PlaceHolder01, 
                                                     penalizeGapLetterMatches = GapPenalty,
                                                     includeTerminalGaps = TerminalPenalty,
                                                     correction = Correction,
                                                     type = "matrix",
                                                     verbose = FALSE)[1, 2]
        } # end similarity scores conditional
      } # end similarity scores loop
      
      ######
      # Go to next matrix position,
      # Assign next list position
      ######
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = Count / Total)
      }
      Count <- Count + 1L
    } # end of columns loop
  } # end of rows loop
  
  DF <- data.frame("TotalCoverage" = unlist(TotalCoverage),
                   "MaxCoverage" = unlist(MaxCoverage),
                   "MinCoverage" = unlist(MinCoverage),
                   "NormDeltaStart" = unlist(NormDeltaStart),
                   "NormDeltaStop" = unlist(NormDeltaStop),
                   "NormGeneDiff" = unlist(NormGeneDiff),
                   "ExactMatch" = unlist(ExactMatch),
                   stringsAsFactors = FALSE)
  if (PIDs) {
    DF <- cbind(DF,
                "PID" = unlist(Scores),
                "PIDType" = unlist(IDType))
  } else {
    # do nothing
  }
  rownames(DF) <- unlist(LabelNames)
  
  ###### -- load in prebuilt, or user specified model -------------------------
  # If link statistics imply that PID will likely give a PID that falls from
  # a distribution of random PIDs, remove that link
  
  if (is.null(Model)) {
    # do nothing
    KeepSet <- rep(TRUE,
                   ncol(DF))
    DF <- cbind(DF,
                "ModelSelect" = KeepSet)
  } else {
    if (Model == "Global") {
      data("GlobalSelect",
           envir = environment(),
           package = "SynExtend")
      KeepSet <- predict(object = GlobalSelect,
                         DF,
                         type = "response")
    } else if (Model == "Local") {
      data("LocalSelect",
           envir = environment(),
           package = "SynExtend")
      KeepSet <- predict(object = LocalSelect,
                         DF,
                         type = "response")
    } else if (Model == "Exact") {
      data("ExactSelect",
           envir = environment(),
           package = "SynExtend")
      KeepSet <- predict(object = ExactSelect,
                         DF,
                         type = "response")
    } else if (!is.character(Model) &
               !is.null(Model)) {
      # User specified model ...
      KeepSet <- predict(object = Model,
                         DF,
                         type = "response")
    }
    DF <- cbind(DF,
                "ModelSelect" = KeepSet >= 0.5)
  }
  
  if (Verbose) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(DF)
}

