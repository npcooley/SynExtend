# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu

PairSummaries <- function(SyntenyLinks,
                          GeneCalls,
                          DBPATH,
                          PIDs = FALSE,
                          IgnoreDefaultStringSet = FALSE,
                          Verbose = FALSE,
                          Model = "Generic",
                          DefaultTranslationTable = "11",
                          ...) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ###### -- Overhead checking -------------------------------------------------
  
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
  
  ###### -- argument passing --------------------------------------------------
  # pass arguments through the ellipsis to either:
  # AlignTranslation,
  # AlignSeqs,
  # or DistanceMatrix
  # if a user supplies an argument that is problematic, such as verbose
  # ignore it
  
  # step one, parse to different functions
  Args <- list(...)
  ArgNames <- names(Args)
  ForbiddenArguments <- c("verbose",
                          "includeTerminalGaps",
                          "type",
                          "myXStringSet",
                          "readingFrame",
                          "sense",
                          "direction",
                          "geneticCode")
  # return(list(Args,
  #             ArgNames))
  
  if (any(ForbiddenArguments %in% ArgNames)) {
    w <- !(ArgNames %in% ForbiddenArguments)
    Args <- Args[w]
    ArgNames <- ArgNames[w]
  }
  
  PossibleATArgs <- unique(c(names(formals(AlignTranslation)),
                             names(formals(AlignSeqs)),
                             names(formals(AlignProfiles))))
  m <- match(x = PossibleATArgs,
             table = ArgNames)
  m <- m[!is.na(m)]
  if (length(m) > 0L) {
    ATArgs <- Args[m]
    rm(PossibleATArgs)
  } else {
    rm(PossibleATArgs)
  }
  PossibleASArgs <- unique(c(names(formals(AlignSeqs)),
                             names(formals(AlignProfiles))))
  m <- match(x = PossibleASArgs,
             table = ArgNames)
  m <- m[!is.na(m)]
  if (length(m) > 0L) {
    # myXStringSet = c(QuerySeqs[m3],
    #                  SubjectSeqs[m3]),
    # verbose = FALSE
    ASArgs <- Args[m]
    rm(PossibleASArgs)
  } else {
    rm(PossibleASArgs)
  }
  PossibleDMArgs <- names(formals(DistanceMatrix))
  m <- match(x = PossibleDMArgs,
             table = ArgNames)
  m <- m[!is.na(m)]
  if (length(m) > 0L) {
    DMArgs <- Args[m]
    rm(PossibleDMArgs)
  } else {
    rm(PossibleDMArgs)
  }
  
  # lifted almost whole cloth from AlignSeqs ...
  # args <- list(...)
  # n <- names(args)
  # m <- character(length(n))
  # for (i in seq_along(n)) {
  #   m[i] <- match.arg(n[i],
  #                     names(c(formals(AlignTranslation),
  #                             formals(AlignSeqs),
  #                             formals(DistanceMatrix))))
  # }
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
  
  ###### -- Deal with Different GeneCall types --------------------------------
  
  L <- length(GeneCalls)
  FeatureRepresentations <- vector(mode = "list",
                                   length = L)
  names(FeatureRepresentations) <- names(GeneCalls)
  
  GCallClasses <- sapply(GeneCalls,
                         function(x) class(x),
                         USE.NAMES = FALSE,
                         simplify = TRUE)
  if (any(GCallClasses == "GRanges")) {
    LimitIndex <- TRUE
    warning("GRanges objects currently only support single contig inputs.")
  }
  
  for (m1 in seq_along(GeneCalls)) {
    if (is(GeneCalls[[m1]],
           "GRanges")) {
      
      if (length(levels(GeneCalls[[m1]]@seqnames)) > 1L) {
        
        warning(paste("GRange object ",
                      m1,
                      " contains more than 1 index and has been truncated.",
                      sep = ""))
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
      
      FeatureRepresentations[[m1]] <- data.frame("Index" = IndexConversion,
                                                 "Strand" = StrandConversion,
                                                 "Start" = StartConversion,
                                                 "Stop" = StopConversion,
                                                 "Lengths" = LengthsConversion,
                                                 "Translation_Table" = rep(NA_character_,
                                                                           length(o)),
                                                 stringsAsFactors = FALSE)
      
      # no synteny object requested, cannot perform this test ...
      # if (any(StopConversion > SyntenyObject[[m1, m1]][1])) {
      #   FeatureRepresentations[[m1]] <- FeatureRepresentations[[m1]][-which(StopConversion > SyntenyObject[[m1, m1]][1]), ]
      # }
      
      rm(list = c("IndexConversion",
                  "StrandConversion",
                  "StartConversion",
                  "StopConversion",
                  "LengthsConversion"))
      
    } else if (is(GeneCalls[[m1]],
                  "DFrame")) {
      
      FeatureRepresentations[[m1]] <- GeneCalls[[m1]]
    } else if (is(GeneCalls[[m1]],
                  "Genes")) {
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
                     "Gene" = ifelse(test = ans[, "Gene"] > 0L,
                                     yes = TRUE,
                                     no = FALSE),
                     "Translation_Table" = rep(NA_character_,
                                               nrow(ans)),
                     "Coding" = ifelse(test = ans[, "Gene"] > 0L,
                                       yes = TRUE,
                                       no = FALSE))
      D <- D[as.vector(ans[, "Gene"]) != 0, ]
      rownames(D) <- NULL
      FeatureRepresentations[[m1]] <- D
      
      rm(list = c("D",
                  "ans",
                  "R"))
    }
  }
  
  GeneCalls <- FeatureRepresentations
  
  rm(list = c("FeatureRepresentations",
              "L"))
  
  ###### -- End Gene call -----------------------------------------------------
  
  ###### -- if a model is specified, load it ----------------------------------
  
  if (!is.null(Model) &
      Model %in% c("Generic")) {
    data("Generic",
         envir = environment(),
         package = "SynExtend")
  }
  
  
  ###### -- Summary stuff -----------------------------------------------------
  
  Size <- dim(SyntenyLinks)[1]
  Total <- (Size^2 - Size) / 2
  Count <- 1L
  PH <- vector(mode = "list",
               length = Total)
  
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
  # QStartDisp == 1
  # QStopDisp == 2
  # SStartDisp == 3
  # SStopDisp == 4
  
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      
      PMatrix <- cbind(SyntenyLinks[[m1, m2]][, 1L],
                       SyntenyLinks[[m1, m2]][, 2L])
      IMatrix <- cbind(SyntenyLinks[[m1, m2]][, 4L],
                       SyntenyLinks[[m1, m2]][, 5L])
      QGeneLength <- GeneCalls[[m1]][PMatrix[, 1L], "Stop"] - GeneCalls[[m1]][PMatrix[, 1L], "Start"] + 1L
      SGeneLength <- GeneCalls[[m2]][PMatrix[, 2L], "Stop"] - GeneCalls[[m2]][PMatrix[, 2L], "Start"] + 1L
      ExactOverLap <- SyntenyLinks[[m1, m2]][, 3L]
      # MinGap <- abs(QGeneLength - SGeneLength)
      TotalKmers <- SyntenyLinks[[m1, m2]][, 11L]
      MaxKmer <- SyntenyLinks[[m1, m2]][, 10L]
      ExteriorMissQuery <- SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L]
      ExteriorMissSubject <- SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L]
      InteriorMissQuery <- QGeneLength - (ExactOverLap + ExteriorMissQuery)
      InteriorMissSubject <- SGeneLength - (ExactOverLap + ExteriorMissSubject)
      QGeneStrand <- GeneCalls[[m1]][PMatrix[, 1L], "Strand"]
      QGeneCoding <- GeneCalls[[m1]][PMatrix[, 1L], "Coding"]
      QGeneTransl <- GeneCalls[[m1]][PMatrix[, 1L], "Translation_Table"]
      SGeneStrand <- GeneCalls[[m2]][PMatrix[, 2L], "Strand"]
      SGeneCoding <- GeneCalls[[m2]][PMatrix[, 2L], "Coding"]
      SGeneTransl <- GeneCalls[[m2]][PMatrix[, 2L], "Translation_Table"]
      
      # collect PIDs if user requests
      # as of the writing of this function extractAt does not recycle x,
      # if in the future it can recycle x, the rep calls interior to extractAt
      # can be removed
      if (PIDs) {
        
        PresQI <- unique(IMatrix[, 1L])
        PresSI <- unique(IMatrix[, 2L])
        QuerySeqs <- vector(mode = "list",
                            length = length(PresQI))
        SubjectSeqs <- vector(mode = "list",
                              length = length(PresSI))
        PresentQRanges <- GeneCalls[[m1]]$Range[PMatrix[, 1L]]
        PresentSRanges <- GeneCalls[[m2]]$Range[PMatrix[, 2L]]
        
        for (m3 in seq_along(PresQI)) {
          CurrentRanges <- PresentQRanges[IMatrix[, 1L] == PresQI[m3]]
          QuerySeqs[[m3]] <- extractAt(x = rep(Genomes[[m1]][PresQI[m3]],
                                               length(CurrentRanges)),
                                       at = CurrentRanges)
          QuerySeqs[[m3]] <- DNAStringSet(sapply(QuerySeqs[[m3]],
                                                 function(x) unlist(x)))
          names(QuerySeqs[[m3]]) <- paste(names(GeneCalls)[m1],
                                          PresQI[m3],
                                          PMatrix[, 1L][IMatrix[, 1L] == PresQI[m3]],
                                          sep = "_")
        }
        QuerySeqs <- do.call(c,
                             QuerySeqs)
        
        for (m3 in seq_along(PresSI)) {
          CurrentRanges <- PresentSRanges[IMatrix[, 2L] == PresSI[m3]]
          SubjectSeqs[[m3]] <- extractAt(x = rep(Genomes[[m2]][PresSI[m3]],
                                                 length(CurrentRanges)),
                                         at = CurrentRanges)
          SubjectSeqs[[m3]] <- DNAStringSet(sapply(SubjectSeqs[[m3]],
                                                   function(x) unlist(x)))
          names(SubjectSeqs[[m3]]) <- paste(names(GeneCalls)[m2],
                                            PresSI[m3],
                                            PMatrix[, 2L][IMatrix[, 2L] == PresSI[m3]],
                                            sep = "_")
        }
        SubjectSeqs <- do.call(c,
                               SubjectSeqs)
        # return SubjectSeqs to the original order
        o <- match(x = paste(names(GeneCalls)[m2],
                                 IMatrix[, 2L],
                                 PMatrix[, 2L],
                                 sep = "_"),
                   table = names(SubjectSeqs))
        SubjectSeqs <- SubjectSeqs[o]
        if (length(QuerySeqs) != length(SubjectSeqs)) {
          stop("Extracted Seqs have differing lengths.")
        }
        
        # reverse complement seqs where necessary
        QuerySeqs[QGeneStrand == 1L] <- reverseComplement(QuerySeqs[QGeneStrand == 1L])
        SubjectSeqs[SGeneStrand == 1L] <- reverseComplement(SubjectSeqs[SGeneStrand == 1L])
        
        # return(list(QuerySeqs, SubjectSeqs, QGeneCoding, SGeneCoding))
        Pident <- vector(mode = "numeric",
                         length = length(QuerySeqs))
        Atype <- vector(mode = "character",
                        length = length(QuerySeqs))
        
        # return(list(QuerySeqs,
        #             SubjectSeqs,
        #             QGeneCoding,
        #             SGeneCoding,
        #             QGeneLength,
        #             SGeneLength,
        #             PMatrix,
        #             IMatrix,
        #             PresQI,
        #             PresSI,
        #             o,
        #             paste(names(GeneCalls)[m2],
        #                   IMatrix[, 2L],
        #                   PMatrix[, 2L],
        #                   sep = "_")))
        if (IgnoreDefaultStringSet) {
          # perform all alignments in nucleotide space
          for (m3 in seq_along(SubjectSeqs)) {
            Atype[m3] <- "NT"
            if ("ASArgs" %in% ls()) {
              CurrentASArgs <- c(list("myXStringSet" = c(QuerySeqs[m3],
                                                         SubjectSeqs[m3]),
                                      "verbose" = FALSE),
                                 ASArgs)
              
            } else {
              CurrentASArgs <- list("myXStringSet" = c(QuerySeqs[m3],
                                                       SubjectSeqs[m3]),
                                    "verbose" = FALSE)
            }
            ph01 <- do.call(what = AlignSeqs,
                            args = CurrentASArgs)
            if ("DMArgs" %in% ls()) {
              CurrentDMArgs <- c(list("myXStringSet" = ph01,
                                      "includeTerminalGaps" = TRUE,
                                      "verbose" = FALSE,
                                      "type" = "matrix"),
                                 DMArgs)
            } else {
              CurrentDMArgs <- list("myXStringSet" = ph01,
                                    "includeTerminalGaps" = TRUE,
                                    "verbose" = FALSE,
                                    "type" = "matrix")
            }
            Pident[m3] <- 1 - do.call(what = DistanceMatrix,
                                      args = CurrentDMArgs)[1, 2]
            # Pident[m3] <- 1 - DistanceMatrix(myXStringSet = AlignSeqs(myXStringSet = c(QuerySeqs[m3],
            #                                                                            SubjectSeqs[m3]),
            #                                                           verbose = FALSE,
            #                                                           ...),
            #                                  includeTerminalGaps = TRUE,
            #                                  verbose = FALSE,
            #                                  type = "matrix",
            #                                  ...)[1, 2]
          }
        } else {
          # perform amino acid alignments where possible
          for (m3 in seq_along(SubjectSeqs)) {
            if (QGeneCoding[m3] &
                SGeneCoding[m3] &
                QGeneLength[m3] %% 3 == 0 &
                SGeneLength[m3] %% 3 == 0) {
              Atype[m3] <- "AA"
              # Pident[m3] <- 1 - DistanceMatrix(myXStringSet = AlignSeqs(myXStringSet = translate(c(QuerySeqs[m3],
              #                                                                                      SubjectSeqs[m3]),
              #                                                                                    if.fuzzy.codon = "solve"),
              #                                                           verbose = FALSE),
              #                                  includeTerminalGaps = TRUE,
              #                                  verbose = FALSE,
              #                                  type = "matrix")[1, 2]
              if ("ATArgs" %in% ls()) {
                gC1 <- if (is.na(QGeneTransl[m3])) {
                  getGeneticCode(id_or_name2 = DefaultTranslationTable,
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                } else {
                  getGeneticCode(id_or_name2 = QGeneTransl[m3],
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                }
                gC2 <- if (is.na(SGeneTransl[m3])) {
                  getGeneticCode(id_or_name2 = DefaultTranslationTable,
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                } else {
                  getGeneticCode(id_or_name2 = SGeneTransl[m3],
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                }
                CurrentATArgs <- c(list("myXStringSet" = c(QuerySeqs[m3],
                                                           SubjectSeqs[m3]),
                                        "readingFrame" = 1,
                                        "sense" = "+",
                                        "direction" = "5' to 3'",
                                        "type" = "DNAStringSet",
                                        "geneticCode" = list(gC1,
                                                             gC2),
                                        "verbose" = FALSE),
                                   ATArgs)
                
              } else {
                gC1 <- if (is.na(QGeneTransl[m3])) {
                  getGeneticCode(id_or_name2 = DefaultTranslationTable,
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                } else {
                  getGeneticCode(id_or_name2 = QGeneTransl[m3],
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                }
                gC2 <- if (is.na(SGeneTransl[m3])) {
                  getGeneticCode(id_or_name2 = DefaultTranslationTable,
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                } else {
                  getGeneticCode(id_or_name2 = SGeneTransl[m3],
                                 full.search = FALSE,
                                 as.data.frame = FALSE)
                }
                CurrentATArgs <- list("myXStringSet" = c(QuerySeqs[m3],
                                                         SubjectSeqs[m3]),
                                      "readingFrame" = 1,
                                      "sense" = "+",
                                      "direction" = "5' to 3'",
                                      "type" = "DNAStringSet",
                                      "geneticCode" = list(gC1,
                                                           gC2),
                                      "verbose" = FALSE)
              }
              # return(CurrentATArgs)
              ph01 <- do.call(what = AlignTranslation,
                              args = CurrentATArgs)
              if ("DMArgs" %in% ls()) {
                CurrentDMArgs <- c(list("myXStringSet" = ph01,
                                        "includeTerminalGaps" = TRUE,
                                        "verbose" = FALSE,
                                        "type" = "matrix"),
                                   DMArgs)
              } else {
                CurrentDMArgs <- list("myXStringSet" = ph01,
                                      "includeTerminalGaps" = TRUE,
                                      "verbose" = FALSE,
                                      "type" = "matrix")
              }
              Pident[m3] <- 1 - do.call(what = DistanceMatrix,
                                        args = CurrentDMArgs)[1, 2]
              # Pident[m3] <- 1 - DistanceMatrix(myXStringSet = AlignTranslation(myXStringSet = c(QuerySeqs[m3],
              #                                                                                   SubjectSeqs[m3]),
              #                                                                  readingFrame = 1,
              #                                                                  sense = "+",
              #                                                                  direction = "5' to 3'",
              #                                                                  type = "DNAStringSet",
              #                                                                  verbose = FALSE,
              #                                                                  ...),
              #                                  includeTerminalGaps = TRUE,
              #                                  verbose = FALSE,
              #                                  type = "matrix",
              #                                  ...)[1, 2]
            } else {
              Atype[m3] <- "NT"
              if ("ASArgs" %in% ls()) {
                CurrentASArgs <- c(list("myXStringSet" = c(QuerySeqs[m3],
                                                           SubjectSeqs[m3]),
                                        "verbose" = FALSE),
                                   ASArgs)
                
              } else {
                CurrentASArgs <- list("myXStringSet" = c(QuerySeqs[m3],
                                                         SubjectSeqs[m3]),
                                      "verbose" = FALSE)
              }
              ph01 <- do.call(what = AlignSeqs,
                              args = CurrentASArgs)
              if ("DMArgs" %in% ls()) {
                CurrentDMArgs <- c(list("myXStringSet" = ph01,
                                        "includeTerminalGaps" = TRUE,
                                        "verbose" = FALSE,
                                        "type" = "matrix"),
                                   DMArgs)
              } else {
                CurrentDMArgs <- list("myXStringSet" = ph01,
                                      "includeTerminalGaps" = TRUE,
                                      "verbose" = FALSE,
                                      "type" = "matrix")
              }
              Pident[m3] <- 1 - do.call(what = DistanceMatrix,
                                        args = CurrentDMArgs)[1, 2]
              # Pident[m3] <- 1 - DistanceMatrix(myXStringSet = AlignSeqs(myXStringSet = c(QuerySeqs[m3],
              #                                                                            SubjectSeqs[m3]),
              #                                                           verbose = FALSE,
              #                                                           ...),
              #                                  includeTerminalGaps = TRUE,
              #                                  verbose = FALSE,
              #                                  type = "matrix",
              #                                  ...)[1, 2]
            }
          }
        }
        
        PH[[Count]] <- data.frame("p1" = names(QuerySeqs),
                                  "p2" = names(SubjectSeqs),
                                  "ExactMatch" = ExactOverLap,
                                  # "MinGap" = MinGap,
                                  "TotalKmers" = TotalKmers,
                                  "MaxKmer" = MaxKmer,
                                  "InteriorQueryMiss" = InteriorMissQuery,
                                  "InteriorSubjectMiss" = InteriorMissSubject,
                                  "ExteriorQueryMiss" = ExteriorMissQuery,
                                  "ExteriorSubjectMiss" = ExteriorMissSubject,
                                  "p1FeatureLength" = QGeneLength,
                                  "p2FeatureLength" = SGeneLength,
                                  "PID" = Pident,
                                  "PIDType" = Atype,
                                  stringsAsFactors = FALSE)
      } else {
        PH[[Count]] <- data.frame("p1" = paste(names(GeneCalls)[m1],
                                               IMatrix[, 1L],
                                               PMatrix[, 1L],
                                               sep = "_"),
                                  "p2" = paste(names(GeneCalls)[m2],
                                               IMatrix[, 2L],
                                               PMatrix[, 2L],
                                               sep = "_"),
                                  "ExactMatch" = ExactOverLap,
                                  # "MinGap" = MinGap,
                                  "TotalKmers" = TotalKmers,
                                  "MaxKmer" = MaxKmer,
                                  "InteriorQueryMiss" = InteriorMissQuery,
                                  "InteriorSubjectMiss" = InteriorMissSubject,
                                  "ExteriorQueryMiss" = ExteriorMissQuery,
                                  "ExteriorSubjectMiss" = ExteriorMissSubject,
                                  "p1FeatureLength" = QGeneLength,
                                  "p2FeatureLength" = SGeneLength,
                                  "PIDType" = ifelse(test = GeneCalls[[m1]][PMatrix[, 1L], "Coding"] &
                                                       GeneCalls[[m2]][PMatrix[, 2L], "Coding"],
                                                     yes = "AA",
                                                     no = "NT"),
                                  stringsAsFactors = FALSE)
      }
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = Count / Total)
      }
      Count <- Count + 1L
    }
  }
  
  DF <- do.call(rbind,
                PH)
  rownames(DF) <- NULL
  
  if (is.null(Model)) {
    # do nothing
  } else if (Model == "Generic") {
    # model constructed from randomly selected prokaryotic pair alignments
    PPids <- predict(object = Generic,
                     DF[])
    DF <- cbind(DF,
                "PPids" = PPids)
  } else if (!is.character(Model) &
             !is.null(Model)) {
    # user supplied model
    PPids <- predict(object = Model,
                     DF)
    DF <- cbind(DF,
                "PPids" = PPids)
  } else if (is.character(Model) &
             !(Model %in% c("Generic"))) {
    cat("\n Selected model is not available.")
  }
  
  
  
  if (Verbose) {
    TimeEnd <- Sys.time()
    cat("\n")
    print(TimeEnd - TimeStart)
  }
  
  return(DF)
}
  
  
  
  
  
  
  
  
  
  