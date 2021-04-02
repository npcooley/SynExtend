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
                          AcceptContigNames = TRUE,
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
    if (Verbose) {
      cat("\nSelecting Genomes.\n")
    }
    Genomes <- vector("list",
                      length = length(GeneCalls))
    for (i in seq_along(Genomes)) {
      Genomes[[i]] <- SearchDB(dbFile = DBPATH,
                               identifier = names(GeneCalls[i]),
                               nameBy = "identifier",
                               verbose = FALSE)
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = i / length(Genomes))
      }
    }
    Genomes <- DNAStringSetList(Genomes)
  }
  
  ###### -- Indices present by name -------------------------------------------
  
  ContigsPresent <- sapply(diag(SyntenyLinks),
                           function(x) names(x),
                           simplify = FALSE)
  
  ###### -- Deal with Different GeneCall types --------------------------------
  
  L <- length(GeneCalls)
  FeatureRepresentations <- vector(mode = "list",
                                   length = L)
  ContigNames <- vector(mode = "list",
                        length = L)
  GCallClasses <- sapply(GeneCalls,
                         function(x) class(x),
                         USE.NAMES = FALSE,
                         simplify = TRUE)
  
  if (Verbose) {
    cat("\nReconciling Genecalls.\n")
  }
  
  for (m1 in seq_along(GeneCalls)) {
    if (is(GeneCalls[[m1]],
           "GRanges")) {
      ContigNames[[m1]] <- seq(length(GeneCalls[[m1]]@seqnames@lengths))
      names(ContigNames[[m1]]) <- unique(as.character(GeneCalls[[m1]]@seqnames))
      
      if (any(is.na(match(x = names(ContigNames[[m1]]),
                          table = ContigsPresent[[m1]])))) {
        stop (paste0("Contig names imply incorrectly matched objects at diag position ",
                     names(GeneCalls)[m1]))
      }
      
      TypePlaceHolder <- as.character(GeneCalls[[m1]]$type)
      GeneCalls[[m1]] <- GeneCalls[[m1]][TypePlaceHolder %in% c("gene",
                                                                "pseudogene"), ]
      StrandConversion <- ifelse(test = as.character(GeneCalls[[m1]]@strand == "+"),
                                 yes = 0L,
                                 no = 1L)
      StartConversion <- GeneCalls[[m1]]@ranges@start
      StopConversion <- GeneCalls[[m1]]@ranges@start + GeneCalls[[m1]]@ranges@width - 1L
      IndexConversion <- rep(unname(ContigNames[[m1]]),
                             GeneCalls[[m1]]@seqnames@lengths)
      LengthsConversion <- StopConversion - StartConversion + 1L
      C.Index <- IndexConversion
      ph1 <- unname(ContigNames[[m1]])
      ph2 <- match(x = names(ContigNames[[m1]]),
                   table = ContigsPresent[[m1]])
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
      StrandMax <- rep(unname(SyntenyLinks[[m1, m1]]),
                       table(IndexConversion))
      
      R <- mapply(function(x, y) IRanges(start = x,
                                         end = y),
                  x = StartConversion,
                  y = StopConversion,
                  SIMPLIFY = FALSE)
      
      # FeatureRepresentations[[m1]] <- data.frame("Index" = IndexConversion,
      #                                            "Strand" = StrandConversion,
      #                                            "Start" = StartConversion,
      #                                            "Stop" = StopConversion,
      #                                            "Lengths" = LengthsConversion,
      #                                            "Translation_Table" = rep(NA_character_,
      #                                                                      length(o)),
      #                                            "Coding" = rep(FALSE,
      #                                                           length(o)),
      #                                            stringsAsFactors = FALSE)
      FeatureRepresentations[[m1]] <- DataFrame("Index" = IndexConversion,
                                                "Strand" = StrandConversion,
                                                "Start" = StartConversion,
                                                "Stop" = StopConversion,
                                                "Lengths" = LengthsConversion,
                                                "Translation_Table" = rep(NA_character_,
                                                                          length(o)),
                                                "Coding" = rep(FALSE,
                                                               length(o)),
                                                "Range" = IRangesList(R))
      FeatureRepresentations[[m1]] <- FeatureRepresentations[[m1]][StopConversion <= StrandMax, ]
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
      
      # FeatureRepresentations[[m1]] <- GeneCalls[[m1]]
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
                     table = ContigsPresent[[m1]])
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
      FeatureRepresentations[[m1]] <- GeneCalls[[m1]][o, ]
      # print(FeatureRepresentations[[m1]])
      
    } else if (is(GeneCalls[[m1]],
                  "Genes")) {
      # convert Erik's gene calls to a temporary DataFrame
      # the column "Gene" assigns whether or not that particular line is the gene
      # that the caller actually picked, calls must be subset to where Gene != 0
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
                     table = ContigsPresent[[m1]])
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
      D <- D[as.vector(ans[, "Gene"]) != 0, ]
      rownames(D) <- NULL
      FeatureRepresentations[[m1]] <- D
      ContigNames[[m1]] <- unique(D$Index)
      names(ContigNames[[m1]]) <- ContigNames[[m1]]
      
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
    cat("\nGeneCalls reconciled.\n")
  }
  
  LinkContigNames <- diag(SyntenyLinks)
  names(FeatureRepresentations) <- names(GeneCalls)
  GeneCalls <- FeatureRepresentations
  
  rm(list = c("FeatureRepresentations",
              "L"))
  
  ###### -- End Gene call -----------------------------------------------------
  
  ###### -- if a model is specified, load it ----------------------------------
  
  if (!is.null(Model) &
      Model %in% c("Generic")) {
    MOD <- get(data(list = "Generic",
                    envir = environment(),
                    package = "SynExtend"))
    
  } else if (!is.null(Model) &
             !(Model %in% c("Generic"))) {
    if (file.exists(Model)) {
      MOD <- get(load(file = Model,
                      verbose = FALSE))
      if (!is(object = MOD,
              class2 = "glm")) {
        stop ("\nUser specified model is not a glm?")
      }
    } else {
      stop ("\nUser specified file does not appear to exist.\n")
    }
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
      
      if (nrow(SyntenyLinks[[m1, m2]]) > 0L) {
        # links table is populated, do whatever
        PMatrix <- cbind(SyntenyLinks[[m1, m2]][, 1L],
                         SyntenyLinks[[m1, m2]][, 2L])
        
        # if (AcceptContigNames) {
        #   ph1 <- SyntenyLinks[[m1, m2]][, 4L]
        #   for (m3 in seq_along(unname(LinkContigNames[[m1]]))) {
        #     ph1 <- replace(x = ph1,
        #                    list = which(ph1 == unname(LinkContigNames[[m1]])[m3]),
        #                    values = unname(ContigNames[[m1]])[m3])
        #   }
        #   ph2 <- SyntenyLinks[[m1, m2]][, 5L]
        #   for (m3 in seq_along(unname(LinkContigNames[[m2]]))) {
        #     ph2 <- replace(x = ph2,
        #                    list = which(ph2 == unname(LinkContigNames[[m2]])[m3]),
        #                    values = unname(ContigNames[[m2]])[m3])
        #   }
        #   IMatrix <- cbind(ph1,
        #                    ph2)
        #   rm(list = c("ph1",
        #               "ph2"))
        # } else {
        #   IMatrix <- cbind(SyntenyLinks[[m1, m2]][, 4L],
        #                    SyntenyLinks[[m1, m2]][, 5L])
        # }
        
        # if (m1 == 4 &
        #     m2 == 8) {
        #   return(list(IMatrix,
        #               cbind("m1" = SyntenyLinks[[m1, m2]][, 1L],
        #                     "m2" = SyntenyLinks[[m1, m2]][, 2L]),
        #               SyntenyLinks[[m1, m2]]))
        # }
        # rownames(IMatrix) <- NULL
        # rownames(PMatrix) <- NULL
        
        IMatrix <- cbind(SyntenyLinks[[m1, m2]][, 4L],
                         SyntenyLinks[[m1, m2]][, 5L])
        
        # index matching
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
          
          # return(list(IMatrix,
          #             PMatrix,
          #             PresQI,
          #             PresSI,
          #             PresentQRanges,
          #             PresentSRanges,
          #             GeneCalls[[m1]],
          #             GeneCalls[[m2]],
          #             Genomes[[m1]],
          #             Genomes[[m2]]))
          
          # print(IMatrix)
          # print(PMatrix)
          
          for (m3 in seq_along(PresQI)) {
            # print(m3)
            # print(c(m3, m2, m1))
            # if (m1 == 1 &
            #     m2 == 4 &
            #     m3 == 2) {
            #   return(list(IMatrix,
            #               PMatrix,
            #               PresQI,
            #               PresentQRanges))
            # }
            CurrentRanges <- PresentQRanges[IMatrix[, 1L] == PresQI[m3]]
            # if (m3 == 3L) {
            #   return(list(PresentQRanges,
            #               PresQI,
            #               IMatrix,
            #               Genomes[[m1]][PresQI[m3]],
            #               CurrentRanges,
            #               GeneCalls[[m1]]))
            # }
            QuerySeqs[[m3]] <- extractAt(x = rep(Genomes[[m1]][PresQI[m3]],
                                                 length(CurrentRanges)),
                                         at = unname(CurrentRanges))
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
            # print(c(m3, m2, m1))
            
            CurrentRanges <- PresentSRanges[IMatrix[, 2L] == PresSI[m3]]
            # return(list(IMatrix,
            #             PMatrix,
            #             PresQI,
            #             PresSI,
            #             PresentQRanges,
            #             PresentSRanges,
            #             Genomes,
            #             CurrentRanges,
            #             GeneCalls[[m1]],
            #             GeneCalls[[m2]],
            #             m3))
            SubjectSeqs[[m3]] <- extractAt(x = rep(Genomes[[m2]][PresSI[m3]],
                                                   length(CurrentRanges)),
                                           at = unname(CurrentRanges))
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
                                    "p1InteriorMiss" = InteriorMissQuery,
                                    "p2InteriorMiss" = InteriorMissSubject,
                                    "p1ExteriorMiss" = ExteriorMissQuery,
                                    "p2ExteriorMiss" = ExteriorMissSubject,
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
                                    "p1InteriorMiss" = InteriorMissQuery,
                                    "p2InteriorMiss" = InteriorMissSubject,
                                    "p1ExteriorMiss" = ExteriorMissQuery,
                                    "p2ExteriorMiss" = ExteriorMissSubject,
                                    "p1FeatureLength" = QGeneLength,
                                    "p2FeatureLength" = SGeneLength,
                                    "PIDType" = ifelse(test = GeneCalls[[m1]][PMatrix[, 1L], "Coding"] &
                                                         GeneCalls[[m2]][PMatrix[, 2L], "Coding"],
                                                       yes = "AA",
                                                       no = "NT"),
                                    stringsAsFactors = FALSE)
        }
      } else {
        # no links in table, leave list position as NULL
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
  
  # return(list(DF,
  #             MOD))
  if (!is.null(DF)) {
    if (is.null(Model)) {
      # do nothing
    } else {
      PPids <- predict(object = MOD,
                       DF,
                       type = "response")
      DF <- cbind(DF,
                  "PredictedPID" = PPids)
    }
  } else {
    DF <- NULL
  }
  
  if (Verbose) {
    TimeEnd <- Sys.time()
    cat("\n")
    print(TimeEnd - TimeStart)
  }
  
  return(DF)
}
  
  
  
  
  
  
  
  
  
  