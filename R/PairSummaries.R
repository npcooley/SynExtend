# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu

PairSummaries <- function(SyntenyLinks,
                          DBPATH,
                          PIDs = FALSE,
                          IgnoreDefaultStringSet = FALSE,
                          Verbose = FALSE,
                          Model = "Generic",
                          DefaultTranslationTable = "11",
                          AcceptContigNames = TRUE,
                          OffSetsAllowed = 2L,
                          ...) {
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ###### -- Overhead checking -------------------------------------------------
  
  # if (length(GeneCalls) != ncol(SyntenyLinks)) {
  #   stop ("LinkedPairs object and gene predictions are not compatible")
  # }
  if (!("DECIPHER" %in% .packages())) {
    stop ("Required package DECIPHER is not loaded")
  }
  if (!is(SyntenyLinks, "LinkedPairs")) {
    stop ("Object is not an LinkedPairs object.")
  }
  if (any(OffSetsAllowed <= 1L)) {
    stop ("Disallowed offsets.")
  }
  # GCallClasses <- sapply(GeneCalls,
  #                        function(x) class(x),
  #                        simplify = TRUE,
  #                        USE.NAMES = FALSE)
  # if (any(GCallClasses == "GRanges")) {
  #   warning("GRanges objects only support Nucleotide Alignments.")
  # }
  
  if (length(OffSetsAllowed) > 0L) {
    AllowGaps <- TRUE
  }
  if (is.null(OffSetsAllowed)) {
    AllowGaps <- FALSE
  }
  
  GeneCalls <- attr(SyntenyLinks, "GeneCalls")
  
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
  
  ###### -- subset gene calls based on the names of the links object ----------
  
  if (length(GeneCalls) != nrow(SyntenyLinks)) {
    GeneCalls <- GeneCalls[match(x = dimnames(SyntenyLinks)[[1]],
                                 table = names(GeneCalls))]
  }
  
  ###### -- extract genomes for stuff -----------------------------------------
  
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
  if (Verbose) {
    cat("\n")
  }
  
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
  if (PIDs) {
    Total <- sum(sapply(SyntenyLinks[upper.tri(SyntenyLinks)],
                        function(x) nrow(x),
                        USE.NAMES = FALSE,
                        simplify = TRUE))
    # read this message after neighbors and k-mer dists ?
    cat("Aligning pairs.\n")
  } else {
    Total <- (Size^2 - Size) / 2
    
    cat("Collecting pairs.\n")
  }
  
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
        
        IMatrix <- cbind(SyntenyLinks[[m1, m2]][, 4L],
                         SyntenyLinks[[m1, m2]][, 5L])
        
        # map neighbors here - generate LN + RN columns
        # loop through possible index combos
        # then remap PMatrix and IMatrix
        # include gap fills here? or later?
        UIM <- unique(IMatrix) # unique index matrix
        IndexKey <- match(x = data.frame(t(IMatrix)),
                          table = data.frame(t(UIM)))
        UIK <- unique(IndexKey)
        
        # return(list(UIM,
        #             UIK,
        #             IndexKey,
        #             IMatrix))
        
        LKey <- RKey <- NeighborMat <- vector(mode = "list",
                                              length = nrow(UIM))
        for (m3 in seq_len(nrow(UIM))) {
          # don't need to bother with subsetting index matrix here
          CIM <- IMatrix[IndexKey == UIK[m3], , drop = FALSE] # current index matrix
          CPM <- PMatrix[IndexKey == UIK[m3], , drop = FALSE] # current index matrix
          
          if (nrow(CPM) > 1L) {
            p1 <- CPM[, 1]
            p2 <- CPM[, 2]
            i1 <- CIM[, 1]
            i2 <- CIM[, 2]
            
            p1e <- rle(p1)
            up1 <- p1e$values
            p2e <- rle(p2)
            up2 <- p2e$values
            
            # if one id is present while being present many times just rep zero
            # else solve
            if (length(up1) > 1L) {
              # successive rightward p1
              rdp1 <- c(up1[2L:length(up1)] - up1[1L:(length(up1) - 1L)], 0L)
              # successive leftward p1
              ldp1 <- abs(c(0L, up1[1L:(length(up1) - 1L)] - up1[2L:length(up1)]))
              # expand unique diffs back to pair list size
              rdp1 <- rep(rdp1,
                          p1e$lengths)
              ldp1 <- rep(ldp1,
                          p1e$lengths)
            } else {
              rdp1 <- rep(0L, p1e$lengths)
              ldp1 <- rep(0L, p1e$lengths)
            }
            
            if (length(up2) > 1L) {
              # successive rightward p2
              rdp2 <- c(up2[2L:length(up2)] - up2[1L:(length(up2) - 1L)], 0L)
              # successive leftward p2
              ldp2 <- c(0L, up2[1L:(length(up2) - 1L)] - up2[2L:length(up2)]) * -1L
              # expand unique diffs back to pair list size
              rdp2 <- rep(rdp2,
                          p2e$lengths)
              ldp2 <- rep(ldp2,
                          p2e$lengths)
            } else {
              rdp2 <- rep(0L, p2e$lengths)
              ldp2 <- rep(0L, p2e$lengths)
            }
            
            # left and right are absolute in p1
            # but left and right are relative to diagonal / anti-diagonal in p2
            # when p2 left or right is negative it is pointing along the anti-diagonal
            # you can have neighbors on either diagonal, but to have a gap, those neighbors must
            # be on the same diagonal
            
            NeighborMat[[m3]] <- cbind("p1" = p1,
                                       "p2" = p2,
                                       "i1" = i1,
                                       "i2" = i2,
                                       "p1rd" = rdp1,
                                       "p1ld" = ldp1,
                                       "p2rd" = rdp2,
                                       "p2ld" = ldp2)
            RKey[[m3]] <- as.integer(NeighborMat[[m3]][, 5L] == 1L & abs(NeighborMat[[m3]][, 7L]) == 1L)
            LKey[[m3]] <- as.integer(NeighborMat[[m3]][, 6L] == 1L & abs(NeighborMat[[m3]][, 8L]) == 1L)
            
          } else if (nrow(CPM) == 1L) {
            # only a single gene appears in this index combo
            LKey[[m3]] <- RKey[[m3]] <- 0L # assign keys as false
            NeighborMat[[m3]] <- cbind("p1" = CPM[, 1L],
                                       "p2" = CPM[, 2L],
                                       "i1" = CIM[, 1L],
                                       "i2" = CIM[, 2L],
                                       "p1rd" = 0L,
                                       "p1ld" = 0L,
                                       "p2rd" = 0L,
                                       "p2ld" = 0L)
          }
        } # end first first pass of neighbor assignment
        
        RKey <- unlist(RKey)
        LKey <- unlist(LKey)
        NeighborMat <- do.call(rbind,
                               NeighborMat)
        
        o1 <- order(NeighborMat[, 3L],
                    NeighborMat[, 1L],
                    NeighborMat[, 2L])
        NeighborMat <- NeighborMat[o1, ]
        
        # return(list(NeighborMat,
        #             RKey,
        #             LKey,
        #             PMatrix,
        #             IMatrix,
        #             o1,
        #             UIK,
        #             IndexKey))
        
        # Create a matrix of gap filled positions
        # if PID calc is specified, calculate them
        # if not, don't bother?
        
        
        # if gaps are allowed fill them
        # the choice is being made here to leave gap fills as neighborless
        # this is intentional
        
        if (AllowGaps &
            nrow(PMatrix) > 1L) {
          
          GapFill <- vector(mode = "list",
                            length = length(OffSetsAllowed))
          
          # diff gives absolute difference to next neighbor
          # already exists in the neighbors matrix
          p1R <- NeighborMat[, 5L]
          p2R <- NeighborMat[, 7L]
          
          # for each gap size allowed:
          for (g1 in seq_along(OffSetsAllowed)) {
            # find where gap size equals allowed size
            w1 <- abs(p1R) == OffSetsAllowed[g1]
            w2 <- abs(p2R) == OffSetsAllowed[g1]
            # check that the gap does not span indices
            
            w3 <- w1 & w2
            if (sum(w3) > 0L) {
              # gap to next is the same distance for both partner columns
              # at least once
              # build vectors of those gene positions
              # and the opposing positions spanning the gap
              # don't overwrite things you need until you're finished with them
              GapFill[[g1]] <- vector(mode = "list",
                                      length = sum(w3))
              w4 <- which(w3) + 1L
              i1l <- NeighborMat[w3, 3L]
              i1r <- NeighborMat[w4, 3L]
              i2l <- NeighborMat[w3, 4L]
              i2r <- NeighborMat[w4, 4L]
              p1l <- NeighborMat[w3, 1L]
              p1r <- NeighborMat[w4, 1L]
              p2l <- NeighborMat[w3, 2L]
              p2r <- NeighborMat[w4, 2L]
              # create new pair partner lines
              # if gap does not span indices
              for (g2 in seq_along(p1l)) {
                if (i1l[g2] == i1r[g2] &
                    i2l[g2] == i2r[g2]) {
                  # gap does not span indices fill in based on gap size
                  # this is already checked earlier and might not be necessary?
                  gp1 <- seq(from = p1l[g2],
                             to = p1r[g2],
                             by = if (p1l[g2] < p1r[g2]) {
                               1L
                             } else {
                               -1L
                             })
                  gp2 <- seq(from = p2l[g2],
                             to = p2r[g2],
                             by = if (p2l[g2] < p2r[g2]) {
                               1L
                             } else {
                               -1L
                             })
                  # if (length(gp1) != length(gp2)) {
                  #   return(list(gp1,
                  #               gp2,
                  #               p1l[g2],
                  #               p1r[g2],
                  #               p2l[g2],
                  #               p2r[g2],
                  #               i1l[g2],
                  #               i1r[g2],
                  #               i2l[g2],
                  #               i2r[g2],
                  #               NeighborMat))
                  # }
                  GapFill[[g1]][[g2]] <- cbind("g1" = rep(m1, OffSetsAllowed[g1] - 1L),
                                               "i1" = rep(i1l[g2], OffSetsAllowed[g1] - 1L),
                                               "p1" = gp1[-c(1, length(gp1))],
                                               "g2" = rep(m2, OffSetsAllowed[g1] - 1L),
                                               "i2" = rep(i2l[g2], OffSetsAllowed[g1] - 1L),
                                               "p2" = gp2[-c(1, length(gp2))])
                  
                }
              }
              # return(GapFill)
              GapFill[[g1]] <- do.call(rbind,
                                       GapFill[[g1]])
              Ins1Str <- GeneCalls[[m1]][GapFill[[g1]][, 3L], "Strand"]
              Ins1Coding <- GeneCalls[[m1]][GapFill[[g1]][, 3L], "Coding"]
              Ins1Transl <- GeneCalls[[m1]][GapFill[[g1]][, 3L], "Translation_Table"]
              Ins2Str <- GeneCalls[[m2]][GapFill[[g1]][, 6L], "Strand"]
              Ins2Coding <- GeneCalls[[m2]][GapFill[[g1]][, 6L], "Coding"]
              Ins2Transl <- GeneCalls[[m2]][GapFill[[g1]][, 6L], "Translation_Table"]
              Ins1GLength <- GeneCalls[[m1]][GapFill[[g1]][, 3L], "Stop"] - GeneCalls[[m1]][GapFill[[g1]][, 3L], "Start"] + 1L
              Ins2GLength <- GeneCalls[[m2]][GapFill[[g1]][, 6L], "Stop"] - GeneCalls[[m2]][GapFill[[g1]][, 6L], "Start"] + 1L
              Ins1IMiss <- Ins1EMiss <- Ins1GLength
              Ins2IMiss <- Ins2EMiss <- Ins2GLength
              InsOv <- InsMax <- InsTot <- rep(0L, nrow(GapFill[[g1]]))
            } else {
              # in this case do ... something?
              # leave list position as null
              # GapFill[[g1]]
            }
          }
          # return(GapFill)
          GapFill <- do.call(rbind,
                             GapFill)
          
          if (!is.null(GapFill)) {
            # gaps were spanned
            # combine and order vectors
            
            pmat1 <- rbind(IMatrix,
                           GapFill[, c(2,5)])
            pmat2 <- rbind(PMatrix,
                           GapFill[, c(3,6)])
            ExactOverLap <- c(SyntenyLinks[[m1, m2]][, 3L],
                              rep(0L, nrow(GapFill)))
            TotalKmers <- c(SyntenyLinks[[m1, m2]][, 11L],
                            rep(0L, nrow(GapFill)))
            MaxKmer <- c(SyntenyLinks[[m1, m2]][, 10L],
                         rep(0L, nrow(GapFill)))
            ExteriorMissQuery <- c(SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L],
                                   GeneCalls[[m1]][GapFill[, 3L], "Stop"] - GeneCalls[[m1]][GapFill[, 3L], "Start"] + 1L)
            ExteriorMissSubject <- c(SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L],
                                     GeneCalls[[m2]][GapFill[, 6L], "Stop"] - GeneCalls[[m2]][GapFill[, 6L], "Start"] + 1L)
            
            o1 <- order(pmat1[, 1L],
                        pmat2[, 1L],
                        pmat2[, 2L])
            ExactOverLap <- ExactOverLap[o1]
            TotalKmers <- TotalKmers[o1]
            MaxKmer <- MaxKmer[o1]
            ExteriorMissQuery <- ExteriorMissQuery[o1]
            ExteriorMissSubject <- ExteriorMissSubject[o1]
            
            IMatrix <- pmat1[o1, ]
            PMatrix <- pmat2[o1, ]
            # regenerate neighbor key 
            RKey <- c(RKey, rep(0L, nrow(GapFill)))
            LKey <- c(LKey, rep(0L, nrow(GapFill)))
            RKey <- RKey[o1]
            LKey <- LKey[o1]
            # return(list(IMatrix,
            #             PMatrix))
            # index matching
            QGeneLength <- GeneCalls[[m1]][PMatrix[, 1L], "Stop"] - GeneCalls[[m1]][PMatrix[, 1L], "Start"] + 1L
            SGeneLength <- GeneCalls[[m2]][PMatrix[, 2L], "Stop"] - GeneCalls[[m2]][PMatrix[, 2L], "Start"] + 1L
            # ExactOverLap <- SyntenyLinks[[m1, m2]][, 3L]
            # MinGap <- abs(QGeneLength - SGeneLength)
            # TotalKmers <- SyntenyLinks[[m1, m2]][, 11L]
            # MaxKmer <- SyntenyLinks[[m1, m2]][, 10L]
            # ExteriorMissQuery <- SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L]
            # ExteriorMissSubject <- SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L]
            InteriorMissQuery <- QGeneLength - (ExactOverLap + ExteriorMissQuery)
            InteriorMissSubject <- SGeneLength - (ExactOverLap + ExteriorMissSubject)
            QGeneStrand <- GeneCalls[[m1]][PMatrix[, 1L], "Strand"]
            QGeneCoding <- GeneCalls[[m1]][PMatrix[, 1L], "Coding"]
            QGeneTransl <- GeneCalls[[m1]][PMatrix[, 1L], "Translation_Table"]
            SGeneStrand <- GeneCalls[[m2]][PMatrix[, 2L], "Strand"]
            SGeneCoding <- GeneCalls[[m2]][PMatrix[, 2L], "Coding"]
            SGeneTransl <- GeneCalls[[m2]][PMatrix[, 2L], "Translation_Table"]
            PairLeft <- LKey
            PairRight <- RKey
            
          } else {
            # do nothing, no gaps discovered
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
            PairLeft <- LKey
            PairRight <- RKey
          }
        } else if (!AllowGaps |
                   nrow(PMatrix) == 1L) {
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
          PairLeft <- LKey
          PairRight <- RKey
        } # End gap checking
        
        # collect PIDs if user requests
        # as of the writing of this function extractAt does not recycle x,
        # if in the future it can recycle x, the rep calls interior to extractAt
        # can be removed
        
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
          CurrentRanges <- PresentSRanges[IMatrix[, 2L] == PresSI[m3]]
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
        
        NucDist <- vector(mode = "numeric",
                          length = length(QuerySeqs))
        nuc1 <- oligonucleotideFrequency(x = QuerySeqs,
                                         width = 4L,
                                         as.prob = TRUE)
        nuc2 <- oligonucleotideFrequency(x = SubjectSeqs,
                                         width = 4L,
                                         as.prob = TRUE)
        for (m3 in seq_along(NucDist)) {
          NucDist[m3] <- sqrt(sum((nuc1[m3, ] - nuc2[m3, ])^2)) / ((sum(nuc1[m3, ]) + sum(nuc2[m3, ])) / 2)
        }
        
        
        if (PIDs) {
          
          Pident <- vector(mode = "numeric",
                           length = length(QuerySeqs))
          Atype <- vector(mode = "character",
                          length = length(QuerySeqs))
          
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
              
              if (Verbose) {
                setTxtProgressBar(pb = pBar,
                                  value = Count / Total)
                Count <- Count + 1L
              }
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
              
              if (Verbose) {
                setTxtProgressBar(pb = pBar,
                                  value = Count / Total)
                Count <- Count + 1L
              }
            } # end m3 loop
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
                                    "Adjacent" = RKey + LKey,
                                    "FourmerDist" = NucDist,
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
                                    "Adjacent" = RKey + LKey,
                                    "FourmerDist" = NucDist,
                                    "PIDType" = ifelse(test = GeneCalls[[m1]][PMatrix[, 1L], "Coding"] &
                                                         GeneCalls[[m2]][PMatrix[, 2L], "Coding"],
                                                       yes = "AA",
                                                       no = "NT"),
                                    stringsAsFactors = FALSE)
        }
      } else {
        # no links in table, leave list position as NULL
      }
      if (Verbose &
          !PIDs) {
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
  attr(DF, "GeneCalls") <- attr(SyntenyLinks, "GeneCalls")
  class(DF) <- c("data.frame", "PairSummaries")
  return(DF)
}
  
  
  
  
  
  
  
  
  
  