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
                          Storage = 1,
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
  if (Storage < 0) {
    stop("Storage must be at least zero.")
  } else {
    Storage <- Storage * 1e9 # conversion to gigabytes
  }
  Size <- dim(SyntenyLinks)[1]
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
                          "pattern",
                          "subject",
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
  # align with only align profiles
  PossibleAPArgs <- unique(c(names(formals(AlignProfiles))))
  m <- match(x = PossibleAPArgs,
             table = ArgNames)
  m <- m[!is.na(m)]
  if (length(m) > 0L) {
    APArgs <- Args[m]
    rm(PossibleAPArgs)
  } else {
    rm(PossibleAPArgs)
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
  
  # load in genomes and ALL extracted features at the top until storage limit
  # is hit
  if (Verbose) {
    cat("\nPreparing overhead data.\n")
  }
  Features <- vector("list",
                     length = length(GeneCalls))
  L <- length(GeneCalls)
  
  Count <- 1L
  while (object.size(Features) < Storage &
         Count <= L) {
    # print(object.size(Features))
    # print(Count)
    Genome <- SearchDB(dbFile = DBPATH,
                       identifier = names(GeneCalls[Count]),
                       nameBy = "description",
                       verbose = FALSE)
    PresentIndices <- unique(GeneCalls[[Count]]$Index)
    
    # return(list(Genome,
    #             PresentIndices,
    #             GeneCalls[[Count]],
    #             Count,
    #             GeneCalls))
    if (length(PresentIndices) > 1L) {
      # many indices, loop through present indices and extract
      # slam together at the end
      Features[[Count]] <- vector(mode = "list",
                                  length = length(PresentIndices))
      for (m3 in seq_along(PresentIndices)) {
        ph <- GeneCalls[[Count]]$Index == PresentIndices[m3]
        Features[[Count]][[m3]] <- extractAt(x = rep(Genome[PresentIndices[m3]],
                                                     sum(ph)),
                                             at = unname(GeneCalls[[Count]]$Range[ph]))
        Features[[Count]][[m3]] <- DNAStringSet(sapply(Features[[Count]][[m3]],
                                                       function(x) unlist(x)))
        FlipMe <- GeneCalls[[Count]]$Strand[ph] == 1L
        if (any(FlipMe)) {
          Features[[Count]][[m3]][FlipMe] <- reverseComplement(Features[[Count]][[m3]][FlipMe])
        }
      }
      Features[[Count]] <- do.call(c,
                                   Features[[Count]])
      
    } else {
      # only 1 index present in gene calls
      Features[[Count]] <- extractAt(x = rep(Genome[PresentIndices],
                                             nrow(GeneCalls[[Count]])),
                                     at = unname(GeneCalls[[Count]]$Range))
      Features[[Count]] <- DNAStringSet(sapply(Features[[Count]],
                                               function(x) unlist(x)))
      FlipMe <- GeneCalls[[Count]]$Strand == 1L
      if (any(FlipMe)) {
        Features[[Count]][FlipMe] <- reverseComplement(Features[[Count]][FlipMe])
      }
      
    }
    names(Features[[Count]]) <- paste(rep(names(GeneCalls)[Count], length(Features[[Count]])),
                                   GeneCalls[[Count]]$Index,
                                   seq(length(Features[[Count]])),
                                   sep = "_")
    Count <- Count + 1L
    # will extract till storage is exceeded
    # will not cap at storage
  }
  if (Verbose) {
    if (Count < L) {
      cat("Overhead is too large to keep entirely in memory.\nPrimary loop will include database lookups.\n")
    } else {
      cat("Overhead complete.\n")
    }
  }
  
  ###### -- if a model is specified, load it ----------------------------------
  
  # if (!is.null(Model) &
  #     Model %in% c("Generic")) {
  #   MOD <- get(data(list = "Generic",
  #                   envir = environment(),
  #                   package = "SynExtend"))
  #   
  # } else if (!is.null(Model) &
  #            !(Model %in% c("Generic"))) {
  #   if (file.exists(Model)) {
  #     MOD <- get(load(file = Model,
  #                     verbose = FALSE))
  #     if (!is(object = MOD,
  #             class2 = "glm")) {
  #       stop ("\nUser specified model is not a glm?\n")
  #     }
  #   } else {
  #     stop ("\nUser specified file does not appear to exist.\n")
  #   }
  # }
  
  if (!is.null(Model)) {
    if (Model %in% "Generic") {
      MOD <- get(data(list = "Generic",
                      envir = environment(),
                      package = "SynExtend"))
    } else {
      if (file.exists(Model)) {
        MOD <- get(load(file = Model,
                        verbose = FALSE))
        if (!is(object = MOD,
                class2 = "glm")) {
          stop ("\nUser specified model is not a glm?\n")
        }
      } else {
        stop ("\nUser specified file does not appear to exist.\n")
      }
    }
  } else {
    # model is null do nothing
  }
  
  ###### -- Summary stuff -----------------------------------------------------
  
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
  # same minus max and total, but for individual linking kmers
  # not all linking kmers
  
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      
      if (nrow(SyntenyLinks[[m1, m2]]) > 0L) {
        # links table is populated, do whatever
        PMatrix <- cbind(SyntenyLinks[[m1, m2]][, 1L],
                         SyntenyLinks[[m1, m2]][, 2L])
        
        IMatrix <- cbind(SyntenyLinks[[m1, m2]][, 4L],
                         SyntenyLinks[[m1, m2]][, 5L])
        # find the concensus of linking ks
        p1l <- GeneCalls[[m1]]$Stop[SyntenyLinks[[m2, m1]][, 1L]] - GeneCalls[[m1]]$Start[SyntenyLinks[[m2, m1]][, 1L]] + 1L
        p2l <- GeneCalls[[m2]]$Stop[SyntenyLinks[[m2, m1]][, 2L]] - GeneCalls[[m2]]$Start[SyntenyLinks[[m2, m1]][, 2L]] + 1L
        a1 <- SyntenyLinks[[m2, m1]][, 6L]
        a2 <- SyntenyLinks[[m2, m1]][, 7L]
        b1 <- SyntenyLinks[[m2, m1]][, 8L]
        b2 <- SyntenyLinks[[m2, m1]][, 9L]
        c1 <- GeneCalls[[m1]]$Start[SyntenyLinks[[m2, m1]][, 1L]]
        c2 <- GeneCalls[[m1]]$Stop[SyntenyLinks[[m2, m1]][, 1L]]
        d1 <- GeneCalls[[m2]]$Start[SyntenyLinks[[m2, m1]][, 2L]]
        d2 <- GeneCalls[[m2]]$Stop[SyntenyLinks[[m2, m1]][, 2L]]
        s1 <- GeneCalls[[m1]]$Strand[SyntenyLinks[[m2, m1]][, 1L]] == 0L
        s2 <- GeneCalls[[m2]]$Strand[SyntenyLinks[[m2, m1]][, 2L]] == 0L
        
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
        
        diff2 <- vector(mode = "numeric",
                        length = nrow(SyntenyLinks[[m1, m2]]))
        for (m3 in seq_along(diff2)) {
          diff2[m3] <- 1 - mean(diff1[SyntenyLinks[[m2, m1]][, 1L] == SyntenyLinks[[m1, m2]][m3, 1L] &
                                        SyntenyLinks[[m2, m1]][, 2L] == SyntenyLinks[[m1, m2]][m3, 2L]])
        }
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
                  
                  if (length(gp1) == length(gp2)) {
                    GapFill[[g1]][[g2]] <- cbind("g1" = rep(m1, OffSetsAllowed[g1] - 1L),
                                                 "i1" = rep(i1l[g2], OffSetsAllowed[g1] - 1L),
                                                 "p1" = gp1[-c(1, length(gp1))],
                                                 "g2" = rep(m2, OffSetsAllowed[g1] - 1L),
                                                 "i2" = rep(i2l[g2], OffSetsAllowed[g1] - 1L),
                                                 "p2" = gp2[-c(1, length(gp2))])
                  }
                }
              }
              # return(GapFill)
              GapFill[[g1]] <- do.call(rbind,
                                       GapFill[[g1]])
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
            diff2 <- c(diff2,
                       rep(0, nrow(GapFill)))
            ExactOverLap <- c(SyntenyLinks[[m1, m2]][, 3L],
                              rep(0L, nrow(GapFill)))
            TotalKmers <- c(SyntenyLinks[[m1, m2]][, 11L],
                            rep(0L, nrow(GapFill)))
            MaxKmer <- c(SyntenyLinks[[m1, m2]][, 10L],
                         rep(0L, nrow(GapFill)))
            # ExteriorMissQuery <- c(SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L],
            #                        GeneCalls[[m1]][GapFill[, 3L], "Stop"] - GeneCalls[[m1]][GapFill[, 3L], "Start"] + 1L)
            # ExteriorMissSubject <- c(SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L],
            #                          GeneCalls[[m2]][GapFill[, 6L], "Stop"] - GeneCalls[[m2]][GapFill[, 6L], "Start"] + 1L)
            
            o1 <- order(pmat1[, 1L],
                        pmat2[, 1L],
                        pmat2[, 2L])
            ExactOverLap <- ExactOverLap[o1]
            TotalKmers <- TotalKmers[o1]
            MaxKmer <- MaxKmer[o1]
            # ExteriorMissQuery <- ExteriorMissQuery[o1]
            # ExteriorMissSubject <- ExteriorMissSubject[o1]
            
            IMatrix <- pmat1[o1, ]
            PMatrix <- pmat2[o1, ]
            diff2 <- diff2[o1]
            # regenerate neighbor key 
            RKey <- c(RKey, rep(0L, nrow(GapFill)))
            LKey <- c(LKey, rep(0L, nrow(GapFill)))
            RKey <- RKey[o1]
            LKey <- LKey[o1]
            # index matching
            QGeneLength <- GeneCalls[[m1]][PMatrix[, 1L], "Stop"] - GeneCalls[[m1]][PMatrix[, 1L], "Start"] + 1L
            SGeneLength <- GeneCalls[[m2]][PMatrix[, 2L], "Stop"] - GeneCalls[[m2]][PMatrix[, 2L], "Start"] + 1L
            # InteriorMissQuery <- QGeneLength - (ExactOverLap + ExteriorMissQuery)
            # InteriorMissSubject <- SGeneLength - (ExactOverLap + ExteriorMissSubject)
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
            TotalKmers <- SyntenyLinks[[m1, m2]][, 11L]
            MaxKmer <- SyntenyLinks[[m1, m2]][, 10L]
            # ExteriorMissQuery <- SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L]
            # ExteriorMissSubject <- SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L]
            # InteriorMissQuery <- QGeneLength - (ExactOverLap + ExteriorMissQuery)
            # InteriorMissSubject <- SGeneLength - (ExactOverLap + ExteriorMissSubject)
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
          TotalKmers <- SyntenyLinks[[m1, m2]][, 11L]
          MaxKmer <- SyntenyLinks[[m1, m2]][, 10L]
          # ExteriorMissQuery <- SyntenyLinks[[m2, m1]][, 1L] + SyntenyLinks[[m2, m1]][, 2L]
          # ExteriorMissSubject <- SyntenyLinks[[m2, m1]][, 3L] + SyntenyLinks[[m2, m1]][, 4L]
          # InteriorMissQuery <- QGeneLength - (ExactOverLap + ExteriorMissQuery)
          # InteriorMissSubject <- SGeneLength - (ExactOverLap + ExteriorMissSubject)
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
        
        # PresQI <- unique(IMatrix[, 1L])
        # PresSI <- unique(IMatrix[, 2L])
        # QuerySeqs <- vector(mode = "list",
        #                     length = length(PresQI))
        # SubjectSeqs <- vector(mode = "list",
        #                       length = length(PresSI))
        # PresentQRanges <- GeneCalls[[m1]]$Range[PMatrix[, 1L]]
        # PresentSRanges <- GeneCalls[[m2]]$Range[PMatrix[, 2L]]
        # 
        # for (m3 in seq_along(PresQI)) {
        #   CurrentRanges <- PresentQRanges[IMatrix[, 1L] == PresQI[m3]]
        #   QuerySeqs[[m3]] <- extractAt(x = rep(Genomes[[m1]][PresQI[m3]],
        #                                        length(CurrentRanges)),
        #                                at = unname(CurrentRanges))
        #   QuerySeqs[[m3]] <- DNAStringSet(sapply(QuerySeqs[[m3]],
        #                                          function(x) unlist(x)))
        #   names(QuerySeqs[[m3]]) <- paste(names(GeneCalls)[m1],
        #                                   PresQI[m3],
        #                                   PMatrix[, 1L][IMatrix[, 1L] == PresQI[m3]],
        #                                   sep = "_")
        # }
        # QuerySeqs <- do.call(c,
        #                      QuerySeqs)
        # for (m3 in seq_along(PresSI)) {
        #   CurrentRanges <- PresentSRanges[IMatrix[, 2L] == PresSI[m3]]
        #   SubjectSeqs[[m3]] <- extractAt(x = rep(Genomes[[m2]][PresSI[m3]],
        #                                          length(CurrentRanges)),
        #                                  at = unname(CurrentRanges))
        #   SubjectSeqs[[m3]] <- DNAStringSet(sapply(SubjectSeqs[[m3]],
        #                                            function(x) unlist(x)))
        #   names(SubjectSeqs[[m3]]) <- paste(names(GeneCalls)[m2],
        #                                     PresSI[m3],
        #                                     PMatrix[, 2L][IMatrix[, 2L] == PresSI[m3]],
        #                                     sep = "_")
        # }
        # SubjectSeqs <- do.call(c,
        #                        SubjectSeqs)
        # # return SubjectSeqs to the original order
        # o <- match(x = paste(names(GeneCalls)[m2],
        #                      IMatrix[, 2L],
        #                      PMatrix[, 2L],
        #                      sep = "_"),
        #            table = names(SubjectSeqs))
        # SubjectSeqs <- SubjectSeqs[o]
        # if (length(QuerySeqs) != length(SubjectSeqs)) {
        #   stop("Extracted Seqs have differing lengths.")
        # }
        
        # ask if features were pulled initially
        # if not grab them
        if (is.null(Features[[m1]])) {
          # hypothetically, we shouldn't go here unless only one 
          # set of features is pulled in at a time,
          # in which case only one set of sequences will be held in memory at a time
          # for every cell in the matrix
          Genome <- SearchDB(dbFile = DBPATH,
                             identifier = names(GeneCalls[m1]),
                             nameBy = "description",
                             verbose = FALSE)
          PresentIndices <- unique(GeneCalls[[m1]]$Index)
          if (length(PresentIndices) > 1L) {
            # many indices, loop through present indices and extract
            # slam together at the end
            Features[[m1]] <- vector(mode = "list",
                                        length = length(PresentIndices))
            for (m3 in seq_along(PresentIndices)) {
              ph <- GeneCalls[[m1]]$Index == PresentIndices[m3]
              Features[[m1]][[m3]] <- extractAt(x = rep(Genome[PresentIndices[m3]],
                                                        sum(ph)),
                                                at = unname(GeneCalls[[m1]]$Range[ph]))
              Features[[m1]][[m3]] <- DNAStringSet(sapply(Features[[m1]][[m3]],
                                                          function(x) unlist(x)))
              FlipMe <- GeneCalls[[m1]]$Strand[ph] == 1L
              if (any(FlipMe)) {
                Features[[m1]][[m3]][FlipMe] <- reverseComplement(Features[[m1]][[m3]][FlipMe])
              }
            }
            Features[[m1]] <- do.call(c,
                                      Features[[m1]])
          } else {
            # only 1 index present in gene calls
            Features[[m1]] <- extractAt(x = rep(Genome[PresentIndices],
                                                nrow(GeneCalls[[m1]])),
                                        at = unname(GeneCalls[[m1]]$Range))
            Features[[m1]] <- DNAStringSet(sapply(Features[[m1]],
                                                  function(x) unlist(x)))
            FlipMe <- GeneCalls[[m1]]$Strand == 1L
            if (any(FlipMe)) {
              Features[[m1]][FlipMe] <- reverseComplement(Features[[m1]][FlipMe])
            }
            
          }
          names(Features[[m1]]) <- paste(rep(names(GeneCalls)[m1], length(Features[[m1]])),
                                         GeneCalls[[m1]]$Index,
                                         seq(length(Features[[m1]])),
                                         sep = "_")
          QuerySeqs <- Features[[m1]][PMatrix[, 1L]]
          Features[m1] <- list(NULL)
        } else {
          QuerySeqs <- Features[[m1]][PMatrix[, 1L]]
        }
        
        if (is.null(Features[[m2]])) {
          Genome <- SearchDB(dbFile = DBPATH,
                             identifier = names(GeneCalls[m2]),
                             nameBy = "description",
                             verbose = FALSE)
          PresentIndices <- unique(GeneCalls[[m2]]$Index)
          if (length(PresentIndices) > 1L) {
            # many indices, loop through present indices and extract
            # slam together at the end
            Features[[m2]] <- vector(mode = "list",
                                     length = length(PresentIndices))
            for (m3 in seq_along(PresentIndices)) {
              ph <- GeneCalls[[m2]]$Index == PresentIndices[m3]
              Features[[m2]][[m3]] <- extractAt(x = rep(Genome[PresentIndices[m3]],
                                                        sum(ph)),
                                                at = unname(GeneCalls[[m2]]$Range[ph]))
              Features[[m2]][[m3]] <- DNAStringSet(sapply(Features[[m2]][[m3]],
                                                          function(x) unlist(x)))
              FlipMe <- GeneCalls[[m2]]$Strand[ph] == 1L
              if (any(FlipMe)) {
                Features[[m2]][[m3]][FlipMe] <- reverseComplement(Features[[m2]][[m3]][FlipMe])
              }
            }
            Features[[m2]] <- do.call(c,
                                      Features[[m2]])
          } else {
            # only 1 index present in gene calls
            Features[[m2]] <- extractAt(x = rep(Genome[PresentIndices],
                                                nrow(GeneCalls[[m2]])),
                                        at = unname(GeneCalls[[m2]]$Range))
            Features[[m2]] <- DNAStringSet(sapply(Features[[m2]],
                                                  function(x) unlist(x)))
            FlipMe <- GeneCalls[[m2]]$Strand == 1L
            if (any(FlipMe)) {
              Features[[m2]][FlipMe] <- reverseComplement(Features[[m2]][FlipMe])
            }
            
          }
          names(Features[[m2]]) <- paste(rep(names(GeneCalls)[m2], length(Features[[m2]])),
                                         GeneCalls[[m2]]$Index,
                                         seq(length(Features[[m2]])),
                                         sep = "_")
          SubjectSeqs <- Features[[m2]][PMatrix[, 2L]]
          PresentFeatures <- unname(sapply(Features,
                                           function(x) !is.null(x),
                                           simplify = TRUE,
                                           USE.NAMES = FALSE))
          if (m1 > 1L) {
            # no replacement of features in memory when m1 == 1L
            # if any of the features that will no longer be visited are not NULL
            # NULL one of those as opposed to what was just loaded in
            if (any(PresentFeatures[1L:(m1 - 1L)])) {
              SubjectSeqs <- Features[[m2]][PMatrix[, 2L]]
              Features[min(which(PresentFeatures))] <- list(NULL)
            } else {
              SubjectSeqs <- Features[[m2]][PMatrix[, 2L]]
              Features[m2] <- list(NULL)
            }
          } else {
            SubjectSeqs <- Features[[m2]][PMatrix[, 2L]]
            Features[m2] <- list(NULL)
          }
          
        } else {
          SubjectSeqs <- Features[[m2]][PMatrix[, 2L]]
        }
        # reverse complement seqs where necessary
        # not necessary anymore, rC was performed at pull
        # QuerySeqs[QGeneStrand == 1L] <- reverseComplement(QuerySeqs[QGeneStrand == 1L])
        # SubjectSeqs[SGeneStrand == 1L] <- reverseComplement(SubjectSeqs[SGeneStrand == 1L])
        
        # return(list(Features,
        #             SubjectSeqs,
        #             QuerySeqs))
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
              if ("APArgs" %in% ls()) {
                CurrentAPArgs <- c(list("pattern" = QuerySeqs[m3],
                                        "subject" = SubjectSeqs[m3]),
                                   APArgs)
                
              } else {
                CurrentAPArgs <- list("pattern" = QuerySeqs[m3],
                                      "subject" = SubjectSeqs[m3])
              }
              ph01 <- do.call(what = AlignProfiles,
                              args = CurrentAPArgs)
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
            }
          } else {
            # perform amino acid alignments where possible
            for (m3 in seq_along(SubjectSeqs)) {
              if (QGeneCoding[m3] &
                  SGeneCoding[m3] &
                  width(SubjectSeqs[m3]) %% 3 == 0 &
                  width(QuerySeqs[m3]) %% 3 == 0 &
                  !grepl(pattern = "[^ATCG]",
                         x = QuerySeqs[m3]) &
                  !grepl(pattern = "[^ATCG]",
                         x = SubjectSeqs[m3])) {
                Atype[m3] <- "AA"
                if ("APArgs" %in% ls()) {
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
                  CurrentAPArgs <- c(list("pattern" = translate(QuerySeqs[m3],
                                                                genetic.code = gC1),
                                          "subject" = translate(SubjectSeqs[m3],
                                                                genetic.code = gC2)),
                                     APArgs)
                  
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
                  CurrentAPArgs <- list("pattern" = translate(QuerySeqs[m3],
                                                              genetic.code = gC1),
                                        "subject" = translate(SubjectSeqs[m3],
                                                              genetic.code = gC2))
                }
                # return(CurrentATArgs)
                ph01 <- do.call(what = AlignProfiles,
                                args = CurrentAPArgs)
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
              } else {
                Atype[m3] <- "NT"
                if ("APArgs" %in% ls()) {
                  CurrentAPArgs <- c(list("pattern" = QuerySeqs[m3],
                                          "subject" = SubjectSeqs[m3]),
                                     APArgs)
                  
                } else {
                  CurrentAPArgs <- list("pattern" = QuerySeqs[m3],
                                        "subject" = SubjectSeqs[m3])
                }
                ph01 <- do.call(what = AlignProfiles,
                                args = CurrentAPArgs)
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
                                    "TotalKmers" = TotalKmers,
                                    "MaxKmer" = MaxKmer,
                                    "Concensus" = diff2,
                                    "p1FeatureLength" = QGeneLength,
                                    "p2FeatureLength" = SGeneLength,
                                    "Adjacent" = RKey + LKey,
                                    "TetDist" = NucDist,
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
                                    "Concensus" = diff2,
                                    "TotalKmers" = TotalKmers,
                                    "MaxKmer" = MaxKmer,
                                    "p1FeatureLength" = QGeneLength,
                                    "p2FeatureLength" = SGeneLength,
                                    "Adjacent" = RKey + LKey,
                                    "TetDist" = NucDist,
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
  
  
  
  
  
  
  
  
  
  