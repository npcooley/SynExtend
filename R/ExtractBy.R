###### -- ExtractBy -----------------------------------------------------------
# author: nicholas cooley
# email: npc19@pitt.edu / npcooley@gmail.com
# extract stringsets from an assembly based on the gene calls
# OR
# extract stringsets from a series of assemblies based on a pairsummaries object
# OR
# extract stringsets from a series of assemblies based on COGs / pairsummaries object

ExtractBy <- function(x,
                      y,
                      z,
                      Verbose = FALSE) {
  
  # overhead things:
  # z does not need to be supplied by the user but under the current scheme must be evaluated
  if (missing(z)) {
    z <- NA_integer_
  }
  # combinations of x, y, and z that are valid:
  # x as a DFrame of genecalls
  # y as a DNAStringSet
  # 
  # x as a PairSummaries object
  # y as a character string supplying the location of a SQLite DB
  # 
  # x as a PairSummaries object
  # y as a character string supplying the location of a SQLite DB
  # z as a list of COGs in the format:
  # [1] "GenomeID_ContigID_FeatureID" "GenomeID_ContigID_FeatureID" etc...
  if (missing(y) |
      missing(x)) {
    stop("x and y must be supplied.")
  }
  
  if (is(object = x,
         class2 = "DFrame") &
      is(object = y,
         class2 = "DNAStringSet")) {
    # method 1: pull sequences from a dna stringset based on a gene calls DFrame
    # x shot, where x is the number of contigs in the assembly that contain gene calls
    # should be relatively fast for a typical bacterial genome
    # grab all nucleotide sequences from the assembly
    # match based on contig names
    
    if (Verbose) {
      TimeStart <- Sys.time()
    }
    AssemblyRef <- unlist(regmatches(x = names(y),
                                     m = gregexpr(pattern = "^[^ ]+",
                                                  text = names(y))))
    GCRef <- unique(x$Contig)
    
    seqs1 <- vector(mode = "list",
                    length = length(GCRef))
    for (m1 in seq_along(GCRef)) {
      # identify the string to access
      w1 <- which(AssemblyRef == GCRef[m1])
      # grab the IRangesList for the current contig
      z1 <- unname(x$Range[x$Contig == GCRef[m1]])
      # create a collapse map
      z2 <- lengths(z1)
      # unlist to an IRanges object
      z1 <- unlist(z1,
                   recursive = FALSE)
      seqs1[[m1]] <- extractAt(x = y[w1][[1L]],
                               at = z1)
      # figure out where collapses are necessary
      collapsecount <- 0L
      w2 <- which(z2 > 1L)
      # if collapsing is necessary, collapse on specified positions
      if (length(w2) > 0L) {
        # for a collapse operation of 2, we need to remove 1 position
        # for an operation of 3, we need to remove 2 positions
        # for an operation of 4, we need to remove 3 positions
        # etc
        removevector <- vector(mode = "integer",
                               length = sum(z2[w2]) - length(w2))
        for (m2 in w2) {
          seqs1[[m1]][[m2 + collapsecount]] <- unlist(seqs1[[m1]][m2:(m2 + z2[m2] - 1L) + collapsecount])
          removevector[(collapsecount + 1L):(collapsecount + z2[m2] - 1L)] <- (m2 + 1L):(m2 + z2[m2] - 1L) + collapsecount
          collapsecount <- collapsecount + z2[m2] - 1L
        } # end m2 loop
        seqs1[[m1]][removevector] <- NULL
      } # end check on collapse operation
      names(seqs1[[m1]]) <- x$ID[x$Contig == GCRef[m1]]
    } # end m1 loop
    seqs1 <- do.call(c,
                     seqs1)
    
    # flip sequences into the forward direction where applicable
    flipseqs <- x$Strand
    if (any(flipseqs == 1L)) {
      seqs1[flipseqs == 1L] <- reverseComplement(seqs1[flipseqs == 1L])
    }
    
    # end method 1
    return(seqs1)
  } else if (is(object = x,
                class2 = "PairSummaries") &
             is(object = y,
                class2 = "character") &
             !is(object = z,
                 class2 = "list")) {
    # method 2: pull sequences from a SQLite DB as a list of paired sequences
    # order n, where n is the number of total contigs that contain both a genecall
    # and an identified pair partner
    
    if (Verbose) {
      TimeStart <- Sys.time()
    }
    # grab unique gene IDs
    u1 <- unique(c(unique(x$p1),
                   unique(x$p2)))
    # build a map to split out the seqs post extraction
    u2 <- matrix(data = c(match(x = x$p1,
                                table = u1),
                          match(x = x$p2,
                                table = u1)),
                 nrow = nrow(x))
    # build a map of which strings are being extracted from
    u3 <- do.call(rbind,
                  strsplit(x = u1,
                           split = "_",
                           fixed = TRUE))
    # all extractions
    u3 <- matrix(data = as.integer(u3),
                 nrow = nrow(u3))
    # all strings to extract
    u4 <- unique(u3[, c(1L, 2L)])
    # all identifiers that are necessary:
    u5 <- unique(u4[, 1L])
    L01 <- nrow(u4)
    L02 <- 0L
    if (Verbose) {
      # pBar target length:
      pBar <- txtProgressBar(style = 1L)
      cat("\nExtracting Sequences:\n")
    }
    Res1 <- vector(mode = "list",
                   length = L01)
    # return(list(u1,
    #             u2,
    #             u3,
    #             u4,
    #             u5,
    #             x))
    # loop through each extraction in the first level
    # in the second level loop through loop through the contigs
    # use name matching to grab the correct contig
    # THIS FORCES THE ASSUMPTION THAT A STANDARD NCBI NAMING CONVENTION HAS BEEN FOLLOWED
    # When the nucleotide overlap object was built the gene calls objects was force-reordered
    GC01 <- attr(x = x,
                 which = "GeneCalls")
    GCN <- names(GC01)
    for (m1 in seq_along(u5)) {
      # ID contigs
      u6 <- u4[u4[, 1L] == u5[m1], 2L]
      # grab assembly
      CurrentAssembly <- SearchDB(dbFile = y,
                                  identifier = as.character(u5[m1]),
                                  verbose = FALSE,
                                  nameBy = "description")
      AssemblyRef <- unlist(regmatches(x = names(CurrentAssembly),
                                       m = gregexpr(pattern = "^[^ ]+",
                                                    text = names(CurrentAssembly))))
      CurrentGeneCalls <- GC01[GCN == as.character(u5[m1])][[1]]
      for (m2 in seq_along(u6)) {
        # identify the string to access
        GCRef <- CurrentGeneCalls$Contig[(CurrentGeneCalls$Index == u6[m2])][1L]
        w1 <- which(AssemblyRef == GCRef)
        # grab the IRangesList for the current contig
        u7 <- u3[u3[, 1L] == u5[m1] & u3[, 2L] == u6[m2], 3L]
        # print(c(u5[m1], u6[m2]))
        z1 <- unname(CurrentGeneCalls$Range[u7])
        # create a collapse map
        z2 <- lengths(z1)
        # unlist to an IRanges object
        z1 <- unlist(z1,
                     recursive = FALSE)
        # return(list(CurrentAssembly,
        #             w1,
        #             z1,
        #             AssemblyRef,
        #             GCRef,
        #             CurrentGeneCalls,
        #             u5,
        #             u6,
        #             u7,
        #             m1,
        #             m2))
        seqs1 <- extractAt(x = CurrentAssembly[w1][[1L]],
                           at = z1)
        # figure out where collapses are necessary
        collapsecount <- 0L
        w2 <- which(z2 > 1L)
        # if collapsing is necessary, collapse on specified positions
        if (length(w2) > 0L) {
          # for a collapse operation of 2, we need to remove 1 position
          # for an operation of 3, we need to remove 2 positions
          # for an operation of 4, we need to remove 3 positions
          # etc
          removevector <- vector(mode = "integer",
                                 length = sum(z2[w2]) - length(w2))
          for (m3 in w2) {
            seqs1[[m3 + collapsecount]] <- unlist(seqs1[m3:(m3 + z2[m3] - 1L) + collapsecount])
            removevector[(collapsecount + 1L):(collapsecount + z2[m3] - 1L)] <- (m3 + 1L):(m3 + z2[m3] - 1L) + collapsecount
            collapsecount <- collapsecount + z2[m3] - 1L
          } # end m3 loop through collapse situations
          seqs1[removevector] <- NULL
        } # end check on collapse operation
        
        # flip sequences into the forward direction where applicable
        flipseqs <- CurrentGeneCalls$Strand[u7]
        if (any(flipseqs == 1L)) {
          seqs1[flipseqs == 1L] <- reverseComplement(seqs1[flipseqs == 1L])
        }
        L02 <- L02 + 1L
        if (Verbose) {
          setTxtProgressBar(pb = pBar,
                            value = L02 / L01)
        } # end verbose logical
        names(seqs1) <- paste(u5[m1],
                              u6[m2],
                              u7,
                              sep = "_")
        Res1[[L02]] <- seqs1
      } # end m2 loop through contigs
    } # end m1 loop through assembly identifiers
    Res1 <- do.call(c,
                    Res1)
    if (Verbose) {
      close(pBar)
      cat("\nArranging Sequences:\n")
    }
    Res2 <- vector(mode = "list",
                   length = nrow(x))
    if (Verbose) {
      L01 <- length(Res2)
      pBar <- txtProgressBar(style = 1L)
    }
    # return(list(Res1,
    #             u2))
    for (m1 in seq_along(Res2)) {
      # call remove gaps to 
      Res2[[m1]] <- RemoveGaps(Res1[names(Res1) %in% u1[u2[m1, ]]])
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / L01)
      }
      
    }
    if (Verbose) {
      close(pBar)
      cat("\n")
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    return(Res2)
    
  } else if (is(object = x,
                class2 = "PairSummaries") &
             is(object = y,
                class2 = "character") &
             is(object = z,
                 class2 = "list")) {
    # method 3: pull sequences from a SQLite DB as a list of identified single-linkage
    # COGs
    # order n, where n is the number of total contigs that contain both a genecall
    # and an identified COG member
    
    if (Verbose) {
      TimeStart <- Sys.time()
    }
    # grab unique gene IDs
    u1 <- unique(unlist(z))
    # build a map to split out the seqs post extraction
    u2 <- unname(sapply(z,
                        function(x) {
                          match(x = x,
                                table = u1)
                        },
                        simplify = FALSE,
                        USE.NAMES = FALSE))
    # build a map of which strings are being extracted from
    u3 <- do.call(rbind,
                  strsplit(x = u1,
                           split = "_",
                           fixed = TRUE))
    # all extractions
    u3 <- matrix(data = as.integer(u3),
                 nrow = nrow(u3))
    # all strings to extract
    u4 <- unique(u3[, c(1L, 2L)])
    # all identifiers that are necessary:
    u5 <- unique(u4[, 1L])
    L01 <- nrow(u4)
    L02 <- 0L
    if (Verbose) {
      # pBar target length:
      pBar <- txtProgressBar(style = 1L)
      cat("\nExtracting Sequences:\n")
    }
    Res1 <- vector(mode = "list",
                   length = L01)
    # return(list(u1,
    #             u2,
    #             u3,
    #             u4,
    #             u5,
    #             z))
    # loop through each extraction in the first level
    # in the second level loop through loop through the contigs
    # use name matching to grab the correct contig
    # THIS FORCES THE ASSUMPTION THAT A STANDARD NCBI NAMING CONVENTION HAS BEEN FOLLOWED
    # When the nucleotide overlap object was built the gene calls objects was force-reordered
    GC01 <- attr(x = x,
                 which = "GeneCalls")
    GCN <- names(GC01)
    for (m1 in seq_along(u5)) {
      # ID contigs
      u6 <- u4[u4[, 1L] == u5[m1], 2L]
      # grab assembly
      CurrentAssembly <- SearchDB(dbFile = y,
                                  identifier = as.character(u5[m1]),
                                  verbose = FALSE,
                                  nameBy = "description")
      AssemblyRef <- unlist(regmatches(x = names(CurrentAssembly),
                                       m = gregexpr(pattern = "^[^ ]+",
                                                    text = names(CurrentAssembly))))
      CurrentGeneCalls <- GC01[GCN == as.character(u5[m1])][[1]]
      for (m2 in seq_along(u6)) {
        # identify the string to access
        GCRef <- CurrentGeneCalls$Contig[(CurrentGeneCalls$Index == u6[m2])][1L]
        w1 <- which(AssemblyRef == GCRef)
        # grab the IRangesList for the current contig
        u7 <- u3[u3[, 1L] == u5[m1] & u3[, 2L] == u6[m2], 3L]
        # print(c(u5[m1], u6[m2]))
        z1 <- unname(CurrentGeneCalls$Range[u7])
        # create a collapse map
        z2 <- lengths(z1)
        # unlist to an IRanges object
        z1 <- unlist(z1,
                     recursive = FALSE)
        # return(list(CurrentAssembly,
        #             w1,
        #             z1,
        #             AssemblyRef,
        #             GCRef,
        #             CurrentGeneCalls,
        #             u5,
        #             u6,
        #             u7,
        #             m1,
        #             m2))
        seqs1 <- extractAt(x = CurrentAssembly[w1][[1L]],
                           at = z1)
        # figure out where collapses are necessary
        collapsecount <- 0L
        w2 <- which(z2 > 1L)
        # if collapsing is necessary, collapse on specified positions
        if (length(w2) > 0L) {
          # for a collapse operation of 2, we need to remove 1 position
          # for an operation of 3, we need to remove 2 positions
          # for an operation of 4, we need to remove 3 positions
          # etc
          removevector <- vector(mode = "integer",
                                 length = sum(z2[w2]) - length(w2))
          for (m3 in w2) {
            seqs1[[m3 + collapsecount]] <- unlist(seqs1[m3:(m3 + z2[m3] - 1L) + collapsecount])
            removevector[(collapsecount + 1L):(collapsecount + z2[m3] - 1L)] <- (m3 + 1L):(m3 + z2[m3] - 1L) + collapsecount
            collapsecount <- collapsecount + z2[m3] - 1L
          } # end m3 loop through collapse situations
          seqs1[removevector] <- NULL
        } # end check on collapse operation
        
        # flip sequences into the forward direction where applicable
        flipseqs <- CurrentGeneCalls$Strand[u7]
        if (any(flipseqs == 1L)) {
          seqs1[flipseqs == 1L] <- reverseComplement(seqs1[flipseqs == 1L])
        }
        L02 <- L02 + 1L
        if (Verbose) {
          setTxtProgressBar(pb = pBar,
                            value = L02 / L01)
        } # end verbose logical
        names(seqs1) <- paste(u5[m1],
                              u6[m2],
                              u7,
                              sep = "_")
        Res1[[L02]] <- seqs1
      } # end m2 loop through contigs
    } # end m1 loop through assembly identifiers
    Res1 <- do.call(c,
                    Res1)
    if (Verbose) {
      close(pBar)
      cat("\nArranging Sequences:\n")
    }
    Res2 <- vector(mode = "list",
                   length = length(z))
    if (Verbose) {
      L01 <- length(Res2)
      pBar <- txtProgressBar(style = 1L)
    }
    # return(list(Res1,
    #             u2))
    for (m1 in seq_along(Res2)) {
      # call remove gaps to 
      Res2[[m1]] <- RemoveGaps(Res1[names(Res1) %in% z[[m1]]])
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / L01)
      }
      
    }
    if (Verbose) {
      close(pBar)
      cat("\n")
      TimeEnd <- Sys.time()
      print(TimeEnd - TimeStart)
    }
    return(Res2)
    
  } else {
    stop("No method available for combination of supplied objects.")
  }
  
}


