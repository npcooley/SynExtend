# Author: Nicholas Cooley
# Maintainer: Nicholas Cooley
# Contact: npc19@pitt.edu
# # to do list:
# add attributes to clusters -- a 3 column dataframe with logicals
# indicating whether it is marked as translatable in the GeneCalls
# and whether it is actually width %% 3L == 0L
# and annotations if they are present in the GeneCalls
# OR the model from FindGenes

ExtractBy <- function(x,
                      y = NULL,
                      DBPATH,
                      Method = "all",
                      DefaultTranslationTable = "11",
                      Translate = TRUE,
                      Storage = 1,
                      Verbose = FALSE) {
  
  ###### -- overhead checking -------------------------------------------------
  # Translation tables not implemented yet
  if (Verbose) {
    pBar <- txtProgressBar(style = 1L)
    TimeStart <- Sys.time()
  }
  
  if (!(Method %in% c("pairs",
                      "clusters",
                      "all",
                      "DataFrame",
                      "columns"))) {
    stop ("Please select a known Method.")
  }
  if (Method == "clusters" &
      is.null(y)) {
    stop ("Please supply valid clusters.")
  }
  if (Storage < 0) {
    stop ("Storage must be at least zero")
  } else {
    # convert storage to bytes
    Storage <- Storage * 1e9 
  }
  
  ###### -- methods are all separate ------------------------------------------
  
  # all - simplest - pull all seqs in PairSummaries and slam them into a single stringset
  # clusters - similar to all, but subset Pairs by a cluster list
  
  if (Method == "all") {
    
    if (Verbose) {
      cat("\nPreparing overhead data.\n")
    }
    
    GeneCalls <- attributes(x)$GeneCalls
    # create the maps necessary to extract seqs as efficiently as possible
    PullTable <- do.call(rbind,
                         strsplit(c(x$p1, x$p2),
                                  split = "_",
                                  fixed = TRUE))
    PullTable <- unique(PullTable)
    PullTable <- matrix(data = as.integer(PullTable),
                        nrow = nrow(PullTable))
    G <- unique(PullTable[, 1L])
    L <- length(G)
    
    # create a list to fill with integers for contigs to pull from
    # and the range positons to pull with
    CList <- LList <- RList <- IList <- SList <- vector(mode = "list",
                                                        length = L)
    w <- vector(mode = "integer",
                length = L)
    for (m1 in seq_along(CList)) {
      # create overhead lookups to pull genes succinctly
      CList[[m1]] <- PullTable[PullTable[, 1L] == G[m1], 2L]
      IList[[m1]] <- PullTable[PullTable[, 1L] == G[m1], 3L]
      w[m1] <- which(names(GeneCalls) == as.character(G[m1]))
      LList[[m1]] <- GeneCalls[[w[m1]]][IList[[m1]], "Coding"]
      RList[[m1]] <- GeneCalls[[w[m1]]][IList[[m1]], "Range"]
      SList[[m1]] <- GeneCalls[[w[m1]]][IList[[m1]], "Strand"]
      
    }
    
    # something else here
    # yay new setup
    # return(list(CList))
    
    
    # create a list to dump stringsets into
    GList <- vector(mode = "list",
                    length = L)
    
    Count <- 1L
    while (object.size(GList) < Storage &
           Count <= L) {
      
      GList[[Count]] <- SearchDB(dbFile = DBPATH,
                              identifier = as.character(G[Count]),
                              nameBy = "description",
                              verbose = FALSE)
      
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
    
    # return(list(GList,
    #             w,
    #             CList,
    #             IList,
    #             RList,
    #             LList))
    
    # initialize res as a list
    # fill with stringsets that will eventually be slammed together
    Res <- vector(mode = "list",
                  length = L)
    
    for (m1 in seq_along(G)) {
      
      # if seqs are in the genome list grab them,
      # if not pull from DB
      if (!is.null(GList[[m1]])) {
        s1 <- GList[[m1]]
      } else {
        s1 <- SearchDB(dbFile = DBPATH,
                       identifier = as.character(G[m1]),
                       nameBy = "description",
                       verbose = FALSE)
      }
      CurrentContigs <- unique(CList[[m1]])
      if (length(CurrentContigs) == 1L) {
        # only one contig present for s1 -- the easiest case
        
        # return(list(s1[CurrentContigs],
        #             length(IList[[m1]]),
        #             RList[[m1]],
        #             CurrentContigs,
        #             IList,
        #             RList,
        #             m1))
        
        Res[[m1]] <- extractAt(x = rep(s1[CurrentContigs],
                                       length(IList[[m1]])),
                               at = RList[[m1]])
        Res[[m1]] <- DNAStringSet(sapply(Res[[m1]],
                                         function(x) unlist(x),
                                         simplify = FALSE,
                                         USE.NAMES = FALSE))
        
        FlipMe <- SList[[m1]] == 1L
        if (any(FlipMe)) {
          Res[[m1]][FlipMe] <- reverseComplement(Res[[m1]][FlipMe])
        }
        # set seq names
        names(Res[[m1]]) <- paste(rep(w[m1],
                                      length(IList[[m1]])),
                                  CList[[m1]],
                                  IList[[m1]],
                                  sep = "_")
      } else {
        # multiple contigs -- in this case seqs don't really need to be ordered
        # as long as they have the correct names
        Res[[m1]] <- vector(mode = "list",
                            length = length(CurrentContigs))
        for (m2 in seq_along(Res[[m1]])) {
          CPos <- CList[[m1]] == CurrentContigs[m2]
          CSum <- sum(CPos)
          Res[[m1]][[m2]] <- extractAt(x = rep(s1[CurrentContigs[m2]],
                                               CSum),
                                       at = RList[[m1]][CPos])
          Res[[m1]][[m2]] <- DNAStringSet(sapply(Res[[m1]][[m2]],
                                                 function(x) unlist(x),
                                                 simplify = FALSE,
                                                 USE.NAMES = FALSE))
          FlipMe <- SList[[m1]][CPos] == 1L
          if (any(FlipMe)) {
            Res[[m1]][[m2]][FlipMe] <- reverseComplement(Res[[m1]][[m2]][FlipMe])
          }
          # set names
          names(Res[[m1]][[m2]]) <- paste(rep(w[m1],
                                              CSum),
                                          CList[[m1]][CPos],
                                          IList[[m1]][CPos],
                                          sep = "_")
        }
        
        Res[[m1]] <- do.call(c,
                             Res[[m1]])
      }
      
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / L)
      }
    }
    
    Res <- do.call(c,
                   Res)
    
    # slam these lists together into a single string set
    if (all(unlist(LList)) &
        all(width(Res) %% 3 == 0) &
        Translate) {
      # all are translatable
      if (any(grepl(pattern = "[^ATCG]",
                    x = Res))) {
        Res <- translate(x = Res,
                         if.fuzzy.codon = "solve")
      } else {
        Res <- translate(x = Res,
                         if.fuzzy.codon = "error")
      }
      
    } else {
      # at least one is not translatable -- leave as DNAStringSet
    }
    
    # end method == all
  } else if (Method == "clusters") {
    # overhead for cluster checking
    # make sure all IDs in the clusters are present in the PairSummaries object
    # and vice-versa
    AllPresentClusters <- unlist(y)
    AllPresentPartners <- unique(c(unique(x$p1),
                                   unique(x$p2)))
    
    if (!all(AllPresentClusters %in% AllPresentPartners) |
        !all(AllPresentPartners %in% AllPresentClusters)) {
      stop ("Discrepency in identifiers.")
    }
    
    # build GList ahead of cluster parsing from GeneCalls attribute
    
    if (Verbose) {
      cat("\nPreparing overhead data.\n")
    }
    
    GeneCalls <- attributes(x)$GeneCalls
    G <- as.integer(names(GeneCalls))
    L <- length(G)
    
    # create a list to dump stringsets into
    GList <- vector(mode = "list",
                    length = L)
    
    Count <- 1L
    while (object.size(GList) < Storage &
           Count <= L) {
      
      GList[[Count]] <- SearchDB(dbFile = DBPATH,
                                 identifier = as.character(G[Count]),
                                 nameBy = "description",
                                 verbose = FALSE)
      
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
    
    Res <- vector(mode = "list",
                  length = length(y))
    
    for (m1 in seq_along(y)) {
      PullTable <- do.call(rbind,
                           strsplit(x = y[[m1]],
                                    split = "_",
                                    fixed = TRUE))
      PullTable <- matrix(data = as.integer(PullTable),
                          nrow = nrow(PullTable))
      # can overwrite G and L ... are now relative to the current Pull Table
      G <- unique(PullTable[, 1L])
      L <- length(G)
      
      # create a list to fill with integers for contigs to pull from
      # and the range positons to pull with
      CList <- LList <- RList <- IList <- SList <- vector(mode = "list",
                                                          length = L)
      w <- vector(mode = "integer",
                  length = L)
      for (m2 in seq_along(CList)) {
        # create overhead lookups to pull genes succinctly
        CList[[m2]] <- PullTable[PullTable[, 1L] == G[m2], 2L]
        IList[[m2]] <- PullTable[PullTable[, 1L] == G[m2], 3L]
        w[m2] <- which(names(GeneCalls) == as.character(G[m2]))
        LList[[m2]] <- GeneCalls[[w[m2]]][IList[[m2]], "Coding"]
        RList[[m2]] <- GeneCalls[[w[m2]]][IList[[m2]], "Range"]
        SList[[m2]] <- GeneCalls[[w[m2]]][IList[[m2]], "Strand"]
        
      }
      
      # pull table and maps initialized
      Res[[m1]] <- vector(mode = "list",
                          length = L)
      
      for (m2 in seq_along(G)) {
        
        # if seqs are in the genome list grab them,
        # if not pull from DB
        if (!is.null(GList[[w[m2]]])) {
          s1 <- GList[[w[m2]]]
        } else {
          s1 <- SearchDB(dbFile = DBPATH,
                         identifier = as.character(G[m2]),
                         nameBy = "description",
                         verbose = FALSE)
        }
        CurrentContigs <- unique(CList[[m2]])
        if (length(CurrentContigs) == 1L) {
          # only one contig present for s1 -- the easiest case
          
          # return(list(s1[CurrentContigs],
          #             length(IList[[m2]]),
          #             RList[[m2]],
          #             CurrentContigs,
          #             IList,
          #             RList,
          #             m2))
          
          Res[[m1]][[m2]] <- extractAt(x = rep(s1[CurrentContigs],
                                               length(IList[[m2]])),
                                       at = RList[[m2]])
          Res[[m1]][[m2]] <- DNAStringSet(sapply(Res[[m1]][[m2]],
                                                 function(x) unlist(x),
                                                 simplify = FALSE,
                                                 USE.NAMES = FALSE))
          
          FlipMe <- SList[[m2]] == 1L
          if (any(FlipMe)) {
            Res[[m1]][[m2]][FlipMe] <- reverseComplement(Res[[m1]][[m2]][FlipMe])
          }
          # set seq names
          names(Res[[m1]][[m2]]) <- paste(rep(w[m2],
                                              length(IList[[m2]])),
                                          CList[[m2]],
                                          IList[[m2]],
                                          sep = "_")
        } else {
          # multiple contigs -- in this case seqs don't really need to be ordered
          # as long as they have the correct names
          Res[[m1]][[m2]] <- vector(mode = "list",
                                    length = length(CurrentContigs))
          for (m3 in seq_along(Res[[m1]][[m2]])) {
            CPos <- CList[[m2]] == CurrentContigs[m3]
            CSum <- sum(CPos)
            Res[[m1]][[m2]][[m3]] <- extractAt(x = rep(s1[CurrentContigs[m3]],
                                                       CSum),
                                               at = RList[[m2]][CPos])
            Res[[m1]][[m2]][[m3]] <- DNAStringSet(sapply(Res[[m1]][[m2]][[m3]],
                                                         function(x) unlist(x),
                                                         simplify = FALSE,
                                                         USE.NAMES = FALSE))
            FlipMe <- SList[[m2]][CPos] == 1L
            if (any(FlipMe)) {
              Res[[m1]][[m2]][[m3]][FlipMe] <- reverseComplement(Res[[m1]][[m2]][[m3]][FlipMe])
            }
            # set names
            names(Res[[m1]][[m2]][[m3]]) <- paste(rep(w[m2],
                                                      CSum),
                                                  CList[[m2]][CPos],
                                                  IList[[m2]][CPos],
                                                  sep = "_")
          }
          
          Res[[m1]][[m2]] <- do.call(c,
                                     Res[[m1]][[m2]])
        }
      }
      
      Res[[m1]] <- do.call(c,
                           Res[[m1]])
      
      # slam these lists together into a single string set
      if (all(unlist(LList)) &
          all(width(Res[[m1]]) %% 3 == 0) &
          Translate) {
        # all are translatable
        if (any(grepl(pattern = "[^ATCG]",
                      x = Res[[m1]]))) {
          Res[[m1]] <- translate(x = Res[[m1]],
                                 if.fuzzy.codon = "solve")
        } else {
          Res[[m1]] <- translate(x = Res[[m1]],
                                 if.fuzzy.codon = "error")
        }
        
      } else {
        # at least one is not translatable -- leave as DNAStringSet
      }
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / length(y))
      }
    } # end m1 loop
  } else if (Method == "DataFrame") {
    stop ("Method currently not implemented.")
  } else if (Method == "pairs") {
    stop ("Method currently not implemented.")
  } else if (Method == "columns") {
    stop ("Method currently not implemented.")
  }
  
  if (Verbose) {
    cat("\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  
  return(Res)
}




