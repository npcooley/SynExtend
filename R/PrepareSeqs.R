###### -- extract a bunch of sequences from some genecalls --------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu

# given an object that has a GeneCalls attribute

PrepareSeqs <- function(SynExtendObject,
                        DataBase,
                        DefaultTranslationTable = "11",
                        Identifiers = NULL,
                        Storage = 1,
                        Verbose = FALSE) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1)
  }
  
  # check object type -- currently only PairSummaries and LinkedPairs
  if (!is(object = SynExtendObject,
          class2 = "PairSummaries") &
      !is(object = SynExtendObject,
          class2 = "LinkedPairs")) {
    stop("'SynExtendObject' must be an object of class 'PairSummaries' or 'LinkedPairs'.")
  }
  # initialize database
  if (is.character(DataBase)) {
    if (!requireNamespace(package = "RSQLite",
                          quietly = TRUE)) {
      stop("Package 'RSQLite' must be installed.")
    }
    if (!("package:RSQLite" %in% search())) {
      print("Eventually character vector access to DECIPHER DBs will be deprecated.")
      require(RSQLite, quietly = TRUE)
    }
    dbConn <- dbConnect(dbDriver("SQLite"), DataBase)
    on.exit(dbDisconnect(dbConn))
  } else {
    dbConn <- DataBase
    if (!dbIsValid(dbConn)) {
      stop("The connection has expired.")
    }
  }
  # check storage
  if (Storage < 0) {
    stop("Storage must be greater than zero.")
  } else {
    Storage <- Storage * 1e9 # conversion to gigabytes
  }
  
  # get the genecalls!
  GeneCalls <- attr(x = SynExtendObject,
                    which = "GeneCalls")
  
  ###### -- subset gene calls based on the names of the links object ----------
  if (is(SynExtendObject,
         class2 = "LinkedPairs")) {
    if (length(GeneCalls) != nrow(SynExtendObject)) {
      GeneCalls <- GeneCalls[match(x = dimnames(SynExtendObject)[[1]],
                                   table = names(GeneCalls))]
    }
  } # else GeneCalls from PairSummaries objects will be guaranteed to be correctly shaped?
  
  if (!is.null(Identifiers)) {
    if (is(object = Identifiers,
           class2 = "character")) {
      GeneCalls <- GeneCalls[names(GeneCalls) %in% Identifiers]
    } else {
      stop ("Identifiers must be supplied as characters.")
    }
  }
  if (length(GeneCalls) < 1) {
    stop("Requested subset not found.")
  }
  
  # load in genomes and ALL extracted features at the top until storage limit
  # is hit
  if (Verbose) {
    cat("Preparing overhead data.\n")
  }
  # load in structure matrices once for PredictHEC
  # these hardcoded values come from Erik ... AlignProfiles maybe?
  MAT1 <- get(data("HEC_MI1",
                   package = "DECIPHER",
                   envir = environment()))
  MAT2 <- get(data("HEC_MI2",
                   package = "DECIPHER",
                   envir = environment()))
  structureMatrix <- matrix(c(0.187, -0.8, -0.873,
                              -0.8, 0.561, -0.979,
                              -0.873, -0.979, 0.221),
                            3,
                            3,
                            dimnames=list(c("H", "E", "C"),
                                          c("H", "E", "C")))
  substitutionMatrix <- matrix(c(1.5, -2.134, -0.739, -1.298,
                                 -2.134, 1.832, -2.462, 0.2,
                                 -0.739, -2.462, 1.522, -2.062,
                                 -1.298, 0.2, -2.062, 1.275),
                               nrow = 4,
                               dimnames = list(DNA_BASES, DNA_BASES))
  Features01 <- Features02 <- AAStruct <- FeatureLengths <- FeatureMods <- FeatureCode <- FeatureCDSCount <- vector("list",
                                                                                                                    length = length(GeneCalls))
  L <- length(GeneCalls)
  
  Count <- 1L
  while (object.size(Features01) < Storage &
         Count <= L) {
    # print(object.size(Features01))
    # print(Count)
    Genome <- SearchDB(dbFile = dbConn,
                       identifier = names(GeneCalls[Count]),
                       nameBy = "description",
                       type = "DNAStringSet",
                       verbose = FALSE)
    PresentIndices <- unique(GeneCalls[[Count]]$Index)
    # Reset any coding & non-translation table features to the default
    # move this somewhere else eventually...
    if (any(is.na(GeneCalls[[Count]]$Translation_Table))) {
      w <- which(is.na(GeneCalls[[Count]]$Translation_Table) &
                   GeneCalls[[Count]]$Coding)
      if (length(w) > 0) {
        GeneCalls[[Count]]$Translation_Table[w] <- DefaultTranslationTable
      }
    }
    if (length(PresentIndices) > 1L) {
      # many indices, loop through present indices and extract
      # slam together at the end
      Features01[[Count]] <- vector(mode = "list",
                                    length = length(PresentIndices))
      for (m3 in seq_along(PresentIndices)) {
        ph <- GeneCalls[[Count]]$Index == PresentIndices[m3]
        # implementation 3 - faster so far
        # set up the succinct extraction
        # build an index of where stringset positions need to be collapsed
        z1 <- unname(GeneCalls[[Count]]$Range[ph])
        z2 <- lengths(z1)
        # convert IRangesList to IRanges object for simple extractAt
        z1 <- unlist(z1,
                     recursive = FALSE)
        Features01[[Count]][[m3]] <- extractAt(x = Genome[[PresentIndices[m3]]],
                                               at = z1)
        CollapseCount <- 0L
        w <- which(z2 > 1L)
        # if no collapsing needs to occur, do not enter loop
        if (length(w) > 0L) {
          # if collapsing must take place build a placeholder of positions to remove
          # once collapsing correct positions has occurred
          remove <- vector(mode = "integer",
                           length = sum(z2[w]) - length(w))
          for (m4 in w) {
            Features01[[Count]][[m3]][[m4 + CollapseCount]] <- unlist(Features01[[Count]][[m3]][m4:(m4 + z2[m4] - 1L) + CollapseCount])
            remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
            CollapseCount <- CollapseCount + z2[m4] - 1L
          }
          Features01[[Count]][[m3]][remove] <- NULL
        }
        FlipMe <- GeneCalls[[Count]]$Strand[ph] == 1L
        if (any(FlipMe)) {
          Features01[[Count]][[m3]][FlipMe] <- reverseComplement(Features01[[Count]][[m3]][FlipMe])
        }
      }
      Features01[[Count]] <- do.call(c,
                                     Features01[[Count]])
      
    } else {
      # implementation 3 - shortest possible collapse loops and fewest copies - so far
      z1 <- unname(GeneCalls[[Count]]$Range)
      z2 <- lengths(z1)
      z1 <- unlist(z1,
                   recursive = FALSE)
      Features01[[Count]] <- extractAt(x = Genome[[PresentIndices]],
                                       at = z1)
      CollapseCount <- 0L
      w <- which(z2 > 1)
      if (length(w) > 0) {
        remove <- vector(mode = "integer",
                         length = sum(z2[w]) - length(w))
        for (m4 in w) {
          Features01[[Count]][[m4 + CollapseCount]] <- unlist(Features01[[Count]][m4:(m4 + z2[m4] - 1L) + CollapseCount])
          remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
          CollapseCount <- CollapseCount + z2[m4] - 1L
        }
        Features01[[Count]][remove] <- NULL
      }
      
      FlipMe <- GeneCalls[[Count]]$Strand == 1L
      if (any(FlipMe)) {
        Features01[[Count]][FlipMe] <- reverseComplement(Features01[[Count]][FlipMe])
      }
      
    }
    names(Features01[[Count]]) <- paste(rep(names(GeneCalls)[Count], length(Features01[[Count]])),
                                        GeneCalls[[Count]]$Index,
                                        seq(length(Features01[[Count]])),
                                        sep = "_")
    FeatureLengths[[Count]] <- width(Features01[[Count]])
    FeatureMods[[Count]] <- (FeatureLengths[[Count]] %% 3L) == 0L
    FeatureCode[[Count]] <- GeneCalls[[Count]]$Coding
    FeatureCDSCount[[Count]] <- lengths(GeneCalls[[Count]]$Range)
    
    # translate all translatable features with as few calls as possible
    ph <- unique(GeneCalls[[Count]]$Translation_Table)
    
    ph <- ph[!is.na(ph)]
    if (length(ph) < 1L) {
      ph <- DefaultTranslationTable
      phkey <- which(GeneCalls[[Count]]$Coding &
                       (GeneCalls[[Count]]$Type == "gene" | GeneCalls[[Count]]$Type == "pseudogene"))
      CurrentGeneticCode <- getGeneticCode(id_or_name2 = ph,
                                           full.search = FALSE,
                                           as.data.frame = FALSE)
      Features02[[Count]] <- suppressWarnings(translate(x = Features01[[Count]][phkey],
                                                        genetic.code = CurrentGeneticCode,
                                                        if.fuzzy.codon = "solve"))
      Features02[[Count]] <- Features02[[Count]][order(phkey)]
      # print(length(Features02[[Count]]))
    } else {
      Features02[[Count]] <- vector(mode = "list",
                                    length = length(ph))
      phkey <- vector(mode = "list",
                      length = length(ph))
      
      for (m4 in seq_along(ph)) {
        matchph <- which(GeneCalls[[Count]]$Translation_Table == ph[m4] &
                           GeneCalls[[Count]]$Coding &
                           (GeneCalls[[Count]]$Type == "gene" | GeneCalls[[Count]]$Type == "pseudogene"))
        phkey[[m4]] <- matchph
        CurrentGeneticCode <- getGeneticCode(id_or_name2 = ph[m4],
                                             full.search = FALSE,
                                             as.data.frame = FALSE)
        Features02[[Count]][[m4]] <- suppressWarnings(translate(x = Features01[[Count]][matchph],
                                                                genetic.code = CurrentGeneticCode,
                                                                if.fuzzy.codon = "solve"))
      }
      Features02[[Count]] <- do.call(c,
                                     Features02[[Count]])
      phkey <- unlist(phkey)
      Features02[[Count]] <- Features02[[Count]][order(phkey)]
      
    }
    # rewrite ph to provide the correct names for the features
    ph <- GeneCalls[[Count]]$Coding & (GeneCalls[[Count]]$Type == "gene" | GeneCalls[[Count]]$Type == "pseudogene")
    
    names(Features02[[Count]]) <- paste(rep(names(GeneCalls)[Count], length(Features02[[Count]])),
                                        GeneCalls[[Count]]$Index[ph],
                                        seq(length(Features01[[Count]]))[ph],
                                        sep = "_")
    
    # generate structures for aa alignments if there are 
    if (!is.null(Features02[[Count]])) {
      AAStruct[[Count]] <- PredictHEC(myAAStringSet = Features02[[Count]],
                                      type = "probabilities",
                                      HEC_MI1 = MAT1,
                                      HEC_MI2 = MAT2)
      names(AAStruct[[Count]]) <- names(Features02[[Count]])
    }
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = Count / L)
    }
    
    Count <- Count + 1L
    # will extract till storage is exceeded
    # will not cap at storage
  } # end while loop
  
  if (Verbose) {
    close(pBar)
    if (Count < L) {
      cat("Overhead is larger than storage request.\nSearch loops will require database lookups.\n")
      # RemoveWhenAble <- TRUE
    } else {
      cat("Complete!\n")
      # RemoveWhenAble <- FALSE
    }
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  
  w1 <- sapply(X = Features01,
               FUN = function(x) {
                 !is.null(x)
               },
               simplify = TRUE)
  
  res <- list("DNA" = Features01[w1],
              "AA" = Features02[w1],
              "Struct" = AAStruct[w1],
              "IDs" = names(GeneCalls)[w1],
              "NTCount" = FeatureLengths[w1],
              "CodingVal1" = FeatureMods[w1],
              "CodingVal2" = FeatureCode[w1],
              "CDSCount" = FeatureCDSCount[w1])
  class(res) <- c("list", "FeatureSeqs")
  return(res)
}



