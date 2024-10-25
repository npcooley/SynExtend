###### -- extract a bunch of sequences from some genecalls --------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu

# TODO
# instead of storing these as an object in the R session, send them to a DB

PrepareSeqs <- function(SynExtendObject,
                        DataBase01,
                        DefaultTranslationTable = "11",
                        Identifiers = NULL,
                        Verbose = FALSE) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1)
  }
  
  # check object type -- currently only PairSummaries and LinkedPairs
  if (!is(object = SynExtendObject,
          class2 = "LinkedPairs")) {
    stop("'SynExtendObject' must be an object of class 'LinkedPairs'.")
  }
  # initialize database01
  if (is.character(DataBase01)) {
    if (!requireNamespace(package = "RSQLite",
                          quietly = TRUE)) {
      stop("Package 'RSQLite' must be installed.")
    }
    if (!("package:RSQLite" %in% search())) {
      print("Eventually character vector access to DECIPHER DBs will be deprecated.")
      requireNamespace(package = "RSQLite",
                       quietly = TRUE)
    }
    dbConn01 <- dbConnect(dbDriver("SQLite"), DataBase01)
    on.exit(dbDisconnect(dbConn01))
  } else {
    dbConn01 <- DataBase01
    if (!dbIsValid(dbConn01)) {
      stop("The connection has expired.")
    }
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
  } else {
    Identifiers <- names(GeneCalls)
  }
  if (length(GeneCalls) < 1) {
    stop("Requested subset not found.")
  }
  
  # load in genomes and ALL extracted features at the top until storage limit
  # is hit
  if (Verbose) {
    cat("Preparing overhead data.\n")
  }
  # it's too expensive to store anything but the seqs in the db
  # load in structure matrices once for PredictHEC
  # these hardcoded values come from Erik ... AlignProfiles maybe?
  # MAT1 <- get(data("HEC_MI1",
  #                  package = "DECIPHER",
  #                  envir = environment()))
  # MAT2 <- get(data("HEC_MI2",
  #                  package = "DECIPHER",
  #                  envir = environment()))
  # structureMatrix <- matrix(c(0.187, -0.8, -0.873,
  #                             -0.8, 0.561, -0.979,
  #                             -0.873, -0.979, 0.221),
  #                           3,
  #                           3,
  #                           dimnames=list(c("H", "E", "C"),
  #                                         c("H", "E", "C")))
  # substitutionMatrix <- matrix(c(1.5, -2.134, -0.739, -1.298,
  #                                -2.134, 1.832, -2.462, 0.2,
  #                                -0.739, -2.462, 1.522, -2.062,
  #                                -1.298, 0.2, -2.062, 1.275),
  #                              nrow = 4,
  #                              dimnames = list(DNA_BASES,
  #                                              DNA_BASES))
  # Features01 <- Features02 <- AAStruct <- FeatureLengths <- FeatureMods <- FeatureCode <- FeatureCDSCount <- vector("list",
  #                                                                                                                   length = length(GeneCalls))
  
  # test <- vector(mode = "list",
  #                length = length(GeneCalls))
  L <- length(GeneCalls)
  Count <- 1L
  while (Count <= L) {
    # print(object.size(Features01))
    # print(Count)
    Genome <- SearchDB(dbFile = dbConn01,
                       tblName = "Seqs",
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
      Features01 <- vector(mode = "list",
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
        Features01[[m3]] <- extractAt(x = Genome[[PresentIndices[m3]]],
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
            Features01[[m3]][[m4 + CollapseCount]] <- unlist(Features01[[m3]][m4:(m4 + z2[m4] - 1L) + CollapseCount])
            remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
            CollapseCount <- CollapseCount + z2[m4] - 1L
          }
          Features01[[m3]][remove] <- NULL
        }
        FlipMe <- GeneCalls[[Count]]$Strand[ph] == 1L
        if (any(FlipMe)) {
          Features01[[m3]][FlipMe] <- reverseComplement(Features01[[m3]][FlipMe])
        }
      }
      Features01 <- do.call(c,
                            Features01)
      
    } else {
      # implementation 3 - shortest possible collapse loops and fewest copies - so far
      z1 <- unname(GeneCalls[[Count]]$Range)
      z2 <- lengths(z1)
      z1 <- unlist(z1,
                   recursive = FALSE)
      Features01 <- extractAt(x = Genome[[PresentIndices]],
                                       at = z1)
      CollapseCount <- 0L
      w <- which(z2 > 1)
      if (length(w) > 0) {
        remove <- vector(mode = "integer",
                         length = sum(z2[w]) - length(w))
        for (m4 in w) {
          Features01[[m4 + CollapseCount]] <- unlist(Features01[m4:(m4 + z2[m4] - 1L) + CollapseCount])
          remove[(CollapseCount + 1L):(CollapseCount + z2[m4] - 1L)] <- (m4 + 1L):(m4 + z2[m4] - 1L) + CollapseCount
          CollapseCount <- CollapseCount + z2[m4] - 1L
        }
        Features01[remove] <- NULL
      }
      
      FlipMe <- GeneCalls[[Count]]$Strand == 1L
      if (any(FlipMe)) {
        Features01[FlipMe] <- reverseComplement(Features01[FlipMe])
      }
      
    }
    names(Features01) <- paste(rep(names(GeneCalls)[Count], length(Features01)),
                                        GeneCalls[[Count]]$Index,
                                        seq(length(Features01)),
                                        sep = "_")
    # FeatureLengths <- width(Features01)
    # FeatureMods <- (FeatureLengths %% 3L) == 0L
    # FeatureCode <- GeneCalls[[Count]]$Coding
    # FeatureCDSCount <- lengths(GeneCalls[[Count]]$Range)
    
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
      Features02 <- suppressWarnings(translate(x = Features01[phkey],
                                                        genetic.code = CurrentGeneticCode,
                                                        if.fuzzy.codon = "solve"))
      Features02 <- Features02[order(phkey)]
      # print(length(Features02))
    } else {
      Features02 <- vector(mode = "list",
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
        Features02[[m4]] <- suppressWarnings(translate(x = Features01[matchph],
                                                                genetic.code = CurrentGeneticCode,
                                                                if.fuzzy.codon = "solve"))
      }
      Features02 <- do.call(c,
                                     Features02)
      phkey <- unlist(phkey)
      Features02 <- Features02[order(phkey)]
      
    }
    # rewrite ph to provide the correct names for the features
    ph <- GeneCalls[[Count]]$Coding & (GeneCalls[[Count]]$Type == "gene" | GeneCalls[[Count]]$Type == "pseudogene")
    
    names(Features02) <- paste(rep(names(GeneCalls)[Count], length(Features02)),
                                        GeneCalls[[Count]]$Index[ph],
                                        seq(length(Features01))[ph],
                                        sep = "_")
    
    Seqs2DB(seqs = Features01,
            dbFile = dbConn01,
            tblName = "NTs",
            type = "XStringSet",
            verbose = FALSE,
            identifier = Identifiers[Count])
    Seqs2DB(seqs = Features02,
            dbFile = dbConn01,
            tblName = "AAs",
            type = "XStringSet",
            verbose = FALSE,
            identifier = Identifiers[Count])
    # new_rows <- dbGetQuery(conn = dbConn01,
    #                        statement = paste("select row_names from NTs where identifier is",
    #                                          Identifiers[Count]))$row_names
    # # print(length(new_rows))
    # # print(head(new_rows))
    # dat01 <- data.frame("len" = FeatureLengths,
    #                     "mod" = as.integer(FeatureMods),
    #                     "code" = as.integer(FeatureCode),
    #                     "cds" = FeatureCDSCount,
    #                     row.names = new_rows)
    # Add2DB(myData = dat01,
    #        dbFile = dbConn01,
    #        tblName = "NTs",
    #        verbose = FALSE)
    # print(Count)
    # print(Identifiers[Count])
    # print(head(dat01))
    # print(dim(dat01))
    # print(length(Features01))
    # print(paste("identifier is",
    #             Identifiers[Count]))
    # test[[Count]] <- list("aa" = Features01,
    #                       "nt" = Features02,
    #                       "dat" = dat01)
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = Count / L)
    }
    
    Count <- Count + 1L
  } # end while loop
  
  if (Verbose) {
    close(pBar)
    cat("Complete!\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  # return(test)
  invisible(Count)
}



