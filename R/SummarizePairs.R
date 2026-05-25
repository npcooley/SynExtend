###### -- summarize pairs from a LinkedPairs object ---------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: Nicholas.Cooley@ul.ie
# given a linked pairs object, and a database connection return a pairsummaries
# object

###### -- NOTES ---------------------------------------------------------------
# this function is getting rewritten yet again
# the goal is to slim down this function to not be such a swiss army knife
# and make workflows involving this function a bit more intuitive
# this is just going to take the LinkedPairs object and slam it into a DF with
# alignment, kmer, match locality, and block context values and stats

# under the hood DB stuff is getting dealt with as well ... the PrepareSeqs function
# will eventually be changed, but the NT table is going to be dropped,
# extraction and munging the genomic stringset is pretty cheap compute wise
# so we will always generate the nucs on the fly, but we WILL be creating an
# AA table if one is not present

###### -- FUNCTION ------------------------------------------------------------

SummarizePairs <- function(SynExtendObject,
                           DataBase01,
                           DefaultTranslationTable = "11",
                           KmerSize = 5,
                           Verbose = FALSE,
                           ShowPlot = FALSE,
                           Processors = 1,
                           Storage = 2,
                           Anchors = c("enforce", "infer", "ignore"),
                           ...) {
  # Verbosity check
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  # overhead checking
  # object types
  if (!is(object = SynExtendObject,
          class2 = "LinkedPairs")) {
    stop ("'SynExtendObject' is not an object of class 'LinkedPairs'.")
  }
  Size <- nrow(SynExtendObject)
  
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
    dbConn <- dbConnect(dbDriver("SQLite"), DataBase01)
    on.exit(dbDisconnect(dbConn))
  } else {
    dbConn <- DataBase01
    if (!dbIsValid(dbConn)) {
      stop("The connection has expired.")
    }
  }
  if (!is.character(DefaultTranslationTable) |
      length(DefaultTranslationTable) > 1) {
    stop("DefaultTranslationTable must be a character of length 1.")
  }
  # check storage
  if (Storage < 0) {
    stop("Storage must be greater than zero.")
  } else {
    Storage <- Storage * 1e9 # conversion to gigabytes
  }
  # deal with Processors, this mimics Erik's error checking
  if (!is.null(Processors) && !is.numeric(Processors)) {
    stop("Processors must be a numeric.")
  }
  if (!is.null(Processors) && floor(Processors) != Processors) {
    stop("Processors must be a whole number.")
  }
  if (!is.null(Processors) && Processors < 1) {
    stop("Processors must be at least 1.")
  }
  if (is.null(Processors)) {
    Processors <- .Call("detectCores",
                        PACKAGE="DECIPHER")
  } else {
    Processors <- as.integer(Processors)
  }
  Anchors <- match.arg(Anchors)
  # deal with user arguments
  # ignore 'verbose'
  UserArgs <- list(...)
  # if (length(UserArgs) > 0) {
  #   UserArgNames <- names(UserArgs)
  #   AlignPairsArgs <- formals(AlignPairs)
  #   AlignPairsArgNames <- names(AlignPairsArgs)
  #   AlignProfilesArgs <- formals(AlignProfiles)
  #   AlignProfilesArgNames <- names(AlignProfilesArgs)
  #   DistanceMatrixArgs <- formals(DistanceMatrix)
  #   DistanceMatrixArgNames <- names(DistanceMatrixArgs)
  #   
  #   APaArgNames <- APrArgNames <- DMArgNames <- vector(mode = "character",
  #                                                    length = length(UserArgNames))
  #   for (a1 in seq_along(UserArgNames)) {
  #     APaArgNames[a1] <- match.arg(arg = UserArgNames[a1],
  #                                  choices = AlignPairsArgNames)
  #     APrArgNames[a1] <- match.arg(arg = UserArgNames[a1],
  #                                  choices = AlignProfilesArgNames)
  #     DMArgNames[a1] <- match.arg(arg = UserArgNames[a1],
  #                                 choices = DistanceMatrixArgNames)
  #   }
  #   # set the user args for AlignProfiles
  #   if (any(APaArgNames)) {
  #     
  #   }
  #   # set the user args for AlignPairs
  #   # set the user args for DistanceMatrix
  # }
  
  GeneCalls <- attr(x = SynExtendObject,
                    which = "GeneCalls")
  GeneCallIDs <- names(GeneCalls)
  ObjectIDs <- rownames(SynExtendObject)
  # when we subset a LinkedPairsObject it doesn't smartly handle the genecalls yet...
  if (!all(ObjectIDs %in% GeneCallIDs)) {
    stop("Function expects all IDs in the SynExtendObject to have supplied GeneCalls.")
  }
  
  # assign some things
  AA_matrix <- DECIPHER:::.getSubMatrix("PFASUM50")
  NT_matrix <- DECIPHER:::.nucleotideSubstitutionMatrix(2L, -1L, 1L)
  
  feature_match <- match(x = ObjectIDs,
                         table = GeneCallIDs)
  # These need to all be switched over to feature_match
  # we're only going to scroll through the cells that are supplied in the object
  # the datapool only needs to be as long as the gene calls object
  DataPool <- vector(mode = "list",
                     length = length(GeneCallIDs))
  
  # AlignSeqs has been dropped as an option, so these is no long necessary
  # MAT1 <- get(data("HEC_MI1",
  #                  package = "DECIPHER",
  #                  envir = environment()))
  # MAT2 <- get(data("HEC_MI2",
  #                  package = "DECIPHER",
  #                  envir = environment()))
  
  # progress bar ticks through each cell
  # technically these alignments have the possibility of not being
  # as good as align profiles, but that's a hit we're willing to take
  Total <- (Size * (Size - 1L)) / 2L
  PH <- Attr_PH <- vector(mode = "list",
                          length = Total)
  Total <- Total * 2L
  if (Verbose) {
    cat("Collecting pairs.\n")
  }
  
  Count <- 1L
  PBCount <- 0L
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
  # 1 == querygene
  # 2 == subjectgene
  # 3 == exactoverlap
  # 4 == queryindex
  # 5 == subjectindex
  # 6 == qleftpos
  # 7 == qrightpos
  # 8 == sleftpos
  # 9 == srightpos
  # 10 == querysubkey
  # 11 == subjectsubkey
  # 12 == syntenicorigin
  block_uid <- 0L
  Prev_a1 <- 0L
  
  # db query stuff
  aa_table_present <- dbExistsTable(conn = dbConn,
                                    name = "AAs")
  if (aa_table_present) {
    ids_in_db <- dbGetQuery(conn = dbConn,
                            statement = "select distinct identifier from AAs")$identifier
    db_id_present <- ObjectIDs %in% ids_in_db
    names(db_id_present) <- ObjectIDs
  } else {
    db_id_present <- rep(FALSE,
                         length(ObjectIDs))
    names(db_id_present) <- ObjectIDs
  }
  
  
  for (a1 in seq_len(Size - 1L)) {
    for (a2 in (a1 + 1L):Size) {
      # print(c(a1, a2))
      # build our data pool first
      # regardless of how we align or if we're including search index or not, we need to prepare and collect
      # the same basic statistics and look aheads
      if (is.null(DataPool[[feature_match[a1]]])) {
        # the pool position is empty, pull from the DB
        # and generate the AAStructures
        # there's no real reason to carry the NT features around,
        # so we regenerate them every time
        # we store the AAs if we make them
        # we grab them if they're already present
        genome_ph <- SearchDB(dbFile = dbConn,
                              identifier = ObjectIDs[a1],
                              verbose = FALSE,
                              nameBy = "description",
                              type = "DNAStringSet")
        seqs_ph <- FeaturesFromDF(Genome = genome_ph,
                                  GeneCalls = GeneCalls[[ObjectIDs[a1]]],
                                  Index = ObjectIDs[a1])
        # names(seqs_ph) <- paste(rep(ObjectIDs[a1],
        #                             length(seqs_ph)),
        #                         GeneCalls[[ObjectIDs[a1]]]$Index,
        #                         seq(length(seqs_ph)),
        #                         sep = "_")
        DataPool[[feature_match[a1]]]$DNA <- seqs_ph
        
        if (db_id_present[ObjectIDs[a1]]) {
          # just pull the AAs, they already exist
          DataPool[[feature_match[a1]]]$AA <- SearchDB(dbFile = dbConn,
                                                       tblName = "AAs",
                                                       identifier = ObjectIDs[a1],
                                                       verbose = FALSE,
                                                       nameBy = "description",
                                                       type = "AAStringSet")
        } else {
          # translate where we can, seqs don't exist currently
          tr_tbl1 <- GeneCalls[[ObjectIDs[a1]]]$Coding
          tr_tbl2 <- GeneCalls[[ObjectIDs[a1]]]$Translation_Table
          tr_tbl3 <- lapply(X = GeneCalls[[ObjectIDs[a1]]]$Range,
                            FUN = function(x) {
                              width(x)
                            })
          tr_tbl3 <- vapply(X = tr_tbl3,
                            FUN = function(x) {
                              sum(x)
                            },
                            FUN.VALUE = vector(mode = "integer",
                                               length = 1L))
          tr_tbl3 <- (tr_tbl3 %% 3) == 0
          
          if (any(is.na(tr_tbl2))) {
            w <- which(is.na(tr_tbl2) &
                         tr_tbl1 &
                         tr_tbl3)
            if (length(w) > 0) {
              tr_tbl2[w] <- DefaultTranslationTable
            }
          }
          
          tr_tbl4 <- unique(tr_tbl2[!is.na(tr_tbl2)])
          aa_ph <- vector(mode = "list",
                          length = length(tr_tbl4))
          w2 <- tr_tbl1 & tr_tbl3
          for (a3 in seq_along(tr_tbl4)) {
            current_genetic_code <- getGeneticCode(id_or_name2 = tr_tbl4[a3],
                                                   full.search = FALSE,
                                                   as.data.frame = FALSE)
            w1 <- tr_tbl2 == tr_tbl4[a3]
            aa_ph[[a3]] <- translate(x = seqs_ph[w1 & w2],
                                     genetic.code = current_genetic_code,
                                     if.fuzzy.codon = "solve")
          }
          if (a3 > 1) {
            # slam the list together and reinforce the order
            aa_seqs <- do.call(c,
                               aa_ph)
            # my ridiculous renaming scheme comes in handy here, if nowhere else
            # return(aa_seqs)
            o1 <- vapply(X = strsplit(x = names(aa_seqs),
                                      split = "_",
                                      fixed = TRUE),
                         FUN = function(x) {
                           x[3]
                         },
                         FUN.VALUE = vector(mode = "character",
                                            length = 1))
            o1 <- order(as.integer(o1))
            aa_seqs <- aa_seqs[o1]
          } else {
            # order from seqs_ph should be retained
            aa_seqs <- aa_ph[[1]]
          }
          DataPool[[feature_match[a1]]]$AA <- aa_seqs
          Seqs2DB(seqs = aa_seqs,
                  dbFile = dbConn,
                  tblName = "AAs",
                  type = "XStringSet",
                  verbose = FALSE,
                  identifier = ObjectIDs[a1])
          # reset this if we have to come back here in a very large object
          db_id_present[ObjectIDs[a1]] <- TRUE
        } # end db ids check
        
        DataPool[[feature_match[a1]]]$len <- width(DataPool[[feature_match[a1]]]$DNA)
        DataPool[[feature_match[a1]]]$mod <- DataPool[[feature_match[a1]]]$len %% 3L == 0
        DataPool[[feature_match[a1]]]$code <- GeneCalls[[feature_match[a1]]]$Coding
        DataPool[[feature_match[a1]]]$cds <- lengths(GeneCalls[[feature_match[a1]]]$Range)
        
        # DataPool[[feature_match[a1]]]$struct <- PredictHEC(myAAStringSet = DataPool[[feature_match[a1]]]$AA,
        #                                                    type = "probabilities",
        #                                                    HEC_MI1 = MAT1,
        #                                                    HEC_MI2 = MAT2)
        DataPool[[feature_match[a1]]]$aa_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[a1]]]$AA)
        DataPool[[feature_match[a1]]]$aa_backgrounds <- DataPool[[feature_match[a1]]]$aa_backgrounds[, colnames(AA_matrix)]
        DataPool[[feature_match[a1]]]$aa_backgrounds <- DataPool[[feature_match[a1]]]$aa_backgrounds / rowSums(DataPool[[feature_match[a1]]]$aa_backgrounds)
        DataPool[[feature_match[a1]]]$dna_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[a1]]]$DNA)
        DataPool[[feature_match[a1]]]$dna_backgrounds <- DataPool[[feature_match[a1]]]$dna_backgrounds[, colnames(NT_matrix)]
        DataPool[[feature_match[a1]]]$dna_backgrounds <- DataPool[[feature_match[a1]]]$dna_backgrounds / rowSums(DataPool[[feature_match[a1]]]$dna_backgrounds)
        DataPool[[feature_match[a1]]]$aa_register <- match(table = which(DataPool[[feature_match[a1]]]$mod &
                                                                           DataPool[[feature_match[a1]]]$code),
                                                           x = seq(length(DataPool[[feature_match[a1]]]$DNA)))
        DataPool[[feature_match[a1]]]$phase <- GeneCalls[[feature_match[a1]]]$Phase
        DataPool[[feature_match[a1]]]$ranges <- GeneCalls[[feature_match[a1]]]$Range
        DataPool[[feature_match[a1]]]$strand <- GeneCalls[[feature_match[a1]]]$Strand
        # DataPool[[feature_match[a1]]]$index <- do.call(what = "IndexSeqs",
        #                                                args = c(list("subject" = DataPool[[feature_match[a1]]]$AA,
        #                                                              "verbose" = FALSE),
        #                                                         IndexParams))
        
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
      } # end a1 data pool check
      
      if (is.null(DataPool[[feature_match[a2]]])) {
        # the pool position is empty
        # collect the AAs if they exist, generate them if they don't
        # generate the NTs, because storage of them is superfluous
        # 
        # the pool position is empty, pull from the DB
        # and generate the AAStructures
        # there's no real reason to carry the NT features around,
        # so we regenerate them every time
        # we store the AAs if we make them
        # we grab them if they're already present
        genome_ph <- SearchDB(dbFile = dbConn,
                              identifier = ObjectIDs[a2],
                              verbose = FALSE,
                              nameBy = "description",
                              type = "DNAStringSet")
        seqs_ph <- FeaturesFromDF(Genome = genome_ph,
                                  GeneCalls = GeneCalls[[ObjectIDs[a2]]],
                                  Index = ObjectIDs[a2])
        # names(seqs_ph) <- paste(rep(ObjectIDs[a2],
        #                             length(seqs_ph)),
        #                         GeneCalls[[ObjectIDs[a2]]]$Index,
        #                         seq(length(seqs_ph)),
        #                         sep = "_")
        DataPool[[feature_match[a2]]]$DNA <- seqs_ph
        
        if (db_id_present[ObjectIDs[a2]]) {
          # just pull the AAs, they already exist
          DataPool[[feature_match[a2]]]$AA <- SearchDB(dbFile = dbConn,
                                                       tblName = "AAs",
                                                       identifier = ObjectIDs[a2],
                                                       verbose = FALSE,
                                                       nameBy = "description",
                                                       type = "AAStringSet")
        } else {
          # translate where we can, seqs don't exist currently
          tr_tbl1 <- GeneCalls[[ObjectIDs[a2]]]$Coding
          tr_tbl2 <- GeneCalls[[ObjectIDs[a2]]]$Translation_Table
          tr_tbl3 <- lapply(X = GeneCalls[[ObjectIDs[a2]]]$Range,
                            FUN = function(x) {
                              width(x)
                            })
          tr_tbl3 <- vapply(X = tr_tbl3,
                            FUN = function(x) {
                              sum(x)
                            },
                            FUN.VALUE = vector(mode = "integer",
                                               length = 1L))
          tr_tbl3 <- (tr_tbl3 %% 3) == 0
          
          if (any(is.na(tr_tbl2))) {
            w <- which(is.na(tr_tbl2) &
                         tr_tbl1)
            if (length(w) > 0) {
              tr_tbl2[w] <- DefaultTranslationTable
            }
          }
          
          tr_tbl4 <- unique(tr_tbl2[!is.na(tr_tbl2)])
          aa_ph <- vector(mode = "list",
                          length = length(tr_tbl4))
          w2 <- tr_tbl1 & tr_tbl3
          for (a3 in seq_along(tr_tbl4)) {
            current_genetic_code <- getGeneticCode(id_or_name2 = tr_tbl4[a3],
                                                   full.search = FALSE,
                                                   as.data.frame = FALSE)
            w1 <- tr_tbl2 == tr_tbl4[a3]
            aa_ph[[a3]] <- translate(x = seqs_ph[w1 & w2],
                                     genetic.code = current_genetic_code,
                                     if.fuzzy.codon = "solve")
          }
          if (a3 > 1) {
            # slam the list together and reinforce the order
            aa_seqs <- do.call(c,
                               aa_ph)
            # my ridiculous renaming scheme comes in handy here, if nowhere else
            o1 <- vapply(X = strsplit(x = names(aa_seqs),
                                      split = "_",
                                      fixed = TRUE),
                         FUN = function(x) {
                           x[3]
                         },
                         FUN.VALUE = vector(mode = "character",
                                            length = 1))
            o1 <- order(as.integer(o1))
            aa_seqs <- aa_seqs[o1]
          } else {
            # order from seqs_ph should be retained
            aa_seqs <- aa_ph[[1]]
          }
          DataPool[[feature_match[a2]]]$AA <- aa_seqs
          Seqs2DB(seqs = aa_seqs,
                  dbFile = dbConn,
                  tblName = "AAs",
                  type = "XStringSet",
                  verbose = FALSE,
                  identifier = ObjectIDs[a2])
          # reset this if we have to come back here in a very large object
          db_id_present[ObjectIDs[a2]] <- TRUE
        } # end db ids check
        
        DataPool[[feature_match[a2]]]$len <- width(DataPool[[feature_match[a2]]]$DNA)
        DataPool[[feature_match[a2]]]$mod <- DataPool[[feature_match[a2]]]$len %% 3L == 0
        DataPool[[feature_match[a2]]]$code <- GeneCalls[[feature_match[a2]]]$Coding
        DataPool[[feature_match[a2]]]$cds <- lengths(GeneCalls[[feature_match[a2]]]$Range)
        
        # DataPool[[feature_match[a2]]]$struct <- PredictHEC(myAAStringSet = DataPool[[feature_match[a2]]]$AA,
        #                                                    type = "probabilities",
        #                                                    HEC_MI1 = MAT1,
        #                                                    HEC_MI2 = MAT2)
        DataPool[[feature_match[a2]]]$aa_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[a2]]]$AA)
        DataPool[[feature_match[a2]]]$aa_backgrounds <- DataPool[[feature_match[a2]]]$aa_backgrounds[, colnames(AA_matrix)]
        DataPool[[feature_match[a2]]]$aa_backgrounds <- DataPool[[feature_match[a2]]]$aa_backgrounds / rowSums(DataPool[[feature_match[a2]]]$aa_backgrounds)
        DataPool[[feature_match[a2]]]$dna_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[a2]]]$DNA)
        DataPool[[feature_match[a2]]]$dna_backgrounds <- DataPool[[feature_match[a2]]]$dna_backgrounds[, colnames(NT_matrix)]
        DataPool[[feature_match[a2]]]$dna_backgrounds <- DataPool[[feature_match[a2]]]$dna_backgrounds / rowSums(DataPool[[feature_match[a2]]]$dna_backgrounds)
        DataPool[[feature_match[a2]]]$aa_register <- match(table = which(DataPool[[feature_match[a2]]]$mod &
                                                                           DataPool[[feature_match[a2]]]$code),
                                                           x = seq(length(DataPool[[feature_match[a2]]]$DNA)))
        DataPool[[feature_match[a2]]]$phase <- GeneCalls[[feature_match[a2]]]$Phase
        DataPool[[feature_match[a2]]]$ranges <- GeneCalls[[feature_match[a2]]]$Range
        DataPool[[feature_match[a2]]]$strand <- GeneCalls[[feature_match[a2]]]$Strand
        # DataPool[[feature_match[a2]]]$index <- do.call(what = "IndexSeqs",
        #                                                args = c(list("subject" = DataPool[[feature_match[a2]]]$AA,
        #                                                              "verbose" = FALSE),
        #                                                         IndexParams))
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
      } # end a2 data pool check
      
      
      if (Prev_a1 != a1) {
        QueryDNA <- DataPool[[feature_match[a1]]]$DNA
        QueryAA <- DataPool[[feature_match[a1]]]$AA
        QNTCount <- DataPool[[feature_match[a1]]]$len
        QMod <- DataPool[[feature_match[a1]]]$mod
        QCode <- DataPool[[feature_match[a1]]]$code
        QCDSCount <- DataPool[[feature_match[a1]]]$cds
        # QueryStruct <- DataPool[[feature_match[a1]]]$struct
        # QueryIndex <- DataPool[[feature_match[a1]]]$index
        Query_Background_AA <- DataPool[[feature_match[a1]]]$aa_backgrounds
        Query_Background_NT <- DataPool[[feature_match[a1]]]$dna_backgrounds
        Q_AA_Register <- DataPool[[feature_match[a1]]]$aa_register
        QPhase <- DataPool[[feature_match[a1]]]$phase
        QRange <- DataPool[[feature_match[a1]]]$ranges
        QStrand <- DataPool[[feature_match[a1]]]$strand
        if (KmerSize < 10L) {
          Q_NucF <- oligonucleotideFrequency(x = QueryDNA,
                                             width = KmerSize,
                                             as.prob = TRUE)
        } else {
          stop ("non-overlapping kmers not implemented")
        }
      } else {
        # do something else?
      }
      
      # a2 never stays the same, it will always change in the current search
      # strategy
      SubjectDNA <- DataPool[[feature_match[a2]]]$DNA
      SubjectAA <- DataPool[[feature_match[a2]]]$AA
      SNTCount <- DataPool[[feature_match[a2]]]$len
      SMod <- DataPool[[feature_match[a2]]]$mod
      SCode <- DataPool[[feature_match[a2]]]$code
      SCDSCount <- DataPool[[feature_match[a2]]]$cds
      # SubjectStruct <- DataPool[[feature_match[a2]]]$struct
      # SubjectIndex <- DataPool[[feature_match[a2]]]$index
      Subject_Background_AA <- DataPool[[feature_match[a2]]]$aa_backgrounds
      Subject_Background_NT <- DataPool[[feature_match[a2]]]$dna_backgrounds
      S_AA_Register <- DataPool[[feature_match[a2]]]$aa_register
      SPhase <- DataPool[[feature_match[a2]]]$phase
      SRange <- DataPool[[feature_match[a2]]]$ranges
      SStrand <- DataPool[[feature_match[a2]]]$strand
      if (KmerSize < 10L) {
        S_NucF <- oligonucleotideFrequency(x = SubjectDNA,
                                           width = KmerSize,
                                           as.prob = TRUE)
      } else {
        stop ("non-overlapping kmers not implemented")
      }
      # data pool management stuff is wrapped up here ...
      
      # no search strategy nonsense, just summarize the pairs that we happen to see here
      # this check should be fine because these get filled with empty matrices
      # during construction before being filled with values so they will not be NULL
      if (nrow(SynExtendObject[[a1, a2]]) > 0L) {
        # links table is populated, do whatever
        PMatrix <- cbind(SynExtendObject[[a1, a2]][, 1L],
                         SynExtendObject[[a1, a2]][, 2L])
        IMatrix <- cbind(SynExtendObject[[a1, a2]][, 4L],
                         SynExtendObject[[a1, a2]][, 5L])
        
        # find the Consensus of each linking kmer hit
        
        diff1 <- HitConsensus(gene1left = GeneCalls[[feature_match[a1]]]$Start[SynExtendObject[[a2, a1]][, 1L]],
                              gene2left = GeneCalls[[feature_match[a2]]]$Start[SynExtendObject[[a2, a1]][, 2L]],
                              gene1right = GeneCalls[[feature_match[a1]]]$Stop[SynExtendObject[[a2, a1]][, 1L]],
                              gene2right = GeneCalls[[feature_match[a2]]]$Stop[SynExtendObject[[a2, a1]][, 2L]],
                              hit1left = SynExtendObject[[a2, a1]][, 6L],
                              hit1right = SynExtendObject[[a2, a1]][, 7L],
                              hit2left = SynExtendObject[[a2, a1]][, 8L],
                              hit2right = SynExtendObject[[a2, a1]][, 9L],
                              strand1 = GeneCalls[[feature_match[a1]]]$Strand[SynExtendObject[[a2, a1]][, 1L]],
                              strand2 = GeneCalls[[feature_match[a2]]]$Strand[SynExtendObject[[a2, a1]][, 2L]])
        # get the mean consensus
        diff2 <- vector(mode = "numeric",
                        length = nrow(SynExtendObject[[a1, a2]]))
        hit_relations <- vector(mode = "integer",
                                length = nrow(SynExtendObject[[a2, a1]]))
        hit_iterator <- 0L
        loop_iterator <- 1L
        h1 <- 0L
        h2 <- 0L
        continue <- TRUE
        while (continue) {
          if (SynExtendObject[[a2, a1]][loop_iterator, 1L] == h1 &
              SynExtendObject[[a2, a1]][loop_iterator, 2L] == h2) {
            hit_relations[loop_iterator] <- hit_iterator
            loop_iterator <- loop_iterator + 1L
          } else {
            hit_iterator <- hit_iterator + 1L
            hit_relations[loop_iterator] <- hit_iterator
            h1 <- SynExtendObject[[a2, a1]][loop_iterator, 1L]
            h2 <- SynExtendObject[[a2, a1]][loop_iterator, 2L]
            loop_iterator <- loop_iterator + 1L
          }
          # setTxtProgressBar(pb = pBar,
          #                   value = loop_iterator / nrow(SynExtendObject[[a2, a1]]))
          if (loop_iterator > nrow(SynExtendObject[[a2, a1]])) {
            continue <- FALSE
          }
        }
        
        diff2 <- 1 - unname(tapply(X = diff1,
                                   INDEX = hit_relations,
                                   FUN = function(x) {
                                     mean(x)
                                   }))
        # max match size
        MatchMax <- SynExtendObject[[a1, a2]][, "MaxKmerSize"]
        # total unique matches
        UniqueMatches <- SynExtendObject[[a1, a2]][, "TotalKmerHits"]
        # total matches at all
        TotalMatch <- SynExtendObject[[a1, a2]][, "ExactOverlap"]
        
        diff3 <- ApproximateBackground(p1 = SynExtendObject[[a1, a2]][, 1L],
                                       p2 = SynExtendObject[[a1, a2]][, 2L],
                                       code1 = QCode[SynExtendObject[[a1, a2]][, 1L]],
                                       code2 = SCode[SynExtendObject[[a1, a2]][, 2L]],
                                       mod1 = QMod[SynExtendObject[[a1, a2]][, 1L]],
                                       mod2 = SMod[SynExtendObject[[a1, a2]][, 2L]],
                                       aa1 = Query_Background_AA,
                                       aa2 = Subject_Background_AA,
                                       nt1 = Query_Background_NT,
                                       nt2 = Subject_Background_NT,
                                       register1 = Q_AA_Register,
                                       register2 = S_AA_Register,
                                       aamat = AA_matrix,
                                       ntmat = NT_matrix)
        # print("e")
        # from here we need to get the kmer differences
        # the PIDs
        # the SCOREs
        
        # if both positions are present in the FeatureSeqs object, do nothing
        # if either or both is missing, they need to be pulled from the DB
        # this functionality needs to be expanded eventually to take in cases where the
        # we're overlaying something with gene calls on something that doesn't have gene calls
        # if (!all(ObjectIDs[c(a1, a2)] %in% FeatureSeqs$IDs)) {
        #   # an object ID does not have a seqs present, pull them
        #   # this is not a priority so we're leaving this blank for a second
        # } else {
        #   TMPSeqs01 <- FALSE
        #   TMPSeqs02 <- FALSE
        # }
        
        # align everyone as AAs who can be, i.e. modulo of 3, is coding, etc
        # then align everyone else as nucs
        # translate the hit locations to the anchor positions
        # every hit is an anchor
        # all hits are stored in the LinkedPairs object in nucleotide space
        # as left/right bounds, orientations will be dependant upon strandedness of the genes
        
        # prepare the kmer distance stuff:
        NucDist <- vector(mode = "numeric",
                          length = nrow(PMatrix))
        
        QueryFeatureLength <- QNTCount[PMatrix[, 1L]]
        SubjectFeatureLength <- SNTCount[PMatrix[, 2L]]
        
        # grab kmer distances ahead of time because AlignPairs doesn't need to loop
        # through anything
        for (a3 in seq_along(NucDist)) {
          it1 <- PMatrix[a3, 1L]
          it2 <- PMatrix[a3, 2L]
          NucDist[a3] <- sqrt(sum((Q_NucF[it1, ] - S_NucF[it2, ])^2)) / ((sum(Q_NucF[it1, ]) + sum(S_NucF[it2, ])) / 2)
          # NucDist[a3] <- sqrt(sum((nuc1[a3, ] - nuc2[a3, ])^2)) / ((sum(nuc1[a3, ]) + sum(nuc2[a3, ])) / 2)
        }
        
        # spit out the subset vectors and logicals to correctly call both AlignPairs calls
        # and both dfs
        AASelect <- PMatrix[, 1L] %in% which(QCode & QMod) & PMatrix[, 2L] %in% which(SCode & SMod)
        NTSelect <- !AASelect
        
        df_aa_reg <- data.frame("Pattern" = Q_AA_Register[PMatrix[AASelect, 1L]],
                                "Subject" = S_AA_Register[PMatrix[AASelect, 2L]])
        df_aa <- data.frame("Pattern" = PMatrix[AASelect, 1L],
                            "Subject" = PMatrix[AASelect, 2L])
        df_nt <- data.frame("Pattern" = PMatrix[NTSelect, 1L],
                            "Subject" = PMatrix[NTSelect, 2L])
        
        # if anchor calls for inference,
        # infer from hits
        inf_fun <- function() {
          stop("inference is not currently implemented")
        }
        # if enforce, just force the terminal anchors
        enf_fun <- function() {
          if (nrow(df_aa_reg) > 0) {
            aa_anchors <- mapply(USE.NAMES = FALSE,
                                 SIMPLIFY = FALSE,
                                 FUN = function(y, z)  {
                                   cbind(matrix(data = 0L,
                                                nrow = 4),
                                         matrix(data = c(y,y,z,z),
                                                nrow = 4))
                                 },
                                 y = width(QueryAA)[df_aa_reg$Pattern] + 1L,
                                 z = width(SubjectAA)[df_aa_reg$Subject] + 1L)
          } else {
            aa_anchors <- NULL
          }
          
          if (nrow(df_nt) > 0) {
            nt_anchors <- mapply(USE.NAMES = FALSE,
                                 SIMPLIFY = FALSE,
                                 FUN = function(y, z)  {
                                   cbind(matrix(data = 0L,
                                                nrow = 4),
                                         matrix(data = c(y,y,z,z),
                                                nrow = 4))
                                 },
                                 y = width(QueryDNA)[df_nt$Pattern] + 1L,
                                 z = width(SubjectDNA)[df_nt$Subject] + 1L)
          } else {
            nt_anchors <- NULL
          }
          return(list(aa_anchors,
                      nt_anchors))
        }
        # if ignore, just don't pass an argument to pairs
        ign_fun <- function() {
          return(list(NULL,
                      NULL))
        }
        
        anchor_coords <- switch(Anchors,
                                infer = inf_fun(),
                                ignore = ign_fun(),
                                enforce = enf_fun())
        
        if (sum(AASelect) > 0) {
          if (length(anchor_coords[[1]]) > 0) {
            df_aa_reg$Position <- anchor_coords[[1]]
          } else {
            # do nothing
          }
          aapairs <- AlignPairs(pattern = QueryAA,
                                subject = SubjectAA,
                                pairs = df_aa_reg,
                                verbose = FALSE,
                                processors = Processors)
          current_local_aa_pids <- aapairs$Matches / aapairs$AlignmentLength
          current_global_aa_pids <- aapairs$Matches / pmax(width(QueryAA)[aapairs$Pattern],
                                                           width(SubjectAA)[aapairs$Subject])
          current_local_aa_scores <- aapairs$Score / aapairs$AlignmentLength
          current_global_aa_scores <- aapairs$Score / pmax(width(QueryAA)[aapairs$Pattern],
                                                           width(SubjectAA)[aapairs$Subject])
          current_global_aa_rawscore <- aapairs$Score
        } else {
          current_local_aa_pids <- numeric()
          current_global_aa_pids <- numeric()
          current_local_aa_scores <- numeric()
          current_global_aa_scores <- numeric()
          current_global_aa_rawscore <- numeric()
        }
        if (Verbose) {
          PBCount <- PBCount + 1L
          setTxtProgressBar(pb = pBar,
                            value = PBCount / Total)
        }
        
        if (sum(NTSelect) > 0) {
          if (length(anchor_coords[[2]]) > 0) {
            df_nt$Position <- anchor_coords[[2]]
          } else {
            # do nothing
          }
          ntpairs <- AlignPairs(pattern = QueryDNA,
                                subject = SubjectDNA,
                                pairs = df_nt,
                                verbose = FALSE,
                                processors = Processors)
          
          current_local_nt_pids <- ntpairs$Matches / ntpairs$AlignmentLength
          current_global_nt_pids <- ntpairs$Matches / pmax(width(QueryDNA)[ntpairs$Pattern],
                                                           width(SubjectDNA)[ntpairs$Subject])
          current_local_nt_scores <- ntpairs$Score / ntpairs$AlignmentLength
          current_global_nt_scores <- ntpairs$Score / pmax(width(QueryDNA)[ntpairs$Pattern],
                                                           width(SubjectDNA)[ntpairs$Subject])
          current_global_nt_rawscore <- ntpairs$Score
        } else {
          current_local_nt_pids <- numeric()
          current_global_nt_pids <- numeric()
          current_local_nt_scores <- numeric()
          current_global_nt_scores <- numeric()
          current_global_nt_rawscore <- numeric()
        }
        if (Verbose) {
          PBCount <- PBCount + 1L
          setTxtProgressBar(pb = pBar,
                            value = PBCount / Total)
        }
        
        vec1 <- vec2 <- vec3 <- vec4 <- vec5 <- vector(mode = "numeric",
                                                       length = nrow(PMatrix))
        vec1[AASelect] <- current_local_aa_pids
        vec1[NTSelect] <- current_local_nt_pids
        vec2[AASelect] <- current_local_aa_scores
        vec2[NTSelect] <- current_local_nt_scores
        vec3[AASelect] <- current_global_aa_pids
        vec3[NTSelect] <- current_global_nt_pids
        vec4[AASelect] <- current_global_aa_scores
        vec4[NTSelect] <- current_global_nt_scores
        vec5[AASelect] <- current_global_aa_rawscore
        vec5[NTSelect] <- current_global_nt_rawscore
        
        blockres <- BlockByRank(index1 = IMatrix[, 1L],
                                partner1 = PMatrix[, 1L],
                                index2 = IMatrix[, 2L],
                                partner2 = PMatrix[, 2L])
        
        # update blocks with the block offset counter that is initialized at
        # the beginning of the function at zero with 'block_uid'
        block_ph <- blockres$blockidmap
        w1 <- block_ph > 0
        if (any(w1)) {
          block_offset <- block_uid
          blockres$blockidmap[blockres$blockidmap > 0] <- blockres$blockidmap[blockres$blockidmap > 0] + block_offset
          block_uid <- max(blockres$blockidmap[blockres$blockidmap > 0])
        }
        
        PH[[Count]] <- data.frame("p1" = names(QueryDNA)[PMatrix[, 1]],
                                  "p2" = names(SubjectDNA)[PMatrix[, 2]],
                                  "Consensus" = diff2,
                                  "p1featurelength" = QueryFeatureLength,
                                  "p2featurelength" = SubjectFeatureLength,
                                  "blocksize" = blockres$absblocksize,
                                  "KDist" = NucDist,
                                  "TotalMatch" = TotalMatch,
                                  "MaxMatch" = MatchMax,
                                  "UniqueMatches" = UniqueMatches,
                                  "Local_PID" = vec1,
                                  "Local_Score" = vec2,
                                  "Approx_Global_PID" = vec3,
                                  "Approx_Global_Score" = vec4,
                                  "Alignment" = ifelse(test = AASelect,
                                                       yes = "AA",
                                                       no = "NT"),
                                  "Block_UID" = blockres$blockidmap,
                                  "Delta_Background" = (vec4 - diff3))
        
      } else {
        # link table is not populated
        
        PH[[Count]] <- data.frame("p1" = character(),
                                  "p2" = character(),
                                  "Consensus" = numeric(),
                                  "p1featurelength" = integer(),
                                  "p2featurelength" = integer(),
                                  "blocksize" = integer(),
                                  "KDist" = numeric(),
                                  "TotalMatch" = integer(),
                                  "MaxMatch" = integer(),
                                  "UniqueMatches" = integer(),
                                  "Local_PID" = numeric(),
                                  "Local_Score" = numeric(),
                                  "Approx_Global_PID" = numeric(),
                                  "Approx_Global_Score" = numeric(),
                                  "Alignment" = character(),
                                  "Block_UID" = integer(),
                                  "Delta_Background" = numeric())
        
        if (Verbose) {
          PBCount <- PBCount + 1L
          setTxtProgressBar(pb = pBar,
                            value = PBCount / Total)
          PBCount <- PBCount + 1L
          setTxtProgressBar(pb = pBar,
                            value = PBCount / Total)
        }
        
      } # end row check
      # Count and PBCount are unlinked,
      # iterate through both separately and correctly
      Count <- Count + 1L
      
      if (object.size(DataPool) > Storage) {
        # ok ... 
        # i need to nuke positions in the pool based on a few things:
        # i can't nuke the current a1,
        # i can't nuke the next a1
        # i can't nuke the next a2
        # return(list("a1" = a1,
        #             "a2" = a2,
        #             "pool" = DataPool))
        sw1 <- sapply(X = DataPool,
                      FUN = function(x) {
                        !is.null(x)
                      })
        sw2 <- which(sw1)
        sw2 <- sw2[!(sw2 %in% c(a1, (a1 + 1L), (a2 + 1L)))]
        if (length(sw2) > 0) {
          # bonk the first one ... this might realistically need to happen in a while loop, but for now we can live with this
          # !!! assigning NULL by default deletes the position, shortening the vector,
          # we can replace the list (the container), with another container (containing NULL) without
          # shortening the list, R inferno reference and relevant examples here:
          # https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
          DataPool[sw2[1L]] <- list(NULL)
        } else {
          # i don't know if this case can happen, but we're putting a print statement here just in case
          print("Please allocate more storage.")
        }
      }
      Prev_a1 <- a1
    } # end a2
  } # end a1
  res <- do.call(rbind,
                 PH)
  if (nrow(res) > 0) {
    # return(res)
    All_UIDs <- unique(res$Block_UID)
    res$Block_UID[res$Block_UID == -1L] <- seq(from = max(All_UIDs) + 1L,
                                               by = 1L,
                                               length.out = sum(res$Block_UID == -1L))
  }
  attr(x = res,
       which = "GeneCalls") <- GeneCalls
  attr(x = res,
       which = "KmerSize") <- KmerSize
  attr(x = res,
       which = "DefaultTranslationTable") <- DefaultTranslationTable
  class(res) <- c("data.frame",
                  "PairSummaries")
  # close pBar and return res
  if (Verbose) {
    TimeEnd <- Sys.time()
    close(pBar)
    print(TimeEnd - TimeStart)
  }
  return(res)
  
} # end function
