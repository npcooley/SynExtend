###### -- Given a database connection create a set of decoy alignments --------
# given a database
# and a genecalls list
# return a PairSummaries object of alignments made between stringsets where one
# stringset is reversed, and the other is not

###### -- NOTES ---------------------------------------------------------------
# this function is not designed to be as generalized as some others,
# so it does not have the size and efficiency controls that summarizepairs
# and nucleotideoverlaps have gone to the trouble to implement

###### -- FUNCTION ------------------------------------------------------------
# overhead checks here are ... modest, but not exhaustive

CreateDecoys <- function(DataBase01,
                         GeneCalls,
                         K_val_01,
                         K_val_02,
                         K_val_03,
                         DefaultTranslationTable = "11",
                         NT_limit = NULL,
                         AA_limit = NULL,
                         TRY_limit = NULL,
                         DecoyMethod = c("reverse", "none"),
                         Verbose = FALSE) {
  
  # DataBase01, GeneCalls, and K_val must all be supplied, error if any are missing
  # names of the genecalls object will be used to plan out the comparisons
  # to be made,
  # scroll through the upper triangle of the present names until it's either done,
  # or an early end has been reached
  
  conv_index_vals <- function(ind_hits) {
    
    hit_lengths <- vapply(X = ind_hits$Position,
                          FUN = function(x) {
                            ncol(x)
                          },
                          FUN.VALUE = vector(mode = "integer",
                                             length = 1L))
    comb_hits <- do.call(cbind,
                         ind_hits$Position)
    
    vec_set <- cbind("pattern" = rep(x = ind_hits$Pattern,
                                     times = hit_lengths),
                     "subject" = rep(x = ind_hits$Subject,
                                     times = hit_lengths))
    
    return(list("positions" = comb_hits,
                "indices" = vec_set,
                "hit_blocking" = hit_lengths))
  }
  
  adhoc_consensus <- function(pattern_left,
                              pattern_right,
                              subject_left,
                              subject_right,
                              pattern_width,
                              subject_width) {
    
    # do hits occur in similar linear positions
    res <- vector(mode = "numeric",
                  length = length(pattern_left))
    for (d1 in seq_along(pattern_left)) {
      ph1 <- pattern_left[d1] / pattern_width[d1]
      ph2 <- pattern_right[d1] / pattern_width[d1]
      
      ph3 <- subject_left[d1] / subject_width[d1]
      ph4 <- subject_right[d1] / subject_width[d1]
      
      ph5 <- ph1 - ph3
      if (ph5 < 0) {
        ph5 <- abs(ph5)
      }
      ph6 <- ph2 - ph4
      if (ph6 < 0) {
        ph6 <- abs(ph6)
      }
      
      res[d1] <- (ph5 + ph6) / 2
      
    }
    
    return(res)
  }
  
  # no overhead checking, alignment background calcs
  # sum of frequencies 1 * frequencies 2
  # the frequencies were collected by subsetting the sequences against the supplied
  # substitution matrix, so our 'checking' should have been accomplished already
  sequence_background <- function(frequencies1,
                                  frequencies2,
                                  index1,
                                  index2,
                                  substitution_matrix) {
    res <- vector(mode = "numeric",
                  length = length(index1))
    for (d1 in seq_along(index1)) {
      res[d1] <- sum(frequencies1[index1[d1], ] * t(frequencies2[index2[d1], ] * substitution_matrix))
    }
    return(res)
  }
  
  if (Verbose) {
    tstart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
    cat("\nCreating decoy alignments ... \n")
  }
  
  # get the linker ids that connects the db and the gene calls
  ObjectIDs <- names(GeneCalls)
  
  DecoyMethod <- match.arg(DecoyMethod)
  
  if (length(ObjectIDs) < 2L) {
    stop ("'GeneCalls' object must be at least length 2.")
  }
  
  # if a database is supplied, just build a generic decoy set
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
      stop("'DataBase01' is not a valid database connection.")
    }
  }
  
  # hardcode some values
  AA_matrix <- DECIPHER:::.getSubMatrix("PFASUM50")
  NT_matrix <- DECIPHER:::.nucleotideSubstitutionMatrix(2L, -1L, 1L)
  
  id_check <- dbGetQuery(conn = dbConn,
                         statement = "select distinct identifier from Seqs")$identifier
  if (!all(ObjectIDs %in% id_check)) {
    stop ("Identifiers do not appear to be shared between the database and the 'GeneCalls' object.")
  }
  
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
  
  # loop limits
  L1 <- length(ObjectIDs)
  pBar <- txtProgressBar(style = 1L)
  if (is.null(TRY_limit)) {
    PBAR <- (L1 * (L1 - 1L)) / 2L
  } else if ((is.numeric(TRY_limit) | is.integer(TRY_limit)) &
             length(TRY_limit) == 1) {
    PBAR <- TRY_limit
  } else {
    stop ("'TRY_limit' must be a numeric or integer of length 1.")
  }
  
  res <- vector(mode = "list",
                length = PBAR)
  
  # other overhead things:
  break_d1 <- FALSE
  count <- 0L
  nt_count <- 0L
  aa_count <- 0L
  
  # print("entering loops")
  # run the whole shenanigan together
  for (d1 in seq(length(ObjectIDs) - 1L)) {
    ph1 <- ObjectIDs[d1]
    # access by name, not index
    ph2 <- GeneCalls[[ph1]]
    
    ph2_nt_w <- lapply(X = ph2$Range,
                       FUN = function(x) {
                         width(x)
                       })
    ph2_nt_w <- vapply(X = ph2_nt_w,
                       FUN = function(x) {
                         sum(x)
                       },
                       FUN.VALUE = vector(mode = "integer",
                                          length = 1))
    
    ph2_coding_subset <- ph2_nt_w %% 3 == 0 & ph2$Coding
    ph2_aa_w <- ph2_nt_w[ph2_coding_subset] / 3
    # if the AAs aren't in the table, do the whole dance and translate
    if (!db_id_present[d1]) {
      curr_d1_genome <- SearchDB(dbFile = dbConn,
                                 nameBy = "description",
                                 identifier = ph1,
                                 verbose = FALSE)
      curr_d1_nucs <- FeaturesFromDF(Genome = curr_d1_genome,
                                     GeneCalls = ph2)
      curr_d1_prots <- GetTranslatedFeatures(Nucs = curr_d1_nucs,
                                             GeneCalls = ph2,
                                             DefaultTranslationTable = DefaultTranslationTable)
      # send these to the database, because why not ...
      Seqs2DB(seqs = curr_d1_prots,
              dbFile = dbConn,
              tblName = "AAs",
              type = "XStringSet",
              verbose = FALSE,
              identifier = ph1)
    } else {
      # else, just pull them
      curr_d1_genome <- SearchDB(dbFile = dbConn,
                                 nameBy = "description",
                                 identifier = ph1,
                                 verbose = FALSE)
      curr_d1_nucs <- FeaturesFromDF(Genome = curr_d1_genome,
                                     GeneCalls = ph2)
      curr_d1_prots <- SearchDB(dbFile = dbConn,
                                tblName = "AAs",
                                identifier = ph1,
                                verbose = FALSE,
                                nameBy = "description")
    }
    # reverse as few times as possible
    # this isn't that expensive, but honestly, why not...
    # curr_d1_prots <- reverse(curr_d1_prots)
    
    # enter inner loop here, same deal with the dimension 2 iterator
    for (d2 in (d1 + 1L):length(ObjectIDs)) {
      ph1 <- ObjectIDs[d2]
      # access by name, not index
      ph3 <- GeneCalls[[ph1]]
      
      ph3_nt_w <- lapply(X = ph3$Range,
                         FUN = function(x) {
                           width(x)
                         })
      ph3_nt_w <- vapply(X = ph3_nt_w,
                         FUN = function(x) {
                           sum(x)
                         },
                         FUN.VALUE = vector(mode = "integer",
                                            length = 1))
      
      ph3_coding_subset <- ph3_nt_w %% 3 == 0 & ph3$Coding
      ph3_aa_w <- ph3_nt_w[ph3_coding_subset] / 3
      # if the AAs aren't in the table, do the whole dance and translate
      if (!db_id_present[d2]) {
        curr_d2_genome <- SearchDB(dbFile = dbConn,
                                   nameBy = "description",
                                   identifier = ph1,
                                   verbose = FALSE)
        curr_d2_nucs <- FeaturesFromDF(Genome = curr_d2_genome,
                                       GeneCalls = ph3,
                                       Index = ph1)
        curr_d2_prots <- GetTranslatedFeatures(Nucs = curr_d2_nucs,
                                               GeneCalls = ph3,
                                               DefaultTranslationTable = DefaultTranslationTable)
        # send these to the database, because why not ...
        Seqs2DB(seqs = curr_d2_prots,
                dbFile = dbConn,
                tblName = "AAs",
                type = "XStringSet",
                verbose = FALSE,
                identifier = ph1)
      } else {
        curr_d2_genome <- SearchDB(dbFile = dbConn,
                                   nameBy = "description",
                                   identifier = ph1,
                                   verbose = FALSE)
        curr_d2_nucs <- FeaturesFromDF(Genome = curr_d2_genome,
                                       GeneCalls = ph3,
                                       Index = ph1)
        curr_d2_prots <- SearchDB(dbFile = dbConn,
                                  tblName = "AAs",
                                  identifier = ph1,
                                  verbose = FALSE,
                                  nameBy = "description")
      }
      # lowest loop level, work occurs here
      # randomly reverse a set before the business
      if (sample(x = 2,
                 size = 1) == 1) {
        # respect method choice, allow for more complex decoy schemes
        # in the future
        if (DecoyMethod == "reverse") {
          subject_prots <- reverse(curr_d1_prots)
          subject_nucs <- reverse(curr_d1_nucs)
        } else if (DecoyMethod == "none") {
          subject_prots <- curr_d1_prots
          subject_nucs <- curr_d1_nucs
        }
        pattern_prots <- curr_d2_prots
        pattern_nucs <- curr_d2_nucs
      } else {
        subject_prots <- curr_d1_prots
        subject_nucs <- curr_d1_nucs
        # respect method choice, allow for more complex decoy schemes
        # in the future
        if (DecoyMethod == "reverse") {
          pattern_prots <- reverse(curr_d2_prots)
          pattern_nucs <- reverse(curr_d2_nucs)
        } else if (DecoyMethod == "none") {
          pattern_prots <- curr_d2_prots
          pattern_nucs <- curr_d2_nucs
        }
      }
      
      # d1 overheads
      subject_aa_freq <- alphabetFrequency(x = subject_prots)
      subject_aa_freq <- subject_aa_freq[, colnames(AA_matrix)]
      subject_aa_freq <- subject_aa_freq / rowSums(subject_aa_freq)
      subject_nt_freq <- alphabetFrequency(x = subject_nucs)
      subject_nt_freq <- subject_nt_freq[, colnames(NT_matrix)]
      subject_nt_freq <- subject_nt_freq / rowSums(subject_nt_freq)
      subject_nt_kmer <- oligonucleotideFrequency(x = subject_nucs,
                                                  width = K_val_01,
                                                  as.prob = TRUE)
      
      # d2 overheads
      pattern_aa_freq <- alphabetFrequency(x = pattern_prots)
      pattern_aa_freq <- pattern_aa_freq[, colnames(AA_matrix)]
      pattern_aa_freq <- pattern_aa_freq / rowSums(pattern_aa_freq)
      pattern_nt_freq <- alphabetFrequency(x = pattern_nucs)
      pattern_nt_freq <- pattern_nt_freq[, colnames(NT_matrix)]
      pattern_nt_freq <- pattern_nt_freq / rowSums(pattern_nt_freq)
      pattern_nt_kmer <- oligonucleotideFrequency(x = pattern_nucs,
                                                  width = K_val_01,
                                                  as.prob = TRUE)
      
      
      # print("overhead completed")
      # search 1 direction
      # eventually we need to shift everything to only descend from these searches
      # when they're populated ...
      vals01 <- IndexSeqs(subject = subject_nucs,
                          K = K_val_03,
                          verbose = FALSE)
      
      vals02 <- IndexSeqs(subject = subject_prots,
                          K = K_val_02,
                          verbose = FALSE)
      
      vals03 <- SearchIndex(pattern = pattern_nucs,
                            invertedIndex = vals01,
                            subject = subject_nucs,
                            verbose = FALSE)
      
      if (nrow(vals03) > 0) {
        # print("entering nucs evals")
        # kmer distances
        nt_K_dists <- vector(mode = "numeric",
                             length = nrow(vals03))
        it1 <- vals03$Subject
        it2 <- vals03$Pattern
        for (d3 in seq_along(nt_K_dists)) {
          nt_K_dists[d3] <- sqrt(sum((subject_nt_kmer[it1[d3], ] - pattern_nt_kmer[it2[d3], ])^2)) / ((sum(subject_nt_kmer[it1[d3], ]) + sum(pattern_nt_kmer[it2[d3], ])) / 2)
        }
        # hits
        hit_sets <- lapply(X = vals03$Position,
                           FUN = function(x) {
                             x[2, ] - x[1, ] + 1L
                           })
        count_hits <- lengths(hit_sets)
        max_hits <- vapply(X = hit_sets,
                           FUN = function(x) {
                             max(x)
                           },
                           FUN.VALUE = vector(mode = "integer",
                                              length = 1L))
        total_hits <- vapply(X = hit_sets,
                             FUN = function(x) {
                               sum(x)
                             },
                             FUN.VALUE = vector(mode = "integer",
                                                length = 1L))
        # print("starting nucs alignments")
        # alignments
        vals05 <- AlignPairs(pattern = pattern_nucs,
                             subject = subject_nucs,
                             pairs = vals03,
                             verbose = FALSE)
        # alignment scores and pids
        nt_local_pids <- vals05$Matches / vals05$AlignmentLength
        nt_global_pids <- vals05$Matches / pmax(width(pattern_nucs)[vals05$Pattern],
                                                width(subject_nucs)[vals05$Subject])
        
        nt_local_scores <- vals05$Score / vals05$AlignmentLength
        nt_global_scores <- vals05$Score / pmax(width(pattern_nucs)[vals05$Pattern],
                                                width(subject_nucs)[vals05$Subject])
        # something else
        vals07 <- conv_index_vals(ind_hits = vals03)
        vals07$widths <- cbind("pattern" = ph3_nt_w[vals07$indices[, "pattern"]],
                               "subject" = ph2_nt_w[vals07$indices[, "subject"]])
        
        # print("getting nucs consensus scores")
        vals09 <- adhoc_consensus(pattern_left = vals07$positions[1, ],
                                  pattern_right = vals07$positions[2, ],
                                  subject_left = vals07$positions[3, ],
                                  subject_right = vals07$positions[4, ],
                                  pattern_width = vals07$widths[, "pattern"],
                                  subject_width = vals07$widths[, "subject"])
        
        # consensus is now pooled
        # ~1 == hits are all in the same relative positions in both sequences
        # ~0 == hits are in as different relative positions in both sequences as is possible
        
        vals09 <- 1 - unname(tapply(X = vals09,
                                    INDEX = rep(x = seq(length(vals07$hit_blocking)),
                                                times = vals07$hit_blocking),
                                    FUN = function(x) {
                                      mean(x)
                                    }))
        
        # print("getting nucs backgrounds")
        vals11 <- sequence_background(frequencies1 = subject_nt_freq,
                                      frequencies2 = pattern_nt_freq,
                                      index1 = vals03$Subject,
                                      index2 = vals03$Pattern,
                                      substitution_matrix = NT_matrix)
        
        nt_delta_background <- nt_global_scores - vals11
        # create the intermediate results
        nt_res <- data.frame("p1" = names(subject_nucs)[vals03$Subject],
                             "p2" = names(pattern_nucs)[vals03$Pattern],
                             "Consensus" = vals09,
                             "p1featurelength" = width(subject_nucs)[vals03$Subject],
                             "p2featurelength" = width(pattern_nucs)[vals03$Pattern],
                             "blocksize" = rep(1,
                                               nrow(vals03)),
                             "KDist" = nt_K_dists,
                             "TotalMatch" = total_hits,
                             "MaxMatch" = max_hits,
                             "UniqueMatches" = count_hits,
                             "Local_PID" = nt_local_pids,
                             "Local_Score" = nt_local_scores,
                             "Approx_Global_PID" = nt_global_pids,
                             "Approx_Global_Score" = nt_global_scores,
                             "Alignment" = rep("NT",
                                               nrow(vals03)),
                             "Block_UID" = seq(nrow(vals03)),
                             "Delta_Background" = nt_delta_background)
        # print("nucs completed")
      } else {
        # no search hits
        nt_res <- data.frame("p1" = character(),
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
      }
      
      vals04 <- SearchIndex(pattern = pattern_prots,
                            invertedIndex = vals02,
                            subject = subject_prots,
                            verbose = FALSE)
      
      if (nrow(vals04) > 0) {
        # print("entering aa evals")
        aa_K_dists <- vector(mode = "numeric",
                             length = nrow(vals04))
        it1 <- vals04$Subject
        it2 <- vals04$Pattern
        string_ph1 <- subject_nt_kmer[ph2_coding_subset, , drop = FALSE]
        string_ph2 <- pattern_nt_kmer[ph3_coding_subset, , drop = FALSE]
        for (d3 in seq_along(aa_K_dists)) {
          aa_K_dists[d3] <- sqrt(sum((string_ph1[it1[d3], ] - string_ph2[it2[d3], ])^2)) / ((sum(string_ph1[it1[d3], ]) + sum(string_ph2[it2[d3], ])) / 2)
        }
        # hits
        hit_sets <- lapply(X = vals04$Position,
                           FUN = function(x) {
                             x[2, ] - x[1, ] + 1L
                           })
        count_hits <- lengths(hit_sets)
        max_hits <- vapply(X = hit_sets,
                           FUN = function(x) {
                             max(x)
                           },
                           FUN.VALUE = vector(mode = "integer",
                                              length = 1L))
        total_hits <- vapply(X = hit_sets,
                             FUN = function(x) {
                               sum(x)
                             },
                             FUN.VALUE = vector(mode = "integer",
                                                length = 1L))
        # print("starting aa alignments")
        vals06 <- AlignPairs(pattern = pattern_prots,
                             subject = subject_prots,
                             pairs = vals04,
                             verbose = FALSE)
        
        # alignment scores and pids
        aa_local_pids <- vals06$Matches / vals06$AlignmentLength
        aa_global_pids <- vals06$Matches / pmax(width(pattern_prots)[vals06$Pattern],
                                                width(subject_prots)[vals06$Subject])
        
        aa_local_scores <- vals06$Score / vals06$AlignmentLength
        aa_global_scores <- vals06$Score / pmax(width(pattern_prots)[vals06$Pattern],
                                                width(subject_prots)[vals06$Subject])
        
        vals08 <- conv_index_vals(ind_hits = vals04)
        vals08$widths <- cbind("pattern" = ph3_aa_w[vals08$indices[, "pattern"]],
                               "subject" = ph2_aa_w[vals08$indices[, "subject"]])
        # print("getting aa consensus scores")
        vals10 <- adhoc_consensus(pattern_left = vals08$positions[1, ],
                                  pattern_right = vals08$positions[2, ],
                                  subject_left = vals08$positions[3, ],
                                  subject_right = vals08$positions[4, ],
                                  pattern_width = vals08$widths[, "pattern"],
                                  subject_width = vals08$widths[, "subject"])
        
        vals10 <- 1 - unname(tapply(X = vals10,
                                    INDEX = rep(x = seq(length(vals08$hit_blocking)),
                                                times = vals08$hit_blocking),
                                    FUN = function(x) {
                                      mean(x)
                                    }))
        # print("getting aa backgrounds")
        vals12 <- sequence_background(frequencies1 = subject_aa_freq,
                                      frequencies2 = pattern_aa_freq,
                                      index1 = vals04$Subject,
                                      index2 = vals04$Pattern,
                                      substitution_matrix = AA_matrix)
        aa_delta_background <- aa_global_scores - vals12
        # create the intermediate results
        # space is in the genomic context, not the AA search space ...
        aa_res <- data.frame("p1" = names(subject_prots)[vals04$Subject],
                             "p2" = names(pattern_prots)[vals04$Pattern],
                             "Consensus" = vals10,
                             "p1featurelength" = width(subject_prots)[vals04$Subject] * 3,
                             "p2featurelength" = width(pattern_prots)[vals04$Pattern] * 3,
                             "blocksize" = rep(1,
                                               nrow(vals04)),
                             "KDist" = aa_K_dists,
                             "TotalMatch" = total_hits * 3,
                             "MaxMatch" = max_hits * 3,
                             "UniqueMatches" = count_hits,
                             "Local_PID" = aa_local_pids,
                             "Local_Score" = aa_local_scores,
                             "Approx_Global_PID" = aa_global_pids,
                             "Approx_Global_Score" = aa_global_scores,
                             "Alignment" = rep("AA",
                                               nrow(vals04)),
                             "Block_UID" = seq(nrow(vals04)),
                             "Delta_Background" = aa_delta_background)
        # print("aas complete")
      } else {
        aa_res <- data.frame("p1" = character(),
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
      }
      
      # drop any nt alignments that are represented in the aa set
      check1 <- paste(nt_res$p1,
                      nt_res$p2,
                      sep = "_")
      check2 <- paste(aa_res$p1,
                      aa_res$p2,
                      sep = "_")
      if (length(check1) > 0 &
          length(check2) > 0) {
        if (any(check1 %in% check2)) {
          nt_res <- nt_res[!(check1 %in% check2), ]
        }
      }
      
      # this can be zero rows!
      int_res <- rbind(nt_res,
                       aa_res)
      res[[count + 1L]] <- int_res
      
      curr_total_counts <- table(factor(x = int_res$Alignment,
                                        levels = c("AA", "NT")))
      
      aa_count <- aa_count + unname(curr_total_counts["AA"])
      nt_count <- nt_count + unname(curr_total_counts["NT"])
      
      # count starts at zero
      # and PBAR is set based on TRY_limit if it is not null
      count <- count + 1L
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = count / PBAR)
      }
      
      if (!is.null(AA_limit)) {
        if (aa_count > AA_limit) {
          break_d1 <- TRUE
          break
        }
      }
      if (!is.null(NT_limit)) {
        if (nt_count > NT_limit) {
          break_d1 <- TRUE
          break
        }
      }
      if (count >= PBAR) {
        break_d1 <- TRUE
        break
      }
    } # end d2 block
    if (break_d1) {
      break
    }
  } # end d1 block
  
  if (Verbose) {
    close(pBar)
  }
  
  null_check <- vapply(X = res,
                       FUN = function(x) {
                         is.null(x)
                       },
                       FUN.VALUE = vector(mode = "logical",
                                          length = 1))
  res <- res[!null_check]
  if (length(res) > 1) {
    res <- do.call(rbind,
                   res)
  } else {
    res <- res[[1]]
  }
  class(res) <- c("data.frame",
                  "PairSummaries")
  attr(x = res,
       which = "GeneCalls") <- GeneCalls
  attr(x = res,
       which = "KmerSize") <- K_val_01
  res$Block_UID <- seq(nrow(res))
  return(res)
  
}


