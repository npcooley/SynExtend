###### -- find co-localized sets ----------------------------------------------

###### -- notes ---------------------------------------------------------------
# given a genecalls list, and a detected sets / clusters object mangled from
# exolabel, return a list of data.frames, where each list position is a
# co-localized group of co-occuring sets
# each column is named after the genome index that the column population is from
# each column is populated by the row indices of the feature members
# each row represents a source cluster identified by exolabel
# clusters *must* be single copy

# if allowed_miss is greater than 0 or allow_subset is TRUE (neither are hooked up yet)
# some positions may be NA, signifying a cluster was co-localized with the other
# co-occuring clusters minus n members, or clusters that had a member that was
# not co-localized at all were allowed to be re-evaluated within the genome index
# grouping where members were observed to be co-localized

FindColocalSets <- function(genecalls,
                            detected_sets,
                            allowed_miss = 0,
                            allow_subset = FALSE,
                            complete_set = FALSE,
                            max_gap = 5000) {
  
  if (allow_subset) {
    warning("co-localized subsets are not yet supported")
  }
  if (allowed_miss > 0) {
    warning("partial co-localization is not yet supported")
  }
  
  if (!is(object = genecalls,
          class2 = "list")) {
    stop("genecalls must be a list of 'DataFrame's")
  }
  if (any(!vapply(X = genecalls,
                  FUN = function(x) {
                    nrow(x) > 0
                  },
                  FUN.VALUE = vector(mode = "logical",
                                     length = 1)))) {
    stop("genecalls must be a list of 'DataFrame's that all have greater than zero rows.")
  }
  if (any(!vapply(X = genecalls,
                  FUN = function(x) {
                    is(object = x,
                       class2 = "DataFrame")
                  },
                  FUN.VALUE = vector(mode = "logical",
                                     length = 1)))) {
    stop("genecalls mus be a list of 'DataFrame's")
  }
  
  # step one, break up sets
  # order by integer value of genome index - this may be a source of fragility
  # if i ever allow non-integer index genome identifiers
  bsets <- lapply(X = detected_sets,
                  FUN = function(x) {
                    y <- strsplit(x = x,
                                  split = "_",
                                  fixed = TRUE)
                    y <- do.call(rbind,
                                 y)
                    y <- matrix(data = as.integer(y),
                                nrow = nrow(y))
                    o1 <- order(y[, 1, drop = TRUE])
                    y <- y[o1, , drop = FALSE]
                    return(y)
                  })
  # true if a value is duplicated!
  validate_sc <- vapply(X = bsets,
                        FUN = function(x) {
                          any(duplicated(x[, 1, drop = TRUE]))
                        },
                        FUN.VALUE = vector(mode = "logical",
                                           length = 1))
  if (any(validate_sc)) {
    stop("supplied sets must represent single-copy sets")
  }
  validate_feature_ids <- vapply(X = bsets,
                                 FUN = function(x) {
                                   ncol(x) == 3L
                                 },
                                 FUN.VALUE = vector(mode = "logical",
                                                    length = 1))
  if (any(!validate_feature_ids)) {
    stop("supplied feature IDs appear malformed")
  }
  # a list where each position is a set of contig names
  contig_pops <- lapply(X = bsets,
                        FUN = function(x) {
                          mapply(SIMPLIFY = TRUE,
                                 USE.NAMES = FALSE,
                                 FUN = function(z, y) {
                                   genecalls[[z]]$Contig[y]
                                 },
                                 z = as.character(x[, 1]),
                                 y = x[, 3])
                        })
  # this is (almost) solely for index matching
  source_pops <- lapply(X = bsets,
                        FUN = function(x) {
                          x[, 1, drop = TRUE]
                        })
  
  # step two manage the genecall names
  u_genome_indices <- list(sort(as.integer(names(genecalls))))
  if (any(is.na(u_genome_indices[[1]]))) {
    stop("workflow currently only supports integer identifiers for genomic data")
  }
  
  consistency_check <- unique(unlist(source_pops))
  if (any(!(consistency_check %in% unlist(u_genome_indices)))) {
    stop("a genomic index is included in 'detected_sets' that does not have an equivalent in the supplied 'genecalls'")
  }
  
  # create our candidate sets
  if (complete_set) {
    # logical set, do you have every genome index present in the user submitted
    # set (once)
    candidate_set <- source_pops %in% u_genome_indices
    # convert to a list
    candidate_set <- list(which(candidate_set))
  } else {
    # match to the first instance of your unique set that appears in the list
    candidate_set <- match(x = source_pops,
                           table = source_pops)
    # convert to a list
    candidate_set <- tapply(X = seq_along(candidate_set),
                            INDEX = candidate_set,
                            FUN = c)
    candidate_set <- candidate_set[lengths(candidate_set) > 1]
    candidate_set <- unname(candidate_set)
  }
  # candidate set is a list of index positions that reference both `bsets` and
  # `contig_pops`
  # return(candidate_set)
  
  # build our search spaces(s)
  set_holder <- vector(mode = "list",
                       length = length(candidate_set))
  
  for (d1 in seq_along(candidate_set)) {
    
    # first check is the contig name -- subset the current candidate set
    # to only candidates that are not just in the correct genome, but in the
    # correct contig,
    # i.e. contigs[[x]][1], contigs[[x+1]][1], and contigs[[x+2]][[1]]
    # must all be on the same contig
    current_candidates <- candidate_set[[d1]]
    current_contig_set <- contig_pops[current_candidates]
    u_contig_sets <- unique(current_contig_set)
    if (length(u_contig_sets) == 1) {
      # nothing to do here
      ph2 <- list(current_candidates)
    } else {
      # split, drop length 1 groupings
      ph1 <- match(x = current_contig_set,
                   table = u_contig_sets)
      ph2 <- unname(split(x = current_candidates,
                          f = ph1))
      ph2 <- ph2[lengths(ph2) > 1]
    }
    set_holder[[d1]] <- ph2
  } # end d1
  set_holder <- unlist(set_holder,
                       recursive = FALSE)
  
  # return(set_holder)
  # each list position now should represent a unique candidate set where
  # candidate colocal sets are now guaranteed to share contigs, so we should
  # now just need to check bounds
  # list positions are populated by integer vals that are index positions
  # in bsets (the lists of broke up feature id values)
  candidate_evaluations <- vector(mode = "list",
                                  length = length(set_holder))
  for (d1 in seq_along(set_holder)) {
    target_positions <- set_holder[[d1]]
    target_gnm <- lapply(X = bsets[target_positions],
                         FUN = function(x) {
                           as.character(x[, 1, drop = TRUE])
                         })
    target_gnm <- unlist(target_gnm)
    target_row <- lapply(X = bsets[target_positions],
                         FUN = function(x) {
                           x[, 3, drop = TRUE]
                         })
    target_row <- unlist(target_row)
    bounds_l <- mapply(USE.NAMES = FALSE,
                       SIMPLIFY = TRUE,
                       FUN = function(x, y) {
                         genecalls[[x]]$Start[y]
                       },
                       x = target_gnm,
                       y = target_row)
    bounds_r <- mapply(USE.NAMES = FALSE,
                       SIMPLIFY = TRUE,
                       FUN = function(x, y) {
                         genecalls[[x]]$Stop[y]
                       },
                       x = target_gnm,
                       y = target_row)
    evaluate_within <- split(x = data.frame("feature_index" = target_row,
                                            "rb" = bounds_r,
                                            "lb" = bounds_l),
                             f = target_gnm)
    localize_key <- vector(mode = "list",
                           length = length(evaluate_within))
    
    # this is *relatively* simplistic and may not handle nested genes as 
    # intended
    for (d2 in seq_along(evaluate_within)) {
      ph1 <- nrow(evaluate_within[[d2]])
      o1 <- order(evaluate_within[[d2]]$rb)
      ph2 <- seq(from = 1L,
                 to = ph1 - 1L,
                 by = 1L)
      ph3 <- seq(from = 2L,
                 to = ph1,
                 by = 1L)
      # this is the gap size from the left bound of x + 1, to the right
      # bound of x
      # calculate the bound gaps on the sorted set
      ph4 <- evaluate_within[[d2]]$lb[o1][ph3] - evaluate_within[[d2]]$rb[o1][ph2] - 1L
      ph5 <- ph4 <= max_gap
      colocal_index <- rep(-1L,
                           nrow(evaluate_within[[d2]]))
      if (any(ph5)) {
        counter <- 1L
        for (d3 in seq_along(ph5)) {
          if (ph5[d3]) {
            colocal_index[c(d3, d3 + 1)] <- counter
          } else {
            counter <- counter + 1L
          }
        }
        # rearrage back to the natural ordering of bsets
        colocal_index[o1] <- colocal_index
      } else {
        # no candidates were co-localized within the co-occuring sets
      }
      localize_key[[d2]] <- colocal_index
    } # end d2
    # localize_key is a list of vectors, each the same length - the length
    # of the candidate evaluation grouping - where any un-localized feature
    # is represented with a -1, and every localized feature set is represented
    # with shared a positive integer
    # each integer grouping is unique to each list position
    # [[1]] -1, 1, 1, -1, 2, 2, 3, 3
    # and
    # [[2]] 1, 1, -1, 2, 2, 2, 3, 3
    # would imply the localized groups, and the combinations would imply where
    # those localizations are c-localizations
    # i.e. the 5th and 6th elements are co-localized, and so are the 7th and 8th
    # localize_key itself is the length of the set within which co-occurance
    # is being tested for co-localization
    localize_key <- do.call(rbind,
                            localize_key)
    
    r1 <- ncol(localize_key)
    r2 <- r1 - 1
    curr_res <- rep(0,
                    r1)
    # allow_miss and allow_subsetting need to be managed here
    # i.e. my set is all co-localized with the exception of ... something
    # or
    # these candidate clusters where all but x cluster member were co-localized
    # can be retested in the genome subset minus the co-occuring only member(s)
    # they are not currently set up though
    
    # this is speculative and not hooked up yet
    leftover <- rep(FALSE,
                    r1)
    
    for (d2 in seq_len(r2)) {
      for (d3 in ((d2 + 1L):r1)) {
        if (all(localize_key[, d2] == localize_key[, d3]) &
            all(localize_key[, d2] > 0)) {
          if (curr_res[d2] == 0) {
            curr_res[c(d2, d3)] <- d2
          } else {
            curr_res[d3] <- curr_res[d2]
          }
        }
      }
    } # end d2, upper triangle scroll
    # this is speculative, and not hooked up yet
    if (any(curr_res == 0)) {
      leftover[curr_res == 0] <- TRUE
    }
    candidate_evaluations[[d1]] <- split(x = set_holder[[d1]],
                                         f = curr_res)
  } # end d1
  candidate_evaluations <- unlist(candidate_evaluations,
                                  recursive = FALSE)
  candidate_evaluations <- candidate_evaluations[!(names(candidate_evaluations) %in% "0")]
  if (length(candidate_evaluations) == 0) {
    return(list())
  } else {
    names(candidate_evaluations) <- NULL
    return(candidate_evaluations)
  }
  
}

