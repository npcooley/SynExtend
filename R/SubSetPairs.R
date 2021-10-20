###### -- SubSetPairs ---------------------------------------------------------
# author: nicholas cooley
# email: npc19@pitt.edu

# take in a PairSummaries object and remove predicted pairs based on user inputs


SubSetPairs <- function(CurrentPairs,
                        UserThresholds,
                        RejectCompetitors = TRUE,
                        RejectionCriteria = "PID",
                        WinnersOnly = TRUE,
                        Verbose = FALSE) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  ###### -- overhead checking -------------------------------------------------
  # CurrentPairs must be a PairSummaries object and must be larger than 1 row
  # user thresholds must be a named vector of numerics that must make sense
  # should be ordered by priority
  if (!missing(UserThresholds)) {
    # user thresholds are present
    if (all(names(UserThresholds) %in% colnames(CurrentPairs))) {
      USETHRESH <- TRUE
    } else {
      stop ("Supplied thresholds do not match pair statistics.")
    }
  } else {
    USETHRESH <- FALSE
  }
  # Competitor rejection left as `TRUE` will reject all competitors
  # can be a numeric between zero and 1, to only reject competitors below a default threshold
  if (is.logical(RejectCompetitors) &
      !is.na(RejectCompetitors)) {
    CompRejection <- RejectCompetitors
    RejectCompetitors <- max(CurrentPairs[, RejectionCriteria])
  } else if (is.integer(RejectCompetitors) |
             is.numeric(RejectCompetitors)) {
    CompRejection <- TRUE
    # do nothing, this is correct
  } else {
    stop ("RejectCompetitors argument must be either a logical to specify default behavior, or a numeric / integer.")
  }
  # RejectionCriteria must be a column name present in CurrentPairs
  if (!(RejectionCriteria %in% colnames(CurrentPairs))) {
    stop ("RejectionCriteria must specify an existing column name.")
  }
  
  
  ###### -- competitor rejection ----------------------------------------------
  
  keep <- rep(TRUE,
              nrow(CurrentPairs))
  
  if (CompRejection) {
    
    # build indices matrices to work through
    if (Verbose) {
      cat("\nConverting node IDs.\n")
    }
    gmat <- cbind(as.integer(unlist(regmatches(x = CurrentPairs$p1,
                                               m = gregexpr(pattern = "(^[0-9]+)(?=_)",
                                                            text = CurrentPairs$p1,
                                                            perl = TRUE)))),
                  as.integer(unlist(regmatches(x = CurrentPairs$p2,
                                               m = gregexpr(pattern = "(^[0-9]+)(?=_)",
                                                            text = CurrentPairs$p2,
                                                            perl = TRUE)))))
    imat <- cbind(as.integer(unlist(regmatches(x = CurrentPairs$p1,
                                               m = gregexpr(pattern = "(?<=_)([0-9]+)(?=_)",
                                                            text = CurrentPairs$p1,
                                                            perl = TRUE)))),
                  as.integer(unlist(regmatches(x = CurrentPairs$p2,
                                               m = gregexpr(pattern = "(?<=_)([0-9]+)(?=_)",
                                                            text = CurrentPairs$p2,
                                                            perl = TRUE)))))
    fmat <- cbind(as.integer(unlist(regmatches(x = CurrentPairs$p1,
                                               m = gregexpr(pattern = "(?<=_)([0-9]+$)",
                                                            text = CurrentPairs$p1,
                                                            perl = TRUE)))),
                  as.integer(unlist(regmatches(x = CurrentPairs$p2,
                                               m = gregexpr(pattern = "(?<=_)([0-9]+$)",
                                                            text = CurrentPairs$p2,
                                                            perl = TRUE)))))
    
    if (Verbose) {
      cat("\nFinding conflicts.\n")
    }
    
    ugm <- unique(gmat)
    L <- nrow(ugm)
    
    # loop priority:
    # first through the unique DB identifier combinations
    # then through unique index combinations within DB identifier combos
    
    for (m1 in seq_len(nrow(ugm))) {
      w1 <- which(gmat[, 1L] == ugm[m1, 1L] & gmat[, 2L] == ugm[m1, 2L])
      CurrentIndices <- imat[w1, ]
      
      uim <- unique(CurrentIndices)
      
      for (m2 in seq_len(nrow(uim))) {
        w2 <- which(imat[w1, 1L] == uim[m2, 1L] & imat[w1, 2L] == uim[m2, 2L])
        
        # ask which feature IDs are repeated in both columns
        repeatedp1 <- table(fmat[w1[w2], 1L])
        repeatedp2 <- table(fmat[w1[w2], 2L])
        checkp1 <- as.integer(names(repeatedp1[repeatedp1 > 1L]))
        checkp2 <- as.integer(names(repeatedp2[repeatedp2 > 1L]))
        
        # if there are p1 indices to resolve, resolve them
        if (length(checkp1) > 0) {
          for (m3 in seq_along(checkp1)) {
            w3 <- which(fmat[w1[w2], 1L] == checkp1[m3])
            w4 <- CurrentPairs[w1[w2[w3]], RejectionCriteria]
            w5a <- which.max(w4)
            w5b <- which(w4 >= RejectCompetitors)
            w6 <- unique(c(w5a, w5b))
            # return(list(w1,
            #             w2,
            #             w3,
            #             w4,
            #             w5a,
            #             w5b,
            #             w6))
            w6 <- w1[w2[w3[-w6]]]
            keep[w6] <- FALSE
          }
        }
        
        # if there are p2 indices to resolve, resolve them
        if (length(checkp2) > 0) {
          for (m3 in seq_along(checkp1)) {
            w3 <- which(fmat[w1[w2], 2L] == checkp2[m3])
            w4 <- CurrentPairs[w1[w2[w3]], RejectionCriteria]
            w5a <- which.max(w4)
            w5b <- which(w4 >= RejectCompetitors)
            w6 <- unique(c(w5a, w5b))
            w6 <- w1[w2[w3[-w6]]]
            keep[w6] <- FALSE
          }
        }
        
      } # end loop through current index combinations
      
      if (Verbose) {
        setTxtProgressBar(pb = pBar,
                          value = m1 / L)
      }
      
    } # end loop through DB identifier combinations
  } # end if statement for competitor rejection
  
  ###### -- user supplied thresholds on link statistics -----------------------
  
  if (USETHRESH) {
    
    cat("\nApplying user thresholds.\n")
    
    for (m1 in seq_along(UserThresholds)) {
      # whoever is below the user specified threshold for the given statistic
      w <- CurrentPairs[, names(UserThresholds)[m1]] < UserThresholds[m1]
      # is set to rejection
      keep[w] <- FALSE
    }
  }
  
  
  ###### -- return object -----------------------------------------------------
  
  EditedPairs <- CurrentPairs[keep, ]
  RejectedPairs <- CurrentPairs[!keep, ]
  
  if (Verbose) {
    cat("\n")
    TimeEnd <- Sys.time()
    print(TimeEnd - TimeStart)
  }
  
  if (WinnersOnly) {
    return(EditedPairs)
  } else {
    return(list(EditedPairs,
                RejectedPairs))
  }
  
  
}
