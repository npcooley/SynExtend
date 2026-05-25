###### -- square a gff by feature types for quick access ----------------------
# take in a gff file, and return a DataFrame with columns that allow quick
# parsing and gene transformations and things

###### -- NOTES ---------------------------------------------------------------
# depends on `FindSets` function in SynExtend
# 
# will always adjust feature bounds that extend over contig bounds
# because the goal is work-able information
# the only other option is to just nuke features that span contig bounds
# current paradigm will sometimes be a problem, but only in cases where folks
# are being loosy-goosy with their feature assignments

SquaregffBy <- function(gff_object,
                        collect_by = c("gene",
                                       "pseudogene",
                                       "ncRNA_gene",
                                       "tRNA_gene"),
                        verbose = FALSE) {
  # overhead checks
  if (!is(object = gff_object,
          class2 = "GRanges")) {
    stop("object must be a GRanges object")
  }
  if (length(gff_object) < 1L) {
    stop("object has no features")
  }
  if (verbose) {
    tstart <- Sys.time()
  }
  
  # data preparation
  # get the parent child relationships
  parent_node <- gff_object$Parent
  child_node <- gff_object$ID
  
  # deal with the CharacterList vector
  # add a self-loop to square off the table
  parent_node[lengths(parent_node) < 1L] <- child_node[lengths(parent_node) < 1L]
  parent_node <- unlist(parent_node)
  
  P1Ids <- unique(parent_node)
  P2Ids <- unique(child_node)
  AllIds <- unique(c(P1Ids,
                     P2Ids))
  F1 <- factor(parent_node,
               levels = AllIds)
  F2 <- factor(child_node,
               levels = AllIds)
  
  FI1 <- as.integer(F1)
  FI2 <- as.integer(F2)
  FC1 <- as.character(F1)
  FC2 <- as.character(F2)
  
  FInts <- c(FI1, FI2)
  FChars <- c(FC1, FC2)
  Key1 <- !duplicated(FInts)
  IntMap <- cbind("FactorRep" = FInts[Key1],
                  "UniqueIDs" = FChars[Key1])
  
  # run FindSets
  IntRes <- FindSets(p1 = FI1,
                     p2 = FI2)
  
  # all nodes
  Nodes <- IntRes[, 1L, drop = TRUE]
  # all representatives
  Reps <- IntRes[, 2L, drop = TRUE]
  # set the order to the unique reps
  AllIds <- IntMap[match(x = IntMap[, 1L],
                         table = Nodes), 2L]
  
  # get the disjoint sets as a vector
  UClusts <- split(x = AllIds,
                   f = Reps)
  # subset the gff to parent features of the appropriate type and their children
  w1 <- gff_object$type %in% collect_by
  w2 <- gff_object$ID[w1]
  if (length(w2) < 1L) {
    stop("no features appear to match those enumerated by 'select_from' argument")
  }
  w3 <- rep(x = seq(length(UClusts)),
            times = lengths(UClusts))
  w4 <- unname(unlist(UClusts))
  
  m1 <- match(x = child_node,
              table = w4)
  m2 <- w3[m1]
  
  # get region max widths:
  w5 <- which(gff_object$type == "region")
  if (length(w5) < 1L) {
    stop("gff does not appear to have expected 'region' features")
  }
  
  contig_bounds <- end(gff_object@ranges[w5])
  names(contig_bounds) <- as.character(gff_object@seqnames[w5])
  
  obj_contigs <- as.character(gff_object@seqnames)
  obj_indices <- as.integer(gff_object@seqnames)
  obj_ranges <- gff_object@ranges
  # need to support:
  # + == positive strand == 0 in decipher scheme
  # - == negative strand == 1 in decipher scheme
  # . == strand not relevant / strand unknown
  # * == both strands relevant
  # ? == strand unknown
  # for now I'm just going to default these to the forward direction because it
  # makes my life easier, this is an expedient decision, but probably not the
  # best decision
  obj_strand <- as.character(gff_object@strand)
  # we're defaulting to positive for ambiguous stuff, so negative is the test
  # and failure returns positive
  obj_strand_rep <- ifelse(test = obj_strand == "-",
                           yes = 1L,
                           no = 0L)
  obj_type <- gff_object$type
  obj_transl <- gff_object$transl_table
  obj_product <- gff_object$product
  obj_pseudo <- gff_object$pseudo
  obj_dbxref <- gff_object$Dbxref
  obj_notes <- gff_object$Notes
  
  u1 <- unique(m2)
  l1 <- length(u1)
  
  # collect the following:
  # Index - x
  # Strand - x
  # Start - x
  # Stop - x
  # Type - x
  # ID - x
  # Range - x
  # Product - x
  # Coding - x
  # Translation Table - x
  # Contig - x
  # and a selection logical
  idx <- strd <- strt <- stp <- vector(mode = "integer",
                                       length = l1)
  typ <- idtfr <- prdct <- contig <- vector(mode = "character",
                                            length = l1)
  trtbl <- rep(NA_character_,
               l1)
  rng <- dbxref <- notes <- vector(mode = "list",
                                   length = l1)
  cdg <- slct <- vector(mode = "logical",
                        length = l1)
  
  if (verbose) {
    pBar <- txtProgressBar(style = 1)
    PBAR <- l1
    cat("parsing parent-child heirarchies...\n")
  }
  
  for (a1 in seq_along(u1)) {
    # print(a1)
    curr_rows <- which(m2 == u1[a1])
    # block 1 is just the subset of the gff object that matches the feature id heirarchy
    # b1 <- gff_object[curr_rows]
    curr_type <- obj_type[curr_rows]
    curr_transl <- obj_transl[curr_rows]
    curr_product <- obj_product[curr_rows]
    curr_pseudo <- obj_pseudo[curr_rows]
    curr_dbxref <- obj_dbxref[curr_rows]
    curr_ids <- child_node[curr_rows]
    curr_notes <- obj_notes[curr_rows]
    
    if (any(curr_type %in% collect_by)) {
      # flip the selection logical
      slct[a1] <- TRUE
    } else {
      # if our block doesn't contain what we want we skip
      next
    }
    # just need the first
    idx[a1] <- obj_indices[curr_rows][1L]
    strd[a1] <- obj_strand_rep[curr_rows][1L]
    contig[a1] <- obj_contigs[curr_rows][1L]
    # block 2 is the iranges of the feature id heirarchy
    b2 <- obj_ranges[curr_rows]
    
    # get the left and right bounds, feature parent should encompass the entire
    # feature
    # Start -- left bound
    strt[a1] <- b2[1]@start
    # Stop -- right bound
    stp[a1] <- strt[a1] + b2[1L]@width - 1L
    
    # Range is either just the left and right bounds of the parent
    # or the IRanges for the CDSs
    if (any(curr_type == "CDS")) {
      u_cds <- unique(curr_ids[curr_type == "CDS"])
      # the iranges for the CDSs
      rng[[a1]] <- b2[curr_ids == u_cds[1L]]
      # reorder everything to be left-right, i did this before, and it's not
      # necessarily convention, but it's how my translation routines later on expect
      # to see this data
      o1 <- order(start(rng[[a1]]))
      rng[[a1]] <- rng[[a1]][o1]
      # the first product line for the CDSs
      prdct[a1] <- curr_product[curr_ids == u_cds[1L]][1L]
      cdg[a1] <- TRUE
      ph <- curr_transl
      if (all(is.na(ph))) {
        # no default setting at this stage
        # trtbl[a1] <- TranslationTable
      } else {
        trtbl[a1] <- ph[!is.na(ph)][1L]
      }
    } else {
      # irange of the parent feature
      rng[[a1]] <- b2[1L]
      # first product line that is not NA
      ph <- curr_product[!is.na(curr_product)]
      if (length(ph) < 1L) {
        prdct[a1] <- ""
      } else {
        prdct[a1] <- ph[1L]
      }
      
    }
    # can grab the first entry, no NAs expected
    idtfr[a1] <- curr_ids[1L]
    # type is asking for pseudo, which is *currently a character*!
    if (all(is.na(curr_pseudo))) {
      typ[a1] <- "gene"
    } else {
      typ[a1] <- "pseudogene"
    }
    
    # grab all the dbxrefs that are available, unique them
    ph <- curr_dbxref
    ph <- unique(unlist(ph))
    dbxref[[a1]] <- ph
    # grab all the notes that are available, unique them
    ph <- curr_notes
    ph <- unique(unlist(ph))
    notes[[a1]] <- ph
    
    if (verbose) {
      setTxtProgressBar(pb = pBar,
                        value = a1 / PBAR)
    }
  } # end a1 loop
  if (verbose) {
    close(pBar)
  }
  
  # return(list(idx,
  #             strd,
  #             strt,
  #             stp,
  #             typ,
  #             idtfr,
  #             rng,
  #             prdct,
  #             cdg,
  #             trtbl,
  #             contig,
  #             slct))
  res <- DataFrame("Index" = idx[slct],
                   "Strand" = strd[slct],
                   "Start" = strt[slct],
                   "Stop" = stp[slct],
                   "Type" = typ[slct],
                   "ID" = idtfr[slct],
                   "Range" = IRangesList(rng[slct]),
                   "Product" = prdct[slct],
                   "Coding" = cdg[slct],
                   "Translation_Table" = trtbl[slct],
                   "Contig" = contig[slct],
                   "Dbxref" = CharacterList(dbxref[slct]),
                   "Notes" = CharacterList(notes[slct]))
  
  if (verbose) {
    pBar <- txtProgressBar(style = 1)
    cat("checking feature / contig bound conflicts...\n")
  }
  PBAR <- sum(slct)
  
  for (a1 in seq_len(PBAR)) {
    # rewrite any ranges where bounds extend over the end of an index
    ph1 <- unlist(end(res$Range[[a1]]))
    ph2 <- unlist(start(res$Range[[a1]]))
    if (any(ph1 > contig_bounds[res$Contig[a1]])) {
      # at least one bound extends over the end of the index
      B_Starts <- ph2 # Bound starts
      B_Ends <- ph1 # Bound ends
      F_Break <- contig_bounds[res$Contig[a1]] # FASTA break
      # do this with a while loop
      RemainingBounds <- TRUE
      C_Counts <- 1L
      N_Counts <- 1L
      N_Starts <- rep(NA_integer_,
                      length(B_Starts) * 2L)
      N_Ends <- rep(NA_integer_,
                    length(B_Ends) * 2L)
      while(RemainingBounds) {
        if (F_Break >= B_Ends[C_Counts] &
            C_Counts == length(B_Ends)) {
          # print("break after, no further evaluations")
          # Case 1:
          # Break occurs after this bound set
          # This bound set is the only remaining set
          # add bounds to new vectors
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- B_Ends[C_Counts]
          RemainingBounds <- FALSE
        } else if (F_Break < B_Starts[C_Counts] &
                   C_Counts == length(B_Starts)) {
          # print("break prior, no further evaluations")
          # Case 2:
          # Break occurs before this bound set
          # This bound set is the only remaining set
          # adjust current bound set
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts] - F_Break
          N_Ends[N_Counts] <- B_Ends[C_Counts] - F_Break
          RemainingBounds <- FALSE
        } else if (F_Break >= B_Starts[C_Counts] &
                   F_Break < B_Ends[C_Counts] &
                   C_Counts == length(B_Starts)) {
          # print("break interior, no further evaluations")
          # Case 3:
          # Break occurs inside current bound set
          # this is the only remaining bound set
          # add set with adjusted current end
          # add new set
          # exit
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- F_Break
          N_Starts[N_Counts + 1L] <- 1L
          N_Ends[N_Counts + 1L] <- B_Ends[C_Counts] - F_Break
          RemainingBounds <- FALSE
        } else if (F_Break >= B_Ends[C_Counts] &
                   C_Counts < length(B_Ends)) {
          # print("break after, further evaluations exist")
          # Case 4:
          # Break occurs after this bound set
          # Bound sets remain to be evaluated
          # add current set whole cloth
          # update counters
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- B_Ends[C_Counts]
          N_Counts <- N_Counts + 1L
          C_Counts <- C_Counts + 1L
        } else if (F_Break < B_Starts[C_Counts] &
                   C_Counts < length(B_Starts)) {
          # print("break prior, further evaluations exist")
          # Case 5:
          # Break occurs before this bound set
          # Bound sets remain to be evaluated
          # Adjust all following bound sets, but not current
          # update counters
          # remaining sets will be visited in future
          B_Starts[(C_Counts + 1L):length(B_Starts)] <- B_Starts[(C_Counts + 1L):length(B_Starts)] - F_Break
          B_Ends[(C_Counts + 1L):length(B_Ends)] <- B_Ends[(C_Counts + 1L):length(B_Ends)] - F_Break
          N_Starts[N_Counts] <- B_Starts[C_Counts] - F_Break
          N_Ends[N_Counts] <- B_Ends[C_Counts] - F_Break
          N_Counts <- N_Counts + 1L
          C_Counts <- C_Counts + 1L
        } else if (F_Break >= B_Starts[C_Counts] &
                   F_Break < B_Ends[C_Counts]) {
          # print("break interior, further evaluations exist")
          # Case 6:
          # Break occurs inside current bound set
          # Bound sets remain to be evaluated
          # Adjust following bound sets
          # add set with adjusted current end
          # add new set
          # update following bounds
          N_Starts[N_Counts] <- B_Starts[C_Counts]
          N_Ends[N_Counts] <- F_Break
          N_Starts[N_Counts + 1L] <- 1L
          N_Ends[N_Counts + 1L] <- B_Ends[C_Counts] - F_Break # check off by 1!
          B_Starts[(C_Counts + 1L):length(B_Starts)] <- B_Starts[(C_Counts + 1L):length(B_Starts)] - F_Break
          B_Ends[(C_Counts + 1L):length(B_Ends)] <- B_Ends[(C_Counts + 1L):length(B_Ends)] - F_Break
          C_Counts <- C_Counts + 1L
          N_Counts <- N_Counts + 2L
        } else {
          print("an unmet condition exists")
          stop("An unmet condition was discovered while adjusting bounds, please contact maintainer!")
        }
      } # while loop
      N_Starts <- N_Starts[!is.na(N_Starts)]
      N_Ends <- N_Ends[!is.na(N_Ends)]
      
      res$Range[[a1]] <- IRanges(start = N_Starts,
                                 end = N_Ends)
    } # conditional to check
    
    if (verbose) {
      setTxtProgressBar(pb = pBar,
                        value = a1 / PBAR)
    }
  } # end a1 loop
  
  # check phases now:
  # offset vals
  # 0 == in phase with reading frame started by the first subfeature
  # i.e. starting at codon position 1
  # 1 == subfeature is starting at codon position 2
  # 2 == subfeature is starting at codon position 3
  tmpvec <- lapply(X = res$Range,
                   FUN = function(x) {
                     if (length(x) == 1) {
                       0
                     } else {
                       c(0, cumsum(width(x[-length(x)])) %% 3)
                     }
                   })
  tmpvec <- IntegerList(tmpvec)
  res$Phase <- tmpvec
  
  if (verbose) {
    close(pBar)
    tend <- Sys.time()
    print(tend - tstart)
  }
  attr(x = res,
       which = "heirarchies") <- UClusts
  # class(res) <- c("data.frame",
  #                 "genecalls")
  return(res)
}




