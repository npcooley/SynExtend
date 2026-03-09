###### -- extract features from DF object -------------------------------------

###### -- NOTES ---------------------------------------------------------------
# low overhead and minimal error checking for now

###### -- FUNCTION ------------------------------------------------------------

FeaturesFromDF <- function(Genome,
                           GeneCalls) {
  if (!is(object = Genome,
          class2 = "DNAStringSet")) {
    stop ("Genome must be a DNAStringSet")
  }
  if (!is(object = GeneCalls,
          class2 = "DataFrame")) {
    stop ("GeneCalls must be a DataFrame")
  }
  # we expect index to be appropriately matched to the fna order now, but
  # non-NCBI/EMBL sources may still cause problems there
  check_cols <- c("Index",
                  "Range",
                  "Strand")
  if (!all(check_cols %in% colnames(GeneCalls))) {
    stop ("GeneCalls columns must contain the columns: ",
          paste(check_cols,
                collapse = ", "))
  }
  # not yet implemented
  # NameBy <- match.arg(NameBy)
  u_indices <- unique(GeneCalls$Index)
  res <- vector(mode = "list",
                length = length(u_indices))
  
  for (a1 in seq_along(u_indices)) {
    curr_index <- GeneCalls$Index == u_indices[a1]
    curr_ranges <- unname(GeneCalls$Range[curr_index])
    # logical, reverseComplement when TRUE
    curr_strand <- GeneCalls$Strand[curr_index] == 1L
    
    curr_lengths <- lengths(curr_ranges)
    curr_ranges <- unlist(curr_ranges,
                          recursive = FALSE)
    curr_feats <- extractAt(x = Genome[[u_indices[a1]]],
                            at = curr_ranges)
    feat_map <- rep(x = seq_along(curr_lengths),
                    times = curr_lengths)
    
    split_feats <- split(x = curr_feats,
                         f = feat_map)
    # return(split_feats)
    split_feats <- lapply(X = split_feats,
                          FUN = function(x) {
                            if (length(x) == 1) {
                              x[[1]]
                            } else {
                              unlist(x)
                            }
                          })
    curr_feats <- DNAStringSet(split_feats)
    if (any(curr_strand)) {
      curr_feats[curr_strand] <- reverseComplement(curr_feats[curr_strand])
    }
    res[[a1]] <- curr_feats
  }
  # if a1 was a loop of length 1, do.call is not needed
  if (a1 > 1L) {
    res <- do.call(c,
                   res)
  } else {
    res <- curr_feats
  }
  
  return(res)
  
}
