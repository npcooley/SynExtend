###### -- get translated features ---------------------------------------------
# author: npc
# maintainer: Nicholas.cooley@ul.ie

###### -- NOTES ---------------------------------------------------------------
# lightweight for now ... what kind of checks do we need here ?

###### -- FUNCTION ------------------------------------------------------------
# take in a DNAStringSet
# a DataFrame of gene calls with at least the columns:
# Coding
# Translation_Table
# Range

# we *should* inherit the names from the names from the DNAStringSet
# and we need to think about the default translation table

GetTranslatedFeatures <- function(Nucs,
                                  GeneCalls,
                                  DefaultTranslationTable = "11") {
  
  # translate where we can
  # where we are indicated to be coding
  # where we fit into a frame
  # figure out if we need to drop in the default translation table
  tr_tbl1 <- GeneCalls$Coding
  tr_tbl2 <- GeneCalls$Translation_Table
  tr_tbl3 <- lapply(X = GeneCalls$Range,
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
    aa_ph[[a3]] <- translate(x = Nucs[w1 & w2],
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
    # order from Nucs should be retained
    aa_seqs <- aa_ph[[1]]
  }
  return(aa_seqs)
  
}
