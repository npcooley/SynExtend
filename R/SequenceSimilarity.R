###### -- Sequence Similarity Score -------------------------------------------
# Written by: Erik Wright
# Contact: ESWright@pitt.edu
# Maintained by: Nicholas Cooley
# Contact: npc19@pitt.edu

# return a numeric between 1 and a large negative value that represents
# the similarity of two aligned sequences given a substitution matrix

SequenceSimilarity <- function(Seqs,
                               SubMat,
                               penalizeGapLetter = TRUE,
                               includeTerminalGaps = TRUE,
                               allowNegative = TRUE) {
  
  # overhead checking
  # Seqs must be of a stringset of length 2
  seqclass <- class(Seqs)
  
  if (!(seqclass %in% c("DNAStringSet", "AAStringSet"))) {
    stop ("Seqs must be a an object of class DNAStringSet or AAStringSet.")
  }
  if (length(Seqs) != 2L) {
    stop ("Function is only designed to work with a single pairwise alignment.")
  }
  # if SubMat argument is not specified, create a generic substitution matrix
  if (missing(SubMat)) {
    if (is(object = Seqs,
           class2 = "DNAStringSet")) {
      # build a generic DNA Substitution matrix where all substitutions are equal
      SubMat <- diag(length(DNA_ALPHABET))
      dimnames(SubMat) <- list(DNA_ALPHABET,
                               DNA_ALPHABET)
    } else if (is(object = Seqs,
                  class2 = "AAStringSet")) {
      # use the 40th PFASUM matrix
      # PFASUM data object loads an object named `PFASUM` but this form ensures
      # that there is a visible binding for check and bioccheck
      PFASUM <- get(data(list = "PFASUM",
                         envir = environment(),
                         package = "DECIPHER"))
      SubMat <- PFASUM[1:20, 1:20, 40]
    }
  }
  
  # ensure that the substition matrix can support all given characters
  PresentCharacters <- alphabetFrequency(x = Seqs)
  PresentCharacters <- names(PresentCharacters[colSums(PresentCharacters) > 0])
  # remove alignment characters
  PresentCharacters <- PresentCharacters[!(PresentCharacters %in% c("-", "+", "."))]
  SubMatChars <- colnames(SubMat)
  # all of the present characters must be accounted for in the substitution
  # matrix
  # return(list(PresentCharacters,
  #             SubMatChars))
  if (!all(PresentCharacters %in% SubMatChars)) {
    stop ("Substitution matrix does not appear to contain all required characters.")
  }
  
  # include terminal gaps, or not i.e. global vs local scoring
  if (!includeTerminalGaps) {
    t <- TerminalChar(Seqs)
    p1 <- max(t[, 1]) + 1L
    p2 <- p1 + min(t[, 3]) - 1L
    Seqs <- subseq(Seqs, p1, p2)
  }
  Seqs <- unname(as.matrix(Seqs))
  
  m1 <- match(Seqs[1,], rownames(SubMat))
  m2 <- match(Seqs[2,], rownames(SubMat))
  w <- is.na(m1) | is.na(m2)
  m1 <- m1[!w]
  m2 <- m2[!w]
  
  s11 <- sum(diag(SubMat) * tabulate(m1, nrow(SubMat)))
  s22 <- sum(diag(SubMat) * tabulate(m2, nrow(SubMat)))
  
  if (s11 > 1e-8 & s22 > 1e-8) {
    # grab the substitution values from the matrix
    # based on the vectors that match alignmed positions
    s12 <- sum(SubMat[cbind(m1, m2)], na.rm = TRUE)
    if (allowNegative) {
      s <- s12/sqrt(s11*s22)
      if (penalizeGapLetter)
        s <- s - sum(w)/length(m1)
    } else {
      if (s12 > 0) {
        if (penalizeGapLetter) {
          s <- s12/sqrt(s11*s22) - sum(w)/length(m1)
          if (s < 0)
            s <- 0
        } else {
          s <- s12/sqrt(s11*s22)
        }
      } else {
        s <- 0
      }
    }
  } else {
    s <- 0
  }
  return(s)
}

