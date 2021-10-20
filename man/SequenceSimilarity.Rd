\name{SequenceSimilarity}
\alias{SequenceSimilarity}
\title{
Return a numeric value that represents the similarity between two aligned sequences as determined by a provided subsitution matrix.
}
\description{
Takes in a \code{DNAStringSet} or \code{AAStringSet} representing a pairwise alignment and a subsitution matrix such as those present in \code{PFASUM}, and return a numeric value representing sequence similarity as defined by the substitution matrix.
}
\usage{
SequenceSimilarity(Seqs,
                   SubMat,
                   penalizeGapLetter = TRUE,
                   includeTerminalGaps = TRUE,
                   allowNegative = TRUE)
}
\arguments{
  \item{Seqs}{
A \code{DNAStringSet} or \code{AAStringSet} of length 2.
}
  \item{SubMat}{
A named matrix representing a substitution matrix. If left ``NULL'' and ``Seqs'' is a \code{AAStringSet}, the 40th ``PFASUM'' matrix is used. If left ``NULL'' and ``Seqs'' is a \code{DNAStringSet}, a matrix with only the diagonal filled with ``1'''s is used.
}
  \item{penalizeGapLetter}{
A logical indicating whether or not to penalize Gap-Letter matches. Defaults to ``TRUE''.
}
  \item{includeTerminalGaps}{
A logical indicating whether or not to penalize terminal matches. Defaults to ``TRUE''.
}
  \item{allowNegative}{
A logical indicating whether or not allow negative scores. Defaults to ``TRUE''. If ``FALSE'' scores that are returned as less than zero are converted to zero.
}
}
\details{
Takes in a \code{DNAStringSet} or \code{AAStringSet} representing a pairwise alignment and a subsitution matrix such as those present in \code{PFASUM}, and return a numeric value representing sequence similarity as defined by the substitution matrix.
}
\value{
Returns a single numeric.
}
\author{
Erik Wright \email{ESWRIGHT@pitt.edu}
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{AlignSeqs}}, \code{\link{AlignProfiles}}, \code{\link{AlignTranslation}}, \code{\link{DistanceMatrix}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, remove="all")
alignedDNA <- AlignSeqs(dna[1:2])

DNAPlaceholder <- diag(15)
dimnames(DNAPlaceholder) <- list(DNA_ALPHABET[1:15],
                                 DNA_ALPHABET[1:15])

SequenceSimilarity(Seqs = alignedDNA,
                   SubMat = DNAPlaceholder,
                   includeTerminalGaps = TRUE,
                   penalizeGapLetter = TRUE,
                   allowNegative = TRUE)
}
