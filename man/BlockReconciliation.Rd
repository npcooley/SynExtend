\name{BlockReconciliation}
\alias{BlockReconciliation}
\title{
Rejection scheme for asyntenic predicted pairs
}
\description{
Take in a \code{PairSummaries} object and reject predicted pairs that conflict with syntenic blocks either locally or globally.
}
\usage{
BlockReconciliation(Pairs,
                    ConservativeRejection = TRUE,
                    AlwaysReject = TRUE,
                    Verbose = FALSE)
}
\arguments{
  \item{Pairs}{
A \code{PairSummaries} object.
}
  \item{ConservativeRejection}{
A logical defaulting to \code{TRUE}. By default only pairs that conflict within a syntenic block will be rejected. When \code{FALSE} any conflict with a cause the rejection of the pair in the smaller block.
}
  \item{AlwaysReject}{
Currently unimplemented. For future implementations of more complex rejection schemes.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
If a given \code{PairSummaries} object contains predicted pairs that conflict, i.e. imply paralogy, or an ``incorrect'' and a ``correct'' ortholog prediction, these predictions will be reconciled. The function scrolls through pairs based on the size of the syntenic block that they are part of, from largest to smallest. When \code{ConservativeRejection} is \code{TRUE} only predicted pairs that exist within the syntenic block ``space'' will be removed, this option leaves room for conflicting predictions to remain if they are non-local to each other, or are on different indices. When \code{ConservativeRejection} is \code{FALSE} any pair that conflicts with a larger syntenic block will be rejected. This option forces only 1-1 feature pairings, for features are part of any syntenic block. Predicted pairs that represent a syntenic block size of 1 feature will not reject other pairs.
} 
\value{
A data.frame of class ``data.frame'' and ``PairSummaries'' of paired genes that are connected by syntenic hits. Contains columns describing the k-mers that link the pair. Columns ``p1'' and ``p2'' give the location ids of the the genes in the pair in the form ``DatabaseIdentifier_ContigIdentifier_GeneIdentifier''. ``ExactMatch'' provides an integer representing the exact number of nucleotides contained in the linking k-mers. ``TotalKmers'' provides an integer describing the number of distinct k-mers linking the pair. ``MaxKmer'' provides an integer describing the largest k-mer that links the pair. A column titled ``Consensus'' provides a value between zero and 1 indicating whether the kmers that link a pair of features are in the same position in each feature, with 1 indicating they are in exactly the same position and 0 indicating they are in as different a position as is possible. The ``Adjacent'' column provides an integer value ranging between 0 and 2 denoting whether a feature pair's direct neighbors are also paired. Gap filled pairs neither have neighbors, or are included as neighbors. The ``TetDist'' column provides the euclidean distance between oligonucleotide - of size 4 - frequences between predicted pairs. ``PIDType'' provides a character vector with values of ``NT'' where either of the pair indicates it is not a translatable sequence or ``AA'' where both sequences are translatable. If users choose to perform pairwise alignments there will be a ``PID'' column providing a numeric describing the percent identity between the two sequences. If users choose to predict PIDs using their own, or a provided model, a ``PredictedPID'' column will be provided.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}, \code{\link{PairSummaries}}
}
\examples{
data("Nitrosocosmicus_Pairs02", package = "SynExtend")
Pairs03 <- BlockReconciliation(Pairs = Nitrosocosmicus_Pairs02,
                               ConservativeRejection = FALSE,
                               Verbose = TRUE)
}
