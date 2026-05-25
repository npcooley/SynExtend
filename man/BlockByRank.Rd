\name{BlockByRank}
\alias{BlockByRank}
\title{
Return simple summaries of blocks of candidate pairs.
}
\description{
This function is designed to work internally to \code{\link{SummarizePairs}} so it works on relatively simple atomic vectors and has little overhead checking. All arguments must be the same length.
}
\usage{
BlockByRank(index1,
            partner1,
            index2,
            partner2)
}
\arguments{
  \item{index1}{
  Integer; references the contigs containing candidate feature partners.
}
  \item{partner1}{
  Integer; references the candidate feature partners by row position in the source DataFrame.
}
  \item{index2}{
  Integer; references the contigs containing candidate feature partners.
}
  \item{partner2}{
  Integer; references the candidate feature partners by row position in the source DataFrame.
}
}
\details{
  \code{BlockByRank} uses the diagonal rank to identify where runs of candidate features are present in sequential blocks. In cases where a candidate feature is part of two competing blocks it is assigned to the larger.
}
\value{
A list with named elements \code{absblocksize} and \code{blockidmap}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}, \code{\link{FindSynteny}}
}
\examples{
data("init_pairs", package = "SynExtend")
x <- paste(init_pairs$p1, init_pairs$p2, sep = "_")
x <- do.call(rbind, strsplit(x = x, split = "_", fixed = TRUE))
x <- matrix(data = as.integer(x), nrow = nrow(x))
y <- BlockByRank(index1 = x[, 2],
                 partner1 = x[, 3],
                 index2 = x[, 5],
                 partner2 = x[, 6])
}
