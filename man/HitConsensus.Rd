\name{HitConsensus}
\alias{HitConsensus}
\title{
Return a numeric measure of whether kmer hits linking two genomic features are in linearly similar locations in both features.
}
\description{
This function is designed to work internally to \code{\link{SummarizePairs}} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
HitConsensus(gene1left,
             gene2left,
             gene1right,
             gene2right,
             strand1,
             strand2,
             hit1left,
             hit1right,
             hit2left,
             hit2right)
}
\arguments{
  \item{gene1left}{
  Integer; feature bound positions in nucleotide space.
}
  \item{gene2left}{
  Integer; feature bound positions in nucleotide space.
}
  \item{gene1right}{
  Integer; feature bound positions in nucleotide space.
}
  \item{gene2right}{
  Integer; feature bound positions in nucleotide space.
}
  \item{strand1}{
  Logical; is feature 1 on the positive or negative strand
}
  \item{strand2}{
  Logical; is feature 2 on the positive or negative strand
}
  \item{hit1left}{
  Integer; kmer hit bound positions in nucleotide space.
}
  \item{hit1right}{
  Integer; kmer hit bound positions in nucleotide space.
}
  \item{hit2left}{
  Integer; kmer hit bound positions in nucleotide space.
}
  \item{hit2right}{
  Integer; kmer hit bound positions in nucleotide space.
}
}
\details{
  \code{HitConsensus} calculates whether the distances between the bounds of a kmer hit and the feature bounds are different between the features linked by the kmer.
}
\value{
  A vector of numerics.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}, \code{\link{FindSynteny}}
}
\examples{
# no current example...
}
