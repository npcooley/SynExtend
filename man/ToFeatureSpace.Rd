\name{ToFeatureSpace}
\alias{ToFeatureSpace}
\title{
Tabulating Features Linked by Syntenic Hits
}
\description{
A function for concisely tabulating where genomic features are connected by syntenic hits.
}
\usage{
ToFeatureSpace(hit_blocks,
               qranges,
               qstrands,
               sranges,
               sstrands)
}
\arguments{
  \item{hit_blocks}{
    A list of \code{data.frame}s created by calling \code{split} on \code{as.data.frame} of the lower triangle of the \code{\link{NucleotideOverlap}} object. Splitting on the interaction of the \code{QueryGene} and \code{SubjectGene} columns.
  }
  \item{qranges}{
    An \code{IRangesList}.
  }
  \item{qstrands}{
    A vector of integers, \code{0} for a feature in the negative strand, \code{1} for the positive strand.
  }
  \item{sranges}{
    An \code{IRangesList}.
  }
  \item{sstrands}{
    A vector of integers, \code{0} for a feature in the negative strand, \code{1} for the positive strand.
  }
}
\details{
  Builds a list of 4-row x n-column matrices from the information provided by \code{\link{NucleotideOverlap}}.
}
\value{
  A list of matrices to be passed to \code{\link{AlignPairs}} as anchors.
}
\author{
  Nicholas Cooley \email{Nicholas.Cooley@ul.ie}
}

\seealso{
  \itemize{
    \item \code{\link{FindSynteny}}
    \item \code{\link{NucleotideOverlap}}
    \item \code{\link{SummarizePairs}}
    \item \code{\link{AlignPairs}}
  }
}
\examples{
data("genecalls", package = "SynExtend")
data("syn", package = "SynExtend")

Links <- NucleotideOverlap(SyntenyObject = syn,
                           GeneCalls = genecalls,
                           LimitIndex = FALSE,
                           Verbose = TRUE)
}
